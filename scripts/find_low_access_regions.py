#!/usr/bin/env python
from __future__ import annotations

import argparse
import pickle
import re
import sys
from pathlib import Path

import numpy as np

from argtest_common import accessible_intervals_from_mu, load_ts, overlap_lengths


def infer_mu_base(ts_stem: str) -> list[str]:
    bases = [ts_stem]
    m = re.match(r"^(.+)\.(\d+)$", ts_stem)
    if m:
        bases.append(m.group(1))
    m = re.match(r"^(.+)[_-](\d+)$", ts_stem)
    if m:
        bases.append(m.group(1))
    dedup = []
    seen = set()
    for b in bases:
        if b not in seen:
            dedup.append(b)
            seen.add(b)
    return dedup


def infer_mu_path(ts_path: Path) -> Path:
    bases = infer_mu_base(ts_path.stem)
    search_dirs = [ts_path.parent, ts_path.parent.parent]
    for d in search_dirs:
        for b in bases:
            p = d / f"{b}.mut_rate.p"
            if p.exists():
                return p
    candidates = []
    for d in search_dirs:
        if not d.exists():
            continue
        for p in d.glob("*.mut_rate.p"):
            nm = p.name
            if any(nm.startswith(f"{b}.") or nm == f"{b}.mut_rate.p" for b in bases):
                candidates.append(p)
    candidates = sorted(set(candidates))
    if len(candidates) == 1:
        return candidates[0]
    if len(candidates) > 1:
        raise RuntimeError(
            f"Ambiguous mutation maps for {ts_path.name}: "
            + ", ".join(str(x) for x in candidates)
        )
    raise FileNotFoundError(
        f"Could not infer mutation map for {ts_path}. Tried bases={bases} in {search_dirs}"
    )


def format_num(x: float) -> str:
    if float(x).is_integer():
        return str(int(x))
    return str(x).replace(".", "p")


def parse_args():
    p = argparse.ArgumentParser(
        description="Identify windows with too few accessible bp and write them as a BED file."
    )
    p.add_argument(
        "--ts-dir",
        required=True,
        type=Path,
        help="Directory containing tree sequence files (.tsz, .ts, .trees).",
    )
    p.add_argument(
        "--window-size",
        required=True,
        type=float,
        help="Window size in bp for accessibility checks.",
    )
    p.add_argument(
        "--cutoff-bp",
        required=True,
        type=float,
        help="Minimum accessible bp in a window to keep it.",
    )
    p.add_argument(
        "--pattern",
        default="*",
        help="Optional glob pattern to filter tree sequence filenames (default: '*').",
    )
    p.add_argument(
        "--out",
        type=Path,
        default=None,
        help="Output BED path (default: <ts-dir>/low_access.ws<window>.accbp<cutoff>.bed).",
    )
    p.add_argument(
        "--log",
        type=Path,
        default=None,
        help="Optional log file path (default: same base as output with .log).",
    )
    return p.parse_args()


def find_tree_files(ts_dir: Path, pattern: str) -> list[Path]:
    if not ts_dir.exists():
        raise FileNotFoundError(f"Tree directory does not exist: {ts_dir}")
    files = sorted(
        [
            p
            for p in ts_dir.glob(pattern)
            if p.is_file() and p.suffix in {".tsz", ".ts", ".trees"}
        ]
    )
    if not files:
        raise RuntimeError(f"No tree files found in {ts_dir} matching pattern '{pattern}'.")
    return files


def default_out_path(ts_dir: Path, window_size: float, cutoff_bp: float) -> Path:
    return ts_dir / f"low_access.ws{format_num(window_size)}.accbp{format_num(cutoff_bp)}.bed"


def main():
    args = parse_args()
    ts_files = find_tree_files(args.ts_dir, args.pattern)
    first_ts = load_ts(ts_files[0])
    mu_path = infer_mu_path(ts_files[0])
    with open(mu_path, "rb") as fh:
        mu = pickle.load(fh)
    sequence_length = float(first_ts.sequence_length)
    windows = np.arange(0, sequence_length + args.window_size, args.window_size, dtype=float)
    if windows[-1] > sequence_length:
        windows[-1] = sequence_length
    acc_bp = overlap_lengths(accessible_intervals_from_mu(mu), windows)

    out_path = args.out or default_out_path(args.ts_dir, args.window_size, args.cutoff_bp)
    lines = []
    chrom = ts_files[0].stem
    for i in range(len(windows) - 1):
        if acc_bp[i] >= args.cutoff_bp:
            continue
        lines.append(
            f"{chrom}\t{int(windows[i])}\t{int(windows[i + 1])}\t{float(acc_bp[i]):.3f}"
        )
    out_path.write_text("\n".join(lines) + ("\n" if lines else ""))
    n_bad = int(np.sum(acc_bp < args.cutoff_bp))
    print(f"Wrote: {out_path}")
    summary = f"Flagged {n_bad} low-accessibility windows out of {len(windows) - 1}"
    print(summary)
    print(summary, file=sys.stderr)

    # Write a simple log
    log_path = args.log or out_path.with_suffix(".log")
    try:
        with open(log_path, "w") as fh:
            fh.write("# find_low_access_regions summary\n")
            fh.write(f"out_path={out_path}\n")
            fh.write(summary + "\n")
    except Exception:
        print(f"WARNING: failed to write log to {log_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
