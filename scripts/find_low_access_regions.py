#!/usr/bin/env python
from __future__ import annotations

import argparse
import pickle
import re
import sys
from pathlib import Path

import numpy as np

from argtest_common import accessible_intervals_from_mu, load_ts, overlap_lengths


def infer_mu_base(ts_stem: str, parent_stem: str | None = None) -> list[str]:
    bases = [ts_stem]
    if parent_stem:
        bases.append(parent_stem)
    m = re.match(r"^(.+)\.(\d+)$", ts_stem)
    if m:
        bases.append(m.group(1))
    m = re.match(r"^(.+)[_-](\d+)$", ts_stem)
    if m:
        bases.append(m.group(1))
    if parent_stem:
        m = re.match(r"^(.+)\.(\d+)$", parent_stem)
        if m:
            bases.append(m.group(1))
        m = re.match(r"^(.+)[_-](\d+)$", parent_stem)
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
    bases = infer_mu_base(ts_path.stem, parent_stem=ts_path.parent.name)
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
        "--ts",
        required=True,
        type=Path,
        help="Tree sequence file (.tsz, .ts, .trees) used to infer sequence length and mutation map.",
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


def default_out_path(ts_path: Path, window_size: float, cutoff_bp: float) -> Path:
    return ts_path.parent / f"low_access.ws{format_num(window_size)}.accbp{format_num(cutoff_bp)}.bed"


def default_log_path(out_path: Path) -> Path:
    return out_path.parent / "logs" / f"{out_path.stem}.log"


def main():
    args = parse_args()
    ts_path = args.ts
    first_ts = load_ts(ts_path)
    mu_path = infer_mu_path(ts_path)
    with open(mu_path, "rb") as fh:
        mu = pickle.load(fh)
    sequence_length = float(first_ts.sequence_length)
    windows = np.arange(0, sequence_length + args.window_size, args.window_size, dtype=float)
    if windows[-1] > sequence_length:
        windows[-1] = sequence_length
    acc_bp = overlap_lengths(accessible_intervals_from_mu(mu), windows)

    out_path = args.out or default_out_path(ts_path, args.window_size, args.cutoff_bp)
    lines = []
    chrom = ts_path.stem
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
    log_path = getattr(args, "log", None) or default_log_path(out_path)
    try:
        log_path.parent.mkdir(parents=True, exist_ok=True)
        with open(log_path, "w") as fh:
            fh.write("# find_low_access_regions summary\n")
            fh.write(f"out_path={out_path}\n")
            fh.write(summary + "\n")
    except Exception:
        print(f"WARNING: failed to write log to {log_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
