#!/usr/bin/env python
from __future__ import annotations

import argparse
import pickle
import sys
from pathlib import Path

import numpy as np

from argtest_common import (
    accessible_intervals_from_mu,
    infer_mu_base,
    infer_mu_path,
    load_ts,
    overlap_lengths,
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
