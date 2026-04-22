#!/usr/bin/env python
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from argtest_common import merge_intervals


def parse_args():
    parser = argparse.ArgumentParser(
        description="Combine one or more BED masks, merge overlapping intervals, and write BED output."
    )
    parser.add_argument("--chrom", required=True, help="Chromosome name written in BED column 1.")
    parser.add_argument("--out", required=True, type=Path, help="Output BED file.")
    parser.add_argument(
        "--inputs",
        nargs="+",
        required=True,
        type=Path,
        help="Input BED files. Missing files are ignored.",
    )
    parser.add_argument(
        "--log",
        type=Path,
        default=None,
        help="Optional log file path (default: <out.parent>/logs/<chrom>.combine_masks.log).",
    )
    return parser.parse_args()


def read_intervals(path: Path):
    intervals = []
    if not path.exists():
        print(f"WARNING: input BED not found, skipping: {path}", file=sys.stderr)
        return intervals
    with open(path, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            start = float(parts[1])
            end = float(parts[2])
            if end > start:
                intervals.append([start, end])
    return intervals


def main():
    args = parse_args()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    log_path = args.log or (args.out.parent / "logs" / f"{args.chrom}.combine_masks.log")
    log_path.parent.mkdir(parents=True, exist_ok=True)
    merged = merge_intervals(
        [
            interval
            for path in args.inputs
            for interval in read_intervals(path)
        ]
    )
    lines = [f"{args.chrom}\t{int(left)}\t{int(right)}" for left, right in merged]
    args.out.write_text("\n".join(lines) + ("\n" if lines else ""))
    with open(log_path, "w") as fh:
        fh.write("# combine_remove_masks summary\n")
        fh.write(f"chrom={args.chrom}\n")
        fh.write(f"out={args.out}\n")
        fh.write(f"inputs={len(args.inputs)}\n")
        fh.write(f"merged_intervals={len(lines)}\n")
    print(f"Wrote merged BED: {args.out}")
    print(f"Wrote log: {log_path}")


if __name__ == "__main__":
    main()
