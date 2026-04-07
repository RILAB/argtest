#!/usr/bin/env python
from __future__ import annotations

import argparse
from pathlib import Path

from argtest_common import collapse_masked_intervals, dump_ts, load_ts, ratemap_from_keep_intervals
from trim_regions import complement_intervals, load_mask_intervals


def parse_args():
    parser = argparse.ArgumentParser(
        description="Trim one tree sequence using one BED mask and write one output tree sequence."
    )
    parser.add_argument("--ts", required=True, type=Path, help="Input tree sequence file.")
    parser.add_argument("--remove", required=True, type=Path, help="BED file with intervals to remove.")
    parser.add_argument("--out", required=True, type=Path, help="Output tree sequence file.")
    parser.add_argument(
        "--log",
        type=Path,
        default=None,
        help="Optional log file path (default: <out.parent>/logs/<ts_stem>.trim_regions.log).",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    ts = load_ts(args.ts)
    log_path = args.log or (args.out.parent / "logs" / f"{args.ts.stem}.trim_regions.log")
    log_path.parent.mkdir(parents=True, exist_ok=True)
    masked = load_mask_intervals(args.remove, ts.sequence_length)
    keep = complement_intervals(masked, ts.sequence_length)
    accessible = ratemap_from_keep_intervals(keep, ts.sequence_length)
    trimmed = collapse_masked_intervals(ts, accessible)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    dump_ts(trimmed, args.out)
    with open(log_path, "w") as fh:
        fh.write("# trim_regions_single summary\n")
        fh.write(f"ts={args.ts}\n")
        fh.write(f"remove={args.remove}\n")
        fh.write(f"out={args.out}\n")
        fh.write(f"sequence_length={ts.sequence_length}\n")
    print(f"Wrote trimmed tree: {args.out}")
    print(f"Wrote log: {log_path}")


if __name__ == "__main__":
    main()
