#!/usr/bin/env python
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from argtest_common import collapse_masked_intervals, dump_ts, load_ts, merge_intervals, ratemap_from_keep_intervals


def parse_args():
    p = argparse.ArgumentParser(
        description="Trim regions specified in a BED file from a directory of tree sequences."
    )
    p.add_argument(
        "--ts-dir",
        required=True,
        type=Path,
        help="Directory containing tree sequence files (.tsz, .ts, .trees).",
    )
    p.add_argument(
        "--remove",
        required=True,
        type=Path,
        help="BED file of regions to remove.",
    )
    p.add_argument(
        "--out-dir",
        type=Path,
        default=None,
        help="Output directory for trimmed tree files (default: <ts-dir>/trimmed).",
    )
    p.add_argument(
        "--pattern",
        default="*",
        help="Optional glob pattern to filter tree sequence filenames (default: '*').",
    )
    p.add_argument(
        "--log",
        type=Path,
        default=None,
        help="Optional log file path (default: <out-dir>/logs/trim_regions.log).",
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


def load_mask_intervals(path: Path, sequence_length: float):
    intervals = []
    with open(path, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 3:
                raise ValueError(f"Invalid BED line in {path}: {line}")
            start = max(0.0, float(parts[1]))
            end = min(float(sequence_length), float(parts[2]))
            if end > start:
                intervals.append([start, end])
    return merge_intervals(intervals)


def complement_intervals(intervals, sequence_length: float):
    keep = []
    last = 0.0
    for left, right in intervals:
        if left > last:
            keep.append([last, left])
        last = max(last, right)
    if last < float(sequence_length):
        keep.append([last, float(sequence_length)])
    return keep


def output_name(ts_path: Path, bed_path: Path) -> str:
    return f"{ts_path.stem}.trimmed.{bed_path.stem}{ts_path.suffix}"


def main():
    args = parse_args()
    ts_files = find_tree_files(args.ts_dir, args.pattern)
    out_dir = args.out_dir or (args.ts_dir / "trimmed")
    out_dir.mkdir(parents=True, exist_ok=True)
    log_path = args.log or (out_dir / "logs" / "trim_regions.log")
    log_path.parent.mkdir(parents=True, exist_ok=True)

    first_ts = load_ts(ts_files[0])
    masked = load_mask_intervals(args.remove, first_ts.sequence_length)
    keep = complement_intervals(masked, first_ts.sequence_length)
    accessible = ratemap_from_keep_intervals(keep, first_ts.sequence_length)

    with open(log_path, "w") as log:
        log.write("# trim_regions\n")
        log.write(f"# ts_dir={args.ts_dir}\n")
        log.write(f"# out_dir={out_dir}\n")
        log.write(f"# remove_bed={args.remove}\n\n")
        log.write(f"reference_ts={ts_files[0]}\n")
        log.write(f"sequence_length={first_ts.sequence_length}\n")
        log.write(f"n_masked_intervals={len(masked)}\n")
        log.write(f"masked_bp={sum(right - left for left, right in masked)}\n\n")

        for ts_file in ts_files:
            ts = load_ts(ts_file)
            if float(ts.sequence_length) != float(first_ts.sequence_length):
                raise RuntimeError(
                    f"Sequence length mismatch for {ts_file}: "
                    f"{ts.sequence_length} != {first_ts.sequence_length}"
                )
            ts2 = collapse_masked_intervals(ts, accessible)
            out_path = out_dir / output_name(ts_file, args.remove)
            dump_ts(ts2, out_path)

            log.write("=" * 80 + "\n")
            log.write(f"ts_file={ts_file}\n")
            log.write(f"out_file={out_path}\n")
            log.write(f"old_L={ts.sequence_length}\n")
            log.write(f"new_L={ts2.sequence_length}\n\n")
            print(f"Wrote: {out_path.name}")

    # Summary printed to stdout and stderr and appended to log
    summary = f"Processed {len(ts_files)} tree files; masked_intervals={len(masked)}; masked_bp={sum(right - left for left, right in masked)}"
    print(summary)
    print(summary, file=sys.stderr)
    print(f"Done. Log: {log_path}")
    try:
        with open(log_path, "a") as log:
            log.write("# Summary\n")
            log.write(summary + "\n")
    except Exception:
        print(f"WARNING: failed to append summary to log {log_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
