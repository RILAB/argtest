#!/usr/bin/env python
from __future__ import annotations

import argparse
from pathlib import Path

from mutload_common import (
    assert_sample_ids_preserved,
    dump_ts,
    load_remove_intervals,
    load_ts,
    name_to_nodes_map,
    trim_ts_by_intervals,
    validate_trimmed_ts,
)


def parse_remove_list(values):
    # Accept comma-separated lists or repeated flags.
    if not values:
        return []
    paths = []
    for value in values:
        parts = [v.strip() for v in value.split(",")]
        paths.extend([p for p in parts if p])
    return paths


def parse_args():
    p = argparse.ArgumentParser(
        description="Remove individuals over BED intervals and write a trimmed tree sequence",
    )
    p.add_argument("ts", help="Tree sequence file (.ts, .trees, or .tsz)")
    p.add_argument(
        "--remove",
        action="append",
        required=True,
        help="BED file(s) of regions to remove per individual (comma-separated or repeated)",
    )
    p.add_argument("--out", help="Output tree sequence path (.ts, .trees, or .tsz)")
    p.add_argument("--suffix-to-strip", default="_anchorwave")
    return p.parse_args()


def main():
    args = parse_args()
    ts_path = Path(args.ts)
    ts = load_ts(ts_path)

    remove_paths = parse_remove_list(args.remove)
    remove_intervals = load_remove_intervals(remove_paths)

    name_to_nodes = name_to_nodes_map(ts, suffix_to_strip=args.suffix_to_strip)
    # Drop edges/mutations for removed individuals across their intervals.
    trimmed_ts = trim_ts_by_intervals(ts, remove_intervals, name_to_nodes)
    assert_sample_ids_preserved(ts, trimmed_ts)
    validate_trimmed_ts(trimmed_ts)

    if args.out:
        out_path = Path(args.out)
    else:
        # Default output to results/ with a trimmed suffix.
        out_dir = Path(__file__).resolve().parent.parent / "results"
        out_dir.mkdir(parents=True, exist_ok=True)
        out_path = out_dir / f"{ts_path.stem}_trimmed.tsz"
    dump_ts(trimmed_ts, out_path)


if __name__ == "__main__":
    main()
