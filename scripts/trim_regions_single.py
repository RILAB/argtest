#!/usr/bin/env python
from __future__ import annotations

import argparse
from pathlib import Path

import msprime
import tskit

from argtest_common import (
    complement_intervals,
    dump_ts,
    load_mask_intervals,
    load_ts,
    ratemap_to_metadata,
    resolve_mu_rate,
)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Trim one tree sequence using one BED mask and write one output tree sequence."
    )
    parser.add_argument("--ts", required=True, type=Path, help="Input tree sequence file.")
    parser.add_argument("--remove", required=True, type=Path, help="BED file with intervals to remove.")
    parser.add_argument("--out", required=True, type=Path, help="Output tree sequence file.")
    parser.add_argument(
        "--mutation-rate",
        type=float,
        default=None,
        help="Positive scalar fallback when no embedded or exact sibling rate map exists.",
    )
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
    trimmed = ts.keep_intervals(keep, simplify=False)
    tables = trimmed.dump_tables()
    metadata = dict(ts.metadata) if isinstance(ts.metadata, dict) else {}
    metadata["kept_intervals"] = [[float(l), float(r)] for l, r in keep]
    mu = resolve_mu_rate(ts, args.ts, scalar_fallback=args.mutation_rate)
    if isinstance(mu, float):
        if mu <= 0:
            raise ValueError("--mutation-rate must be > 0")
        mu = msprime.RateMap(position=[0.0, float(ts.sequence_length)], rate=[mu])
    metadata.update(ratemap_to_metadata(mu))
    tables.metadata_schema = tskit.MetadataSchema({"codec": "json"})
    tables.metadata = metadata
    trimmed = tables.tree_sequence()

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
