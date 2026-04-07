#!/usr/bin/env python
from __future__ import annotations

import argparse
import re
from collections import defaultdict
from pathlib import Path

from argtest_common import dump_ts, load_ts


VALID_SUFFIXES = {".tree", ".trees", ".tsz"}


def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Merge chromosome-specific tree sequence files by replicate. "
            "Inputs must be named <base>.<chromosome>.<replicate><suffix>."
        )
    )
    p.add_argument(
        "--ts-dir",
        required=True,
        type=Path,
        help="Directory containing chromosome-specific tree sequence files.",
    )
    p.add_argument(
        "--out-dir",
        type=Path,
        default=None,
        help="Output directory for merged tree sequences (default: <ts-dir>/combined).",
    )
    p.add_argument(
        "--pattern",
        default="*",
        help="Optional glob pattern to filter input filenames (default: '*').",
    )
    p.add_argument(
        "--out-suffix",
        choices=[".tree", ".trees", ".tsz"],
        default=None,
        help="Optional output suffix. Defaults to the suffix of the first file in each group.",
    )
    return p.parse_args()


def natural_key(text: str):
    return [int(part) if part.isdigit() else part.lower() for part in re.split(r"(\d+)", text)]


def parse_input_name(path: Path):
    if path.suffix not in VALID_SUFFIXES:
        return None
    parts = path.stem.rsplit(".", 2)
    if len(parts) != 3:
        return None
    base, chrom, replicate = parts
    if not base or not chrom or not replicate:
        return None
    return base, chrom, replicate


def find_tree_files(ts_dir: Path, pattern: str):
    if not ts_dir.exists():
        raise FileNotFoundError(f"Tree directory does not exist: {ts_dir}")
    if not ts_dir.is_dir():
        raise NotADirectoryError(f"Expected a directory for --ts-dir: {ts_dir}")
    files = []
    for path in sorted(ts_dir.glob(pattern)):
        if path.is_file() and path.suffix in VALID_SUFFIXES:
            files.append(path)
    return files


def group_tree_files(paths):
    grouped = defaultdict(list)
    skipped = []
    for path in paths:
        parsed = parse_input_name(path)
        if parsed is None:
            skipped.append(path)
            continue
        base, chrom, replicate = parsed
        grouped[(base, replicate)].append((chrom, path))
    return grouped, skipped


def merge_group(paths):
    chrom_paths = sorted(paths, key=lambda item: natural_key(item[0]))
    merged = load_ts(chrom_paths[0][1])
    for _, path in chrom_paths[1:]:
        merged = merged.concatenate(load_ts(path))
    return merged, chrom_paths


def main():
    args = parse_args()
    out_dir = args.out_dir if args.out_dir is not None else args.ts_dir / "combined"
    out_dir.mkdir(parents=True, exist_ok=True)

    files = find_tree_files(args.ts_dir, args.pattern)
    grouped, skipped = group_tree_files(files)
    if skipped:
        skipped_str = ", ".join(path.name for path in skipped)
        raise ValueError(
            "All input tree files must be named <base>.<chromosome>.<replicate>"
            f"<suffix>. Skipped: {skipped_str}"
        )
    if not grouped:
        raise ValueError(f"No matching tree files found in {args.ts_dir}")

    for (base, replicate), chrom_paths in sorted(grouped.items()):
        merged, ordered = merge_group(chrom_paths)
        out_suffix = args.out_suffix if args.out_suffix is not None else ordered[0][1].suffix
        out_path = out_dir / f"{base}.combined.{replicate}{out_suffix}"
        dump_ts(merged, out_path)


if __name__ == "__main__":
    main()
