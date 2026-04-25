#!/usr/bin/env python
from __future__ import annotations

import argparse
import re
from collections import defaultdict
from pathlib import Path

import tskit

from argtest_common import dump_ts, load_ts, merge_ratemaps, ratemap_from_metadata, ratemap_to_metadata


VALID_SUFFIXES = {".ts", ".trees", ".tsz"}


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
        "--layout",
        choices=["auto", "flat", "nested"],
        default="auto",
        help=(
            "Input layout: 'flat' expects <base>.<chromosome>.<replicate><suffix> files in --ts-dir; "
            "'nested' expects --ts-dir/<chromosome>/<replicate><suffix>; "
            "'auto' detects either layout (default)."
        ),
    )
    p.add_argument(
        "--base-name",
        default=None,
        help=(
            "Base name for nested layout outputs (default: --ts-dir directory name). "
            "Ignored for flat layout."
        ),
    )
    p.add_argument(
        "--out-suffix",
        choices=[".ts", ".trees", ".tsz"],
        default=None,
        help="Optional output suffix. Defaults to the suffix of the first file in each group.",
    )
    p.add_argument(
        "--replicate",
        default=None,
        help=(
            "If set, only write the merged output for this replicate "
            "(default: process all replicates)."
        ),
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
    for path in sorted(ts_dir.rglob(pattern)):
        if path.is_file() and path.suffix in VALID_SUFFIXES:
            files.append(path)
    return files


def parse_nested_input_name(path: Path, ts_dir: Path, base_name: str | None):
    rel = path.relative_to(ts_dir)
    if len(rel.parts) != 2:
        return None
    chrom, filename = rel.parts
    replicate = Path(filename).stem
    if not chrom or not replicate:
        return None
    base = base_name or ts_dir.name
    return base, chrom, replicate


def group_tree_files(paths, ts_dir: Path | None = None, layout: str = "flat", base_name: str | None = None):
    grouped = defaultdict(list)
    skipped = []
    mode_seen = set()
    for path in paths:
        parsed = None
        mode = None
        if ts_dir is None:
            parsed = parse_input_name(path)
            mode = "flat" if parsed is not None else None
        elif layout == "flat":
            rel = path.relative_to(ts_dir)
            if len(rel.parts) == 1:
                parsed = parse_input_name(path)
            mode = "flat" if parsed is not None else None
        elif layout == "nested":
            parsed = parse_nested_input_name(path, ts_dir=ts_dir, base_name=base_name)
            mode = "nested" if parsed is not None else None
        else:
            rel = path.relative_to(ts_dir)
            if len(rel.parts) == 1:
                parsed = parse_input_name(path)
                mode = "flat" if parsed is not None else None
            if parsed is None:
                parsed = parse_nested_input_name(path, ts_dir=ts_dir, base_name=base_name)
                mode = "nested" if parsed is not None else None
        if parsed is None:
            skipped.append(path)
            continue
        mode_seen.add(mode)
        base, chrom, replicate = parsed
        grouped[(base, replicate)].append((chrom, path))
    if layout == "auto" and len(mode_seen) > 1:
        raise ValueError(
            "Detected mixed input layouts (both flat and nested). "
            "Please use a single layout or set --layout explicitly."
        )
    detected = next(iter(mode_seen)) if mode_seen else None
    return grouped, skipped, detected


def merge_group(paths):
    chrom_paths = sorted(paths, key=lambda item: natural_key(item[0]))
    tseqs = [load_ts(p) for _, p in chrom_paths]

    merged = tseqs[0]
    for ts in tseqs[1:]:
        merged = merged.concatenate(ts)

    # tskit's concatenate keeps only the first ts's top-level metadata, so any
    # coordinate-shifted fields (ratemap, kept_intervals) must be re-merged here
    # against the cumulative per-chromosome offsets.
    offsets = []
    cumulative = 0.0
    for ts in tseqs:
        offsets.append(cumulative)
        cumulative += float(ts.sequence_length)

    extra: dict = {}

    ratemaps = [ratemap_from_metadata(ts.metadata or {}) for ts in tseqs]
    if all(mu is not None for mu in ratemaps):
        extra.update(ratemap_to_metadata(merge_ratemaps(ratemaps, offsets)))

    kept_lists = [(ts.metadata or {}).get("kept_intervals") for ts in tseqs]
    if all(k is not None for k in kept_lists):
        merged_kept: list = []
        for off, intervals in zip(offsets, kept_lists):
            merged_kept.extend([[float(l) + off, float(r) + off] for l, r in intervals])
        extra["kept_intervals"] = merged_kept

    if extra:
        tables = merged.dump_tables()
        existing = merged.metadata if isinstance(merged.metadata, dict) else {}
        tables.metadata_schema = tskit.MetadataSchema({"codec": "json"})
        tables.metadata = {**existing, **extra}
        merged = tables.tree_sequence()

    return merged, chrom_paths


def main():
    args = parse_args()
    out_dir = args.out_dir if args.out_dir is not None else args.ts_dir / "combined"
    out_dir.mkdir(parents=True, exist_ok=True)

    files = find_tree_files(args.ts_dir, args.pattern)
    grouped, skipped, detected_layout = group_tree_files(
        files,
        ts_dir=args.ts_dir,
        layout=args.layout,
        base_name=args.base_name,
    )
    if skipped:
        skipped_str = ", ".join(str(path.relative_to(args.ts_dir)) for path in skipped)
        if args.layout == "flat":
            raise ValueError(
                "All input tree files must be named <base>.<chromosome>.<replicate>"
                f"<suffix> directly inside --ts-dir. Skipped: {skipped_str}"
            )
        if args.layout == "nested":
            raise ValueError(
                "All input tree files must be in --ts-dir/<chromosome>/<replicate><suffix>. "
                f"Skipped: {skipped_str}"
            )
        raise ValueError(
            "Could not parse all input files as either flat "
            "(<base>.<chromosome>.<replicate><suffix>) or nested "
            "(<chromosome>/<replicate><suffix>) layout. "
            f"Skipped: {skipped_str}"
        )
    if not grouped:
        raise ValueError(f"No matching tree files found in {args.ts_dir}")
    if args.layout == "auto" and detected_layout is None:
        raise ValueError(f"No parseable tree files found in {args.ts_dir}")

    if args.replicate is not None:
        grouped = {k: v for k, v in grouped.items() if k[1] == args.replicate}
        if not grouped:
            raise ValueError(
                f"No matching tree files found for replicate '{args.replicate}' in {args.ts_dir}"
            )

    for (base, replicate), chrom_paths in sorted(grouped.items()):
        merged, ordered = merge_group(chrom_paths)
        out_suffix = args.out_suffix if args.out_suffix is not None else ordered[0][1].suffix
        out_path = out_dir / f"{base}.combined.{replicate}{out_suffix}"
        dump_ts(merged, out_path)


if __name__ == "__main__":
    main()
