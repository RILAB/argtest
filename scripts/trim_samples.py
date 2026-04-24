#!/usr/bin/env python
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

from argtest_common import (
    dump_ts,
    load_remove_intervals,
    load_ts,
    name_to_nodes_map,
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


def parse_individuals(values):
    # Parse comma-separated IDs into a list.
    if not values:
        return []
    return [v.strip() for v in values.split(",") if v.strip()]


def merge_intervals(base, extra):
    # Merge per-individual intervals, keeping them sorted.
    merged = {k: {"starts": list(v["starts"]), "ends": list(v["ends"])} for k, v in base.items()}
    for name, spans in extra.items():
        entry = merged.setdefault(name, {"starts": [], "ends": []})
        entry["starts"].extend(spans["starts"])
        entry["ends"].extend(spans["ends"])
    for name, spans in merged.items():
        paired = sorted(zip(spans["starts"], spans["ends"]))
        spans["starts"] = [s for s, _ in paired]
        spans["ends"] = [e for _, e in paired]
    return merged


def parse_args():
    p = argparse.ArgumentParser(
        description="Remove individuals over BED intervals and write a trimmed tree sequence",
    )
    p.add_argument("ts", help="Tree sequence file (.ts, .trees, or .tsz)")
    p.add_argument(
        "--individuals",
        help="Comma-separated individual IDs to remove across the entire sequence",
    )
    p.add_argument(
        "--remove",
        action="append",
        help=(
            "BED file of per-individual intervals to remove. Column 4 should "
            "contain the sample ID (comma-separated for multiple IDs sharing "
            "the same interval); if column 4 is absent the BED filename stem "
            "is used. Can be repeated or given as a comma-separated list to "
            "supply multiple BED files."
        ),
    )
    p.add_argument(
        "--out",
        help="Output tree sequence path (default: results/<ts_stem>_trimmed.tsz).",
    )
    p.add_argument(
        "--suffix-to-strip",
        default="",
        help='Suffix removed from sample names before matching (default: "").',
    )
    p.add_argument(
        "--log",
        type=Path,
        default=None,
        help="Optional log file path (default: <out.parent>/logs/<ts_stem>_trim_samples.log).",
    )
    return p.parse_args()


def remove_ancestry(ts, samples, left, right):
    """
    Removes the ancestry for `samples` over `[left, right)`, by:
    1. Split all edges intersecting "left" and "right"
    2. Remove singleton edges above nodes for which we don't want ancestry over [left, right]
    3. Throw into simplify, which will remove this ancestry
    4. Squash edges to "join" previously split edges
    """

    def split_edges_at(tables, position):
        for i, edge in enumerate(tables.edges):
            if edge.left < position < edge.right:
                tables.edges[i] = edge.replace(right=position)
                tables.edges.append(edge.replace(left=position))

    tables = ts.dump_tables()
    # Split edges so we can drop the target interval exactly.
    split_edges_at(tables, left)
    split_edges_at(tables, right)
    drop_edges = np.logical_and.reduce(
        [
            np.isin(tables.edges.child, samples),
            tables.edges.left >= left,
            tables.edges.right <= right,
        ]
    )
    tables.edges.keep_rows(~drop_edges)
    tables.sort()
    tables.edges.drop_metadata()
    # Simplify drops the removed ancestry and may renumber nodes.
    tables.simplify()
    tables.edges.squash()
    tables.build_index()
    tables.compute_mutation_parents()
    return tables.tree_sequence()


def main():
    args = parse_args()
    ts_path = Path(args.ts)
    ts = load_ts(ts_path)
    default_out_dir = ts_path.parent / "trimmed"

    remove_intervals = {}
    if args.remove:
        remove_paths = parse_remove_list(args.remove)
        remove_intervals = load_remove_intervals(remove_paths)

    individuals = parse_individuals(args.individuals)
    if individuals:
        # Expand full-length removals to [0, sequence_length).
        full = {
            name: {"starts": [0.0], "ends": [float(ts.sequence_length)]}
            for name in individuals
        }
        remove_intervals = merge_intervals(remove_intervals, full)

    if not args.remove and not args.individuals:
        raise SystemExit("ERROR: provide --individuals and/or --remove")


    # Track what we removed for a brief summary.
    trimmed_ts = ts
    names_removed = set()
    intervals_applied = 0
    sample_nodes_removed = 0
    for name, intervals in remove_intervals.items():
        # Rebuild the name->nodes map after each simplify.
        name_to_nodes = name_to_nodes_map(trimmed_ts, suffix_to_strip=args.suffix_to_strip)
        samples = name_to_nodes.get(name, [])
        if not samples:
            continue
        names_removed.add(name)
        sample_nodes_removed += len(samples)
        intervals_applied += len(intervals["starts"])
        for left, right in zip(intervals["starts"], intervals["ends"]):
            trimmed_ts = remove_ancestry(trimmed_ts, samples, left, right)
    validate_trimmed_ts(trimmed_ts)

    if args.out:
        out_path = Path(args.out)
    else:
        # Default output to a sibling trimmed/ directory with a trimmed suffix.
        default_out_dir.mkdir(parents=True, exist_ok=True)
        out_path = default_out_dir / f"{ts_path.stem}_trimmed.tsz"
    dump_ts(trimmed_ts, out_path)

    # Summary to stdout/stderr and optional log
    summary = (
        f"Trimmed: individuals_specified={len(parse_individuals(args.individuals) if args.individuals else [])} "
        f"names_removed={len(names_removed)} intervals_applied={intervals_applied} sample_nodes_removed={sample_nodes_removed} -> out={out_path}"
    )
    print(summary)
    print(summary, file=sys.stderr)
    log_path = args.log or (out_path.parent / "logs" / f"{ts_path.stem}_trim_samples.log")
    try:
        log_path.parent.mkdir(parents=True, exist_ok=True)
        with open(log_path, "w") as fh:
            fh.write("# trim_samples summary\n")
            fh.write(summary + "\n")
    except Exception:
        print(f"WARNING: failed to write log to {log_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
