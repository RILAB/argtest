#!/usr/bin/env python
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import tskit

from argtest_common import (
    dump_ts,
    load_remove_intervals,
    load_ts,
    merge_intervals,
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


def merge_full_length(base, extra):
    # Merge per-individual full-length removals into the BED-derived intervals.
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


def _intersect_edges_with_keep(left, right, parent, keep_left, keep_right):
    """Intersect a set of (disjoint) edges of one child with the child's
    ``keep`` intervals, carrying the parent through.

    ``keep_left``/``keep_right`` are the complement of the child's removal
    intervals: sorted, disjoint, half-open ``[keep_left, keep_right)``. Returns
    the surviving edge fragments as ``(left, right, parent)`` arrays. Fully
    vectorized: each edge contributes the (contiguous) run of keep intervals it
    overlaps, clipped to the edge span.
    """
    if keep_left.size == 0 or left.size == 0:
        empty = np.empty(0, dtype=np.float64)
        return empty, empty, np.empty(0, dtype=parent.dtype)

    # Keep interval i overlaps edge j iff keep_right[i] > left[j] and
    # keep_left[i] < right[j]. Because keep intervals are sorted and disjoint,
    # the overlapping i for a given edge form the contiguous range [lo, hi).
    lo = np.searchsorted(keep_right, left, side="right")
    hi = np.searchsorted(keep_left, right, side="left")
    counts = np.maximum(hi - lo, 0)
    total = int(counts.sum())
    if total == 0:
        empty = np.empty(0, dtype=np.float64)
        return empty, empty, np.empty(0, dtype=parent.dtype)

    edge_idx = np.repeat(np.arange(left.size), counts)
    # Offset of each fragment within its edge's run of keep intervals (0,1,...).
    within = np.arange(total) - np.repeat(np.cumsum(counts) - counts, counts)
    keep_idx = lo[edge_idx] + within

    frag_left = np.maximum(left[edge_idx], keep_left[keep_idx])
    frag_right = np.minimum(right[edge_idx], keep_right[keep_idx])
    frag_parent = parent[edge_idx]
    return frag_left, frag_right, frag_parent


def _filter_trimmed_sample_mutations(tables, ts, merged_by_node) -> None:
    """Drop mutations carried by trimmed target sample nodes inside their trim
    intervals.

    Only mutations on *leaf* target nodes are dropped. A target node that acts
    as a parent anywhere in the edge table still has descendant samples attached
    after trimming (only ``child == target`` edges are cut), so a mutation on it
    is inherited by untrimmed samples and must be kept to avoid corrupting their
    genotypes. Sample nodes in this pipeline are leaves, so this guard leaves the
    common case unchanged and only protects the internal-sample-node edge case.

    ``merged_by_node`` maps each target node to its already-merged (sorted,
    disjoint) removal intervals, reused from the edge-trimming pass.
    """
    if not merged_by_node:
        return

    # Restrict dropping to leaf targets: a node that never appears as an edge
    # parent has no descendants, so its mutations are private to it.
    parent_nodes = np.unique(ts.tables.edges.parent)
    target_arr = np.fromiter(merged_by_node.keys(), dtype=np.int64)
    leaf_targets = target_arr[~np.isin(target_arr, parent_nodes)]
    drop_by_node = {
        int(node): merged_by_node[int(node)]
        for node in leaf_targets
        if merged_by_node[int(node)]
    }
    if not drop_by_node:
        return

    mutations = tables.mutations
    num_mutations = mutations.num_rows
    if num_mutations == 0:
        return

    mut_node = mutations.node
    mut_site = mutations.site
    site_position = tables.sites.position
    drop_mutation = np.zeros(num_mutations, dtype=bool)

    leaf_targets = np.fromiter(drop_by_node.keys(), dtype=mut_node.dtype)
    target_mutations = np.flatnonzero(np.isin(mut_node, leaf_targets))
    if target_mutations.size == 0:
        return

    target_order = np.argsort(mut_node[target_mutations], kind="stable")
    target_mutations = target_mutations[target_order]
    target_nodes = mut_node[target_mutations]
    group_starts = np.concatenate(([0], np.flatnonzero(np.diff(target_nodes)) + 1))
    group_ends = np.concatenate((group_starts[1:], [target_mutations.size]))

    for start, end in zip(group_starts, group_ends):
        node = int(target_nodes[start])
        intervals = np.asarray(drop_by_node[node], dtype=np.float64)
        if intervals.size == 0:
            continue
        idx = target_mutations[start:end]
        positions = site_position[mut_site[idx]]
        interval_idx = np.searchsorted(intervals[:, 0], positions, side="right") - 1
        in_bounds = interval_idx >= 0
        in_interval = np.zeros(idx.size, dtype=bool)
        in_interval[in_bounds] = positions[in_bounds] < intervals[interval_idx[in_bounds], 1]
        drop_mutation[idx[in_interval]] = True

    if not np.any(drop_mutation):
        return

    keep_mutation = ~drop_mutation
    keep_site = np.zeros(tables.sites.num_rows, dtype=bool)
    keep_site[mut_site[keep_mutation]] = True

    # Null mutation parents before row deletion. Some kept mutations may have
    # had a dropped mutation parent; compute_mutation_parents below rebuilds the
    # parent column after sorting/simplifying, matching the previous behavior.
    mutations.set_columns(
        site=mutations.site,
        node=mutations.node,
        time=mutations.time,
        derived_state=mutations.derived_state,
        derived_state_offset=mutations.derived_state_offset,
        parent=np.full(num_mutations, tskit.NULL, dtype=mutations.parent.dtype),
        metadata=mutations.metadata,
        metadata_offset=mutations.metadata_offset,
    )
    mutations.keep_rows(keep_mutation)
    site_map = tables.sites.keep_rows(keep_site)
    mutations.set_columns(
        site=site_map[mutations.site],
        node=mutations.node,
        time=mutations.time,
        derived_state=mutations.derived_state,
        derived_state_offset=mutations.derived_state_offset,
        parent=mutations.parent,
        metadata=mutations.metadata,
        metadata_offset=mutations.metadata_offset,
    )


def trim_samples_single_pass(ts, remove_intervals, suffix_to_strip=""):
    """Remove the ancestry of each named individual over its intervals, in a
    single pass.

    Equivalent to dropping, for every sample node, all parent (leaf) edges that
    fall within the union of that individual's removal intervals, then a single
    ``simplify``. Because each per-interval removal is idempotent and
    order-independent, this matches the previous per-interval loop (which
    simplified once per interval) but does O(edges) numpy work plus one
    simplify, instead of one full simplify per interval.

    Returns ``(trimmed_ts, summary)`` where summary has names_removed,
    intervals_applied and sample_nodes_removed counts.
    """
    name_to_nodes = name_to_nodes_map(ts, suffix_to_strip=suffix_to_strip)
    seq_len = float(ts.sequence_length)

    # Map each target node -> its removal intervals (collected across names),
    # mirroring the old code's `samples = name_to_nodes[name]` targeting.
    node_intervals: dict[int, list] = {}
    names_removed = set()
    intervals_applied = 0
    sample_nodes_removed = 0
    for name, spans in remove_intervals.items():
        nodes = name_to_nodes.get(name, [])
        if not nodes:
            continue
        names_removed.add(name)
        sample_nodes_removed += len(nodes)
        intervals_applied += len(spans["starts"])
        pairs = list(zip(spans["starts"], spans["ends"]))
        for node in nodes:
            node_intervals.setdefault(int(node), []).extend(pairs)

    tables = ts.dump_tables()

    if node_intervals:
        # Merge each target's spans once; reused for edge trimming and for
        # dropping the targets' private mutations below.
        merged_by_node = {
            node: merge_intervals(pairs)  # sorted, disjoint [start, end]
            for node, pairs in node_intervals.items()
        }
        edges = tables.edges
        e_left = edges.left
        e_right = edges.right
        e_parent = edges.parent
        e_child = edges.child

        target_children = np.fromiter(node_intervals.keys(), dtype=e_child.dtype)
        is_target = np.isin(e_child, target_children)

        # Non-target edges (child not being trimmed) pass through untouched.
        out_left = [e_left[~is_target]]
        out_right = [e_right[~is_target]]
        out_parent = [e_parent[~is_target]]
        out_child = [e_child[~is_target]]

        for node, merged in merged_by_node.items():
            sel = np.flatnonzero(e_child == node)
            if sel.size == 0:
                continue
            if not merged:
                # No removal spans for this node -> keep all its edges verbatim.
                out_left.append(e_left[sel])
                out_right.append(e_right[sel])
                out_parent.append(e_parent[sel])
                out_child.append(e_child[sel])
                continue
            rem = np.asarray(merged, dtype=np.float64)
            rem_starts = rem[:, 0]
            rem_ends = rem[:, 1]
            # keep = complement of removal within [0, seq_len)
            keep_left = np.concatenate(([0.0], rem_ends))
            keep_right = np.concatenate((rem_starts, [seq_len]))
            nonempty = keep_left < keep_right
            keep_left = keep_left[nonempty]
            keep_right = keep_right[nonempty]

            f_left, f_right, f_parent = _intersect_edges_with_keep(
                e_left[sel], e_right[sel], e_parent[sel], keep_left, keep_right
            )
            out_left.append(f_left)
            out_right.append(f_right)
            out_parent.append(f_parent)
            out_child.append(np.full(f_left.size, node, dtype=e_child.dtype))

        edges.set_columns(
            left=np.concatenate(out_left),
            right=np.concatenate(out_right),
            parent=np.concatenate(out_parent).astype(e_parent.dtype),
            child=np.concatenate(out_child).astype(e_child.dtype),
        )
        _filter_trimmed_sample_mutations(tables, ts, merged_by_node)

    # One canonical pass, mirroring the old per-interval tail.
    tables.sort()
    tables.simplify()
    tables.edges.squash()
    tables.build_index()
    tables.compute_mutation_parents()

    summary = {
        "names_removed": names_removed,
        "intervals_applied": intervals_applied,
        "sample_nodes_removed": sample_nodes_removed,
    }
    return tables.tree_sequence(), summary


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
        remove_intervals = merge_full_length(remove_intervals, full)

    if not args.remove and not args.individuals:
        raise SystemExit("ERROR: provide --individuals and/or --remove")

    trimmed_ts, summary = trim_samples_single_pass(
        ts, remove_intervals, suffix_to_strip=args.suffix_to_strip
    )
    validate_trimmed_ts(trimmed_ts)

    if args.out:
        out_path = Path(args.out)
    else:
        # Default output to a sibling trimmed/ directory with a trimmed suffix.
        default_out_dir.mkdir(parents=True, exist_ok=True)
        out_path = default_out_dir / f"{ts_path.stem}_trimmed.tsz"
    dump_ts(trimmed_ts, out_path)

    # Summary to stdout/stderr and optional log
    summary_line = (
        f"Trimmed: individuals_specified={len(individuals)} "
        f"names_removed={len(summary['names_removed'])} "
        f"intervals_applied={summary['intervals_applied']} "
        f"sample_nodes_removed={summary['sample_nodes_removed']} -> out={out_path}"
    )
    print(summary_line)
    print(summary_line, file=sys.stderr)
    log_path = args.log or (out_path.parent / "logs" / f"{ts_path.stem}_trim_samples.log")
    try:
        log_path.parent.mkdir(parents=True, exist_ok=True)
        with open(log_path, "w") as fh:
            fh.write("# trim_samples summary\n")
            fh.write(summary_line + "\n")
    except Exception:
        print(f"WARNING: failed to write log to {log_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
