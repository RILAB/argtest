#!/usr/bin/env python
from __future__ import annotations

import sys
from bisect import bisect_right
from pathlib import Path

import numpy as np
import tskit

try:
    import tszip
except Exception:  # pragma: no cover - optional dependency
    tszip = None


def load_ts(path: Path) -> tskit.TreeSequence:
    # Handle compressed tree sequences if requested.
    if path.suffix == ".tsz":
        if tszip is None:
            raise RuntimeError("tszip is required to load .tsz files")
        return tszip.load(str(path))
    return tskit.load(str(path))


def dump_ts(ts: tskit.TreeSequence, out_path: Path) -> None:
    # Mirror input format when writing if extension is .tsz.
    if out_path.suffix == ".tsz":
        if tszip is None:
            raise RuntimeError("tszip is required to write .tsz files")
        tszip.compress(ts, out_path)
        return
    ts.dump(str(out_path))


def get_individual_name(ind, suffix_to_strip="_anchorwave") -> str:
    # Prefer metadata id when present; otherwise fall back to a stable synthetic name.
    nm = None
    try:
        if isinstance(ind.metadata, dict):
            nm = ind.metadata.get("id")
    except Exception:
        nm = None
    if nm is None:
        nm = f"ind{ind.id}"
    if isinstance(nm, bytes):
        nm = nm.decode()
    return nm.replace(suffix_to_strip, "")


def sample_names(ts: tskit.TreeSequence, suffix_to_strip="_anchorwave"):
    # Map each sample node to its individual name (or node id if missing).
    names = []
    for u in ts.samples():
        ind_id = ts.node(u).individual
        if ind_id != tskit.NULL:
            nm = get_individual_name(ts.individual(ind_id), suffix_to_strip=suffix_to_strip)
        else:
            nm = f"node{u}"
        names.append(nm)
    return names


def aggregate_by_individual(load, names):
    # Collapse per-sample loads into per-individual loads by name.
    unique = []
    idx_map = {}
    for i, nm in enumerate(names):
        if nm not in idx_map:
            idx_map[nm] = len(unique)
            unique.append(nm)

    if load.ndim == 1:
        agg = np.zeros(len(unique), dtype=float)
        for i, nm in enumerate(names):
            agg[idx_map[nm]] += load[i]
        return agg, unique

    agg = np.zeros((load.shape[0], len(unique)), dtype=float)
    for i, nm in enumerate(names):
        agg[:, idx_map[nm]] += load[:, i]
    return agg, unique


def name_to_nodes_map(ts: tskit.TreeSequence, suffix_to_strip="_anchorwave"):
    # Build lookup from individual name to all node ids for that individual.
    mapping = {}
    for ind_id, ind in enumerate(ts.individuals()):
        nm = get_individual_name(ind, suffix_to_strip=suffix_to_strip)
        nodes = list(ts.individual(ind_id).nodes)
        mapping[nm] = nodes
    return mapping


def mutational_load(
    ts: tskit.TreeSequence,
    windows: np.ndarray | None = None,
    remove_intervals=None,
    name_to_nodes=None,
) -> np.ndarray:
    # Compute derived mutation counts per sample, optionally per window and with removals.
    genome_windows = np.array([0, ts.sequence_length]) if windows is None else windows
    assert genome_windows[0] == 0 and genome_windows[-1] == ts.sequence_length
    load = np.zeros((genome_windows.size - 1, ts.num_samples))

    segments = [(0.0, ts.sequence_length, frozenset())]
    if remove_intervals and name_to_nodes:
        # Convert BED intervals into contiguous genome segments with active drop sets.
        segments = build_removal_segments(remove_intervals, ts.sequence_length)

    for left, right, drop_names in segments:
        if right <= left:
            continue
        drop_nodes = set()
        for nm in drop_names:
            drop_nodes.update(name_to_nodes.get(nm, []))

        sub = ts.keep_intervals([(left, right)], simplify=False)
        if drop_nodes:
            # Simplify to retained samples but preserve original ids for load accumulation.
            keep = [u for u in ts.samples() if u not in drop_nodes]
            sub, node_map = sub.simplify(samples=keep, keep_unary=True, map_nodes=True)
            rev_map = np.full(sub.num_nodes, tskit.NULL, dtype=int)
            for u in keep:
                new_id = node_map[u]
                if new_id != tskit.NULL:
                    rev_map[new_id] = u
        else:
            rev_map = None

        site_windows = None
        for tree in sub.trees(sample_lists=True):
            for s in tree.sites():
                if site_windows is None:
                    # Cache site->window indices once per sub-TS.
                    site_windows = np.digitize(sub.sites_position, genome_windows) - 1
                window = int(site_windows[s.id])
                if window < 0 or window >= genome_windows.size - 1:
                    continue
                for m in s.mutations:
                    if m.edge != tskit.NULL:
                        samples = list(tree.samples(m.node))
                        if rev_map is not None:
                            samples = [rev_map[u] for u in samples if rev_map[u] != tskit.NULL]
                        load[window, samples] += 1.0
    return load.squeeze(0) if windows is None else load


def build_removal_segments(remove_intervals, sequence_length):
    # Sweep line to build segments where a set of names is "dropped".
    events = {}
    for name, intervals in remove_intervals.items():
        for start, end in zip(intervals["starts"], intervals["ends"]):
            events.setdefault(start, []).append((name, 1))
            events.setdefault(end, []).append((name, -1))

    positions = sorted(set([0.0, float(sequence_length)] + list(events.keys())))
    active = {}
    segments = []
    for i, pos in enumerate(positions[:-1]):
        for name, delta in events.get(pos, []):
            active[name] = active.get(name, 0) + delta
            if active[name] <= 0:
                active.pop(name, None)
        next_pos = positions[i + 1]
        if next_pos > pos:
            segments.append((pos, next_pos, frozenset(active.keys())))
    return segments


def build_segments_with_drop_nodes(remove_intervals, name_to_nodes, sequence_length):
    # Convert name-based segments into node-id drop sets for trimming.
    segments = build_removal_segments(remove_intervals, sequence_length)
    out = []
    for left, right, drop_names in segments:
        drop_nodes = set()
        for nm in drop_names:
            drop_nodes.update(name_to_nodes.get(nm, []))
        out.append((left, right, drop_nodes))
    return out


def trim_ts_by_intervals(ts: tskit.TreeSequence, remove_intervals, name_to_nodes):
    # Remove edges and mutations that overlap dropped individuals in each segment.
    segments = build_segments_with_drop_nodes(remove_intervals, name_to_nodes, ts.sequence_length)
    seg_lefts = [s[0] for s in segments]
    seg_rights = [s[1] for s in segments]

    tables = ts.dump_tables()
    tables.edges.clear()
    tables.sites.clear()
    tables.mutations.clear()

    # Edges: split by segments, remove if parent/child dropped in that segment
    for edge in ts.edges():
        i = bisect_right(seg_rights, edge.left)
        while i < len(segments) and seg_lefts[i] < edge.right:
            left, right, drop_nodes = segments[i]
            seg_l = max(edge.left, left)
            seg_r = min(edge.right, right)
            if seg_r > seg_l:
                if edge.parent not in drop_nodes and edge.child not in drop_nodes:
                    tables.edges.add_row(
                        left=seg_l,
                        right=seg_r,
                        parent=edge.parent,
                        child=edge.child,
                        metadata=edge.metadata,
                    )
            i += 1

    # Sites and mutations: drop mutations whose node is removed in the segment
    for site in ts.sites():
        i = bisect_right(seg_rights, site.position)
        if i >= len(segments) or seg_lefts[i] > site.position:
            i = bisect_right(seg_lefts, site.position) - 1
        if i < 0 or i >= len(segments):
            continue
        _, _, drop_nodes = segments[i]
        kept = []
        for mut in site.mutations:
            if mut.node in drop_nodes:
                continue
            kept.append(mut)
        if not kept:
            continue
        new_site_id = tables.sites.add_row(
            position=site.position,
            ancestral_state=site.ancestral_state,
            metadata=site.metadata,
        )
        for mut in kept:
            tables.mutations.add_row(
                site=new_site_id,
                node=mut.node,
                derived_state=mut.derived_state,
                parent=tskit.NULL,
                time=mut.time,
                metadata=mut.metadata,
            )

    tables.sort()
    tables.build_index()
    tables.compute_mutation_parents()
    return tables.tree_sequence()


def assert_sample_ids_preserved(orig_ts: tskit.TreeSequence, trimmed_ts: tskit.TreeSequence):
    # Ensure we only removed edges/mutations, not samples or ids.
    orig_samples = orig_ts.samples()
    trimmed_samples = trimmed_ts.samples()
    if len(orig_samples) != len(trimmed_samples):
        msg = "Trimmed tree sequence changed the number of samples"
        print(f"ERROR: {msg}", file=sys.stderr)
        raise RuntimeError(msg)
    if not np.array_equal(orig_samples, trimmed_samples):
        msg = "Trimmed tree sequence changed sample IDs/order"
        print(f"ERROR: {msg}", file=sys.stderr)
        raise RuntimeError(msg)


def validate_trimmed_ts(ts: tskit.TreeSequence):
    # Ensure internal indexes and topology are consistent
    if hasattr(ts, "check_index"):
        try:
            ts.check_index()
        except Exception as e:
            print(f"ERROR: trimmed ts failed check_index: {e}", file=sys.stderr)
            raise
    if hasattr(ts, "validate"):
        try:
            ts.validate()
        except Exception as e:
            print(f"ERROR: trimmed ts failed validate: {e}", file=sys.stderr)
            raise


def load_remove_intervals(paths):
    # Load BED intervals into per-name start/end arrays.
    remove = {}
    for path in paths:
        p = Path(path)
        if not p.exists():
            raise FileNotFoundError(f"Remove BED not found: {p}")
        with open(p, "r") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split()
                if len(parts) < 3:
                    raise ValueError(f"Invalid BED line in {p}: {line}")
                start = float(parts[1])
                end = float(parts[2])
                if end <= start:
                    continue
                if len(parts) >= 4:
                    raw = parts[3]
                    names = [n.strip() for n in raw.split(",") if n.strip()]
                else:
                    names = [p.stem]
                for name in names:
                    remove.setdefault(name, []).append((start, end))

    intervals = {}
    for name, spans in remove.items():
        spans.sort()
        starts = [s for s, _ in spans]
        ends = [e for _, e in spans]
        intervals[name] = {"starts": starts, "ends": ends}
    return intervals
