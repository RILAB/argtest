#!/usr/bin/env python
"""Drop low-sample intervals from a tree sequence (pipeline step between step 5
and merge).

Counts the non-isolated *retained sample nodes* (haploids, not individuals) as
a function of genome position: a sample node is non-isolated where it is the
``child`` of a covering edge. This is computed with a vectorized sweep over
sample-edge endpoints (see ``_sample_edge_coverage``) rather than a per-tree,
per-sample loop; spans with no covering sample edge count 0. Any interval whose
count is below ``--min-samples`` is removed with
``ts.delete_intervals(...)``, which excises those intervals' content but
PRESERVES sequence coordinates (dropped spans become empty gaps; no ``trim()``
/ coordinate compaction). The existing ``kept_intervals`` metadata is then
intersected with the complement of the dropped spans so downstream
accessibility is not overestimated.

Locked design decisions (min_samples_filter_plan.md, 2026-06-26) and the
defaults chosen here for the still-open minor questions:

- Unit is sample NODES (haploids). The ``individuals`` unit mode is deferred.
- Only the ``delete_intervals`` behavior; there is no ``mask`` mode.
- The diagnostic BED of dropped low-sample intervals is emitted ALONGSIDE the
  filtered tree sequence (it is not merged into any existing remove mask).
- Threshold is an absolute ``N`` only (no fraction flag).
- The filter drops at tree-interval granularity; downstream window stats simply
  inherit the resulting coordinate gaps (no special handling).

CLI mirrors trim_samples.py / trim_regions_single.py conventions.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import tskit

from argtest_common import dump_ts, load_ts, merge_intervals, validate_trimmed_ts


def _sample_edge_coverage(ts: tskit.TreeSequence):
    """Coverage step-function of the non-isolated sample count across the genome.

    A sample node is non-isolated at a position iff it is the ``child`` of an
    edge covering that position. Because a child has a unique parent at any
    position, each sample contributes at most one covering edge there, so the
    number of *sample-child* edges covering a position equals the number of
    non-isolated samples. This computes that coverage as an O(E log E) sweep
    over edge endpoints instead of a per-tree, per-sample Python loop.

    Returns ``(seg_left, seg_right, coverage)`` numpy arrays partitioning
    ``[0, sequence_length)`` into maximal constant-coverage segments (spans with
    no covering sample edge have coverage 0).
    """
    seq_len = float(ts.sequence_length)
    edges = ts.tables.edges
    is_sample = np.zeros(ts.num_nodes, dtype=bool)
    is_sample[ts.samples()] = True
    m = is_sample[edges.child]
    left = edges.left[m].astype(float)
    right = edges.right[m].astype(float)

    if left.size == 0:
        return np.array([0.0]), np.array([seq_len]), np.array([0], dtype=int)

    pos = np.concatenate([left, right])
    delta = np.concatenate([np.ones(left.size), -np.ones(right.size)])
    order = np.argsort(pos, kind="mergesort")
    pos = pos[order]
    csum = np.cumsum(delta[order])
    # Collapse ties: coverage on [pos, next_pos) is the cumsum after the LAST
    # delta at that position (a right endpoint and a coinciding left endpoint
    # both apply before the segment to their right).
    is_last = np.empty(pos.size, dtype=bool)
    is_last[-1] = True
    is_last[:-1] = pos[1:] != pos[:-1]
    bp = pos[is_last]
    cov_after = csum[is_last]

    # Partition [0, seq_len) at every breakpoint, plus a leading gap before the
    # first edge and a trailing gap after the last (both coverage 0).
    seg_edges = np.unique(np.concatenate([[0.0, seq_len], bp]))
    seg_edges = seg_edges[(seg_edges >= 0.0) & (seg_edges <= seq_len)]
    seg_left = seg_edges[:-1]
    idx = np.searchsorted(bp, seg_left, side="right") - 1
    coverage = np.where(idx >= 0, cov_after[idx], 0.0).astype(int)

    # Merge adjacent segments carrying the same coverage.
    keep = np.empty(coverage.size, dtype=bool)
    keep[0] = True
    keep[1:] = coverage[1:] != coverage[:-1]
    m_left = seg_left[keep]
    m_cov = coverage[keep]
    m_right = np.append(m_left[1:], seq_len)
    return m_left, m_right, m_cov


def retained_sample_intervals(ts: tskit.TreeSequence):
    """Return ``[(left, right, retained_count), ...]`` over the genome.

    ``retained_count`` is the number of non-isolated sample nodes, given as
    maximal constant-count segments (adjacent trees with the same count are
    merged; spans with no covering sample edge count 0). Computed via the
    vectorized edge-endpoint sweep in ``_sample_edge_coverage``.
    """
    seg_left, seg_right, coverage = _sample_edge_coverage(ts)
    return [
        (float(l), float(r), int(c))
        for l, r, c in zip(seg_left, seg_right, coverage)
    ]


def low_sample_intervals(tree_counts, min_samples):
    """Merge adjacent tree intervals whose retained count is below the threshold.

    Returns a list of merged ``[left, right]`` half-open intervals to drop.
    """
    below = [[l, r] for (l, r, c) in tree_counts if c < min_samples]
    return merge_intervals(below)


def _min_median_max(counts):
    if not counts:
        return (0, 0.0, 0)
    arr = np.asarray(counts, dtype=float)
    return (int(arr.min()), float(np.median(arr)), int(arr.max()))


def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Drop intervals with too few retained sample nodes from a tree "
            "sequence via delete_intervals (coordinates preserved) and update "
            "kept_intervals metadata."
        )
    )
    p.add_argument("--ts", required=True, type=Path, help="Input tree sequence file.")
    p.add_argument(
        "--min-samples",
        required=True,
        type=int,
        help="Minimum non-isolated retained sample nodes per local tree; "
        "intervals below this are dropped.",
    )
    p.add_argument("--out", "--out-ts", dest="out", required=True, type=Path,
                   help="Output (filtered) tree sequence file.")
    p.add_argument(
        "--out-mask",
        dest="out_mask",
        type=Path,
        default=None,
        help="Diagnostic BED of dropped low-sample intervals (columns: chrom "
        "start end retained_samples min_samples). Default: "
        "<out.parent>/<ts_stem>.low_sample.bed.",
    )
    p.add_argument(
        "--chrom",
        default=None,
        help="Chromosome label for the diagnostic BED (default: ts stem).",
    )
    p.add_argument(
        "--log",
        type=Path,
        default=None,
        help="Optional log path (default: <out.parent>/logs/<ts_stem>.filter_min_samples.log).",
    )
    return p.parse_args()


def _intersect_keep_with_drops(keep, dropped, seq_len):
    """Return keep intervals intersected with the complement of dropped spans."""
    keep = merge_intervals([[float(l), float(r)] for l, r in keep])
    dropped = merge_intervals([[float(l), float(r)] for l, r in dropped])
    out = []
    drop_i = 0
    for kl, kr in keep:
        cur = kl
        while drop_i < len(dropped) and dropped[drop_i][1] <= cur:
            drop_i += 1
        j = drop_i
        while j < len(dropped):
            dl, dr = dropped[j]
            if dl >= kr:
                break
            if dl > cur:
                out.append([cur, min(dl, kr)])
            cur = max(cur, dr)
            if cur >= kr:
                break
            j += 1
        if cur < kr:
            out.append([cur, kr])
    return [[l, r] for l, r in out if r > l]


def _min_counts_for_dropped_intervals(dropped, tree_counts):
    """Return the minimum retained count overlapped by each dropped interval."""
    mins = []
    count_i = 0
    for dl, dr in dropped:
        while count_i < len(tree_counts) and tree_counts[count_i][1] <= dl:
            count_i += 1
        j = count_i
        rc = None
        while j < len(tree_counts):
            l, r, c = tree_counts[j]
            if l >= dr:
                break
            if r > dl:
                rc = c if rc is None else min(rc, c)
            j += 1
        mins.append(0 if rc is None else rc)
    return mins


def filter_min_samples(ts: tskit.TreeSequence, min_samples: int):
    """Drop low-sample intervals from ``ts``.

    Returns ``(filtered_ts, dropped_intervals, tree_counts, dropped_bp)`` where
    ``dropped_intervals`` is a list of merged ``[left, right]`` spans and
    ``tree_counts`` is the per-tree ``(left, right, count)`` list.
    """
    tree_counts = retained_sample_intervals(ts)
    dropped = low_sample_intervals(tree_counts, min_samples)
    dropped_bp = float(sum(r - l for l, r in dropped))

    if not dropped:
        return ts, dropped, tree_counts, dropped_bp

    seq_len = float(ts.sequence_length)
    filtered = ts.delete_intervals(dropped, simplify=False)

    # Intersect existing kept_intervals with the complement of the dropped
    # spans so downstream accessibility is not overestimated. If no
    # kept_intervals metadata is present, fall back to the whole sequence.
    tables = filtered.dump_tables()
    existing_md = ts.metadata if ts.metadata is not None else {}
    if isinstance(existing_md, dict) and existing_md.get("kept_intervals"):
        prior_keep = [[float(l), float(r)] for l, r in existing_md["kept_intervals"]]
    else:
        prior_keep = [[0.0, seq_len]]

    new_keep = _intersect_keep_with_drops(prior_keep, dropped, seq_len)

    metadata = dict(existing_md) if isinstance(existing_md, dict) else {}
    metadata["kept_intervals"] = [[float(l), float(r)] for l, r in new_keep]
    tables.metadata_schema = tskit.MetadataSchema({"codec": "json"})
    tables.metadata = metadata
    filtered = tables.tree_sequence()
    return filtered, dropped, tree_counts, dropped_bp


def main():
    args = parse_args()
    if args.min_samples < 0:
        raise SystemExit("ERROR: --min-samples must be >= 0")

    ts = load_ts(args.ts)
    chrom = args.chrom or args.ts.stem

    filtered, dropped, tree_counts, dropped_bp = filter_min_samples(ts, args.min_samples)
    validate_trimmed_ts(filtered)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    dump_ts(filtered, args.out)

    out_mask = args.out_mask or (args.out.parent / f"{args.ts.stem}.low_sample.bed")
    out_mask.parent.mkdir(parents=True, exist_ok=True)
    # For BED annotation, report the smallest retained count spanned by each
    # dropped (merged) interval — the worst offending tree it covers.
    min_counts = _min_counts_for_dropped_intervals(dropped, tree_counts)
    with open(out_mask, "w") as fh:
        for (dl, dr), rc in zip(dropped, min_counts):
            fh.write(f"{chrom}\t{int(dl)}\t{int(dr)}\t{rc}\t{args.min_samples}\n")

    seq_len = float(ts.sequence_length)
    cmin, cmed, cmax = _min_median_max([c for _, _, c in tree_counts])
    frac = dropped_bp / seq_len if seq_len else 0.0
    summary_line = (
        f"filter_min_samples: min_samples={args.min_samples} "
        f"trees={len(tree_counts)} dropped_intervals={len(dropped)} "
        f"dropped_bp={dropped_bp:.0f} dropped_frac={frac:.4f} "
        f"retained_samples(min/median/max)={cmin}/{cmed:.1f}/{cmax} -> out={args.out}"
    )
    print(summary_line)
    print(summary_line, file=sys.stderr)

    log_path = args.log or (args.out.parent / "logs" / f"{args.ts.stem}.filter_min_samples.log")
    try:
        log_path.parent.mkdir(parents=True, exist_ok=True)
        with open(log_path, "w") as fh:
            fh.write("# filter_min_samples summary\n")
            fh.write(f"ts={args.ts}\n")
            fh.write(f"min_samples={args.min_samples}\n")
            fh.write(f"out={args.out}\n")
            fh.write(f"out_mask={out_mask}\n")
            fh.write(f"sequence_length={seq_len}\n")
            fh.write(f"dropped_bp={dropped_bp:.0f}\n")
            fh.write(f"dropped_frac={frac:.6f}\n")
            fh.write(
                f"retained_samples_min={cmin} retained_samples_median={cmed} "
                f"retained_samples_max={cmax}\n"
            )
            fh.write(summary_line + "\n")
    except Exception:
        print(f"WARNING: failed to write log to {log_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
