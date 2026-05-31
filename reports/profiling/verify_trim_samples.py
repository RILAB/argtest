#!/usr/bin/env python
"""Verify the vectorized single-pass trim_samples against the original
per-interval algorithm on a small ARG with variable sample counts across trees.

The reference implementation below is the ORIGINAL trim_samples code (per
interval: split all edges, drop covered leaf edges, simplify) copied verbatim,
so we compare new-vs-old behavior directly.
"""
import sys
from pathlib import Path

import numpy as np
import msprime
import tskit

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))
from trim_samples import trim_samples_single_pass  # the NEW implementation
from argtest_common import name_to_nodes_map


# ---------------------------------------------------------------------------
# ORIGINAL reference implementation (verbatim from the pre-speedup script)
# ---------------------------------------------------------------------------
def _ref_remove_ancestry(ts, samples, left, right):
    def split_edges_at(tables, position):
        for i, edge in enumerate(tables.edges):
            if edge.left < position < edge.right:
                tables.edges[i] = edge.replace(right=position)
                tables.edges.append(edge.replace(left=position))

    tables = ts.dump_tables()
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
    tables.simplify()
    tables.edges.squash()
    tables.build_index()
    tables.compute_mutation_parents()
    return tables.tree_sequence()


def ref_trim(ts, remove_intervals, suffix_to_strip=""):
    trimmed = ts
    for name, intervals in remove_intervals.items():
        name_to_nodes = name_to_nodes_map(trimmed, suffix_to_strip=suffix_to_strip)
        samples = name_to_nodes.get(name, [])
        if not samples:
            continue
        for left, right in zip(intervals["starts"], intervals["ends"]):
            trimmed = _ref_remove_ancestry(trimmed, samples, left, right)
    return trimmed


# ---------------------------------------------------------------------------
# Build a small ARG with recombination + mutations
# ---------------------------------------------------------------------------
def build_ts():
    ts = msprime.sim_ancestry(
        samples=6,                  # 6 diploids -> 12 sample nodes, individuals ind0..ind5
        ploidy=2,
        sequence_length=100_000,
        recombination_rate=1e-7,
        population_size=10_000,
        random_seed=42,
    )
    ts = msprime.sim_mutations(ts, rate=1e-7, random_seed=7)
    return ts


def make_remove_intervals():
    # Half-open [start, end). Deliberately: overlapping spans for one individual,
    # spans crossing many tree boundaries, multi-individual sharing, and a
    # near-full-genome removal -> variable sample counts across the genome.
    raw = {
        "ind0": [(10_000.0, 30_000.0), (25_000.0, 40_000.5)],   # overlapping -> union
        "ind1": [(0.0, 55_000.3)],                              # big chunk
        "ind2": [(70_000.0, 90_000.0)],                         # shared interval (see ind4)
        "ind4": [(70_000.0, 90_000.0), (5_000.7, 15_000.0)],
        "ind5": [(0.0, 99_999.0)],                              # almost fully isolated
    }
    intervals = {}
    for name, spans in raw.items():
        spans = sorted(spans)
        intervals[name] = {
            "starts": [s for s, _ in spans],
            "ends": [e for _, e in spans],
        }
    return intervals


# ---------------------------------------------------------------------------
# Semantic comparison (robust to node renumbering + isolated samples)
# ---------------------------------------------------------------------------
def pairwise_mrca_times(ts, pos):
    n = ts.num_samples
    samp = ts.samples()
    tree = ts.at(pos)
    out = np.full((n, n), np.inf)
    for a in range(n):
        for b in range(a + 1, n):
            m = tree.mrca(samp[a], samp[b])
            out[a, b] = tree.time(m) if m != tskit.NULL else np.inf
    return out


def isolated_sample_count(ts, pos):
    samp = ts.samples()
    tree = ts.at(pos)
    return int(sum(tree.parent(u) == tskit.NULL and tree.num_children(u) == 0 for u in samp))


def main():
    ts = build_ts()
    print(f"base ts: n_samples={ts.num_samples} trees={ts.num_trees} "
          f"edges={ts.num_edges} sites={ts.num_sites} L={ts.sequence_length:.0f}")

    remove_intervals = make_remove_intervals()

    ts_old = ref_trim(ts, remove_intervals)
    ts_new, summary = trim_samples_single_pass(ts, remove_intervals)
    print(f"trimmed: old trees={ts_old.num_trees} new trees={ts_new.num_trees} | "
          f"old edges={ts_old.num_edges} new edges={ts_new.num_edges}")
    print(f"summary: {summary}")

    ok = True

    def check(cond, msg):
        nonlocal ok
        print(("  PASS " if cond else "  FAIL ") + msg)
        ok = ok and bool(cond)

    # 1. Basic structure
    check(ts_old.sequence_length == ts_new.sequence_length, "sequence_length equal")
    check(ts_old.num_samples == ts_new.num_samples, "num_samples equal")
    check(np.array_equal(ts_old.samples(), ts_new.samples()), "sample node ids equal")

    # 2. Identical tree breakpoints -> identical recombination structure
    bp_old = ts_old.breakpoints(as_array=True)
    bp_new = ts_new.breakpoints(as_array=True)
    check(bp_old.shape == bp_new.shape and np.allclose(bp_old, bp_new),
          f"breakpoints equal (old={len(bp_old)} new={len(bp_new)})")

    # 3. Mutations. Accepted semantics: the single-pass version drops mutations
    #    whose carriers are all masked out (orphans pinned to isolated samples by
    #    the old sequential simplify). So: new sites are a SUBSET of old, genotypes
    #    agree on every shared site, and each old-only site is a genuine orphan
    #    (its derived carriers are all isolated at that position in the old ts).
    pos_old = ts_old.tables.sites.position
    pos_new = ts_new.tables.sites.position
    only_new = np.setdiff1d(pos_new, pos_old)
    only_old = np.setdiff1d(pos_old, pos_new)
    check(only_new.size == 0, f"new sites are a subset of old (new-only={only_new.size})")

    common = np.intersect1d(pos_old, pos_new)
    idx_old = {p: i for i, p in enumerate(pos_old)}
    idx_new = {p: i for i, p in enumerate(pos_new)}
    G_old, G_new = ts_old.genotype_matrix(), ts_new.genotype_matrix()
    geno_mismatch = sum(
        not np.array_equal(G_old[idx_old[p]], G_new[idx_new[p]]) for p in common
    )
    check(geno_mismatch == 0, f"genotypes identical on all {len(common)} shared sites")

    # every old-only site must be a masked-region orphan
    orphan_ok = True
    for p in only_old:
        row = G_old[idx_old[p]]
        tree = ts_old.at(p)
        derived = [s for s in ts_old.samples() if row[list(ts_old.samples()).index(s)] > 0]
        if not all(tree.parent(s) == tskit.NULL for s in derived):
            orphan_ok = False
            break
    check(orphan_ok,
          f"all {only_old.size} old-only sites are masked-region orphans "
          f"(derived allele only on isolated samples)")

    # 4. Genealogy: pairwise MRCA times at the midpoint of every tree
    mids = (bp_new[:-1] + bp_new[1:]) / 2.0
    mismatch = 0
    iso_old_list, iso_new_list = [], []
    for x in mids:
        mo = pairwise_mrca_times(ts_old, x)
        mn = pairwise_mrca_times(ts_new, x)
        finite = np.isfinite(mo) & np.isfinite(mn)
        if not np.array_equal(np.isinf(mo), np.isinf(mn)) or not np.allclose(mo[finite], mn[finite]):
            mismatch += 1
        iso_old_list.append(isolated_sample_count(ts_old, x))
        iso_new_list.append(isolated_sample_count(ts_new, x))
    check(mismatch == 0, f"pairwise MRCA times match at all {len(mids)} tree midpoints")

    # 5. Variable-n is actually exercised, and matches between old and new
    iso_old = np.array(iso_old_list)
    iso_new = np.array(iso_new_list)
    check(np.array_equal(iso_old, iso_new), "isolated-sample counts match per tree")
    check(iso_old.max() > 0 and iso_old.min() != iso_old.max(),
          f"variable n across trees exercised (isolated samples range {iso_old.min()}..{iso_old.max()})")

    # 6. Output validity
    from argtest_common import validate_trimmed_ts
    try:
        validate_trimmed_ts(ts_new)
        check(True, "new ts passes validate_trimmed_ts")
    except Exception as e:
        check(False, f"new ts validate_trimmed_ts: {e}")

    print("\nRESULT:", "ALL CHECKS PASSED" if ok else "FAILURES PRESENT")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
