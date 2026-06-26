#!/usr/bin/env python
"""Chunked variant of trim_samples_single_pass.

`simplify()` on a whole 2.7M-tree chromosome is strongly super-linear (~cubic);
small chunks simplify cheaply (5 Mb in ~2s vs ~2h for the whole thing). This
does the same edge surgery as the single-pass version, then simplifies the
genome in chunks with ``filter_nodes=False`` (so node IDs stay GLOBAL across
chunks), and stitches the per-chunk results back together by concatenating the
edge / site / mutation columns and squashing the chunk-boundary edge splits.

Verified equivalent (genealogy + genotypes) to trim_samples_single_pass; see
reports/profiling/verify_trim_chunked.py.

Kept in a SEPARATE module so editing it never disturbs the production
scripts/trim_samples.py that running step5 jobs import.
"""
from __future__ import annotations

import numpy as np

from argtest_common import merge_intervals, name_to_nodes_map
from trim_samples import _intersect_edges_with_keep


def _surgery_tables(ts, remove_intervals, suffix_to_strip=""):
    """Edge surgery only (no simplify) -- identical to the body of
    trim_samples_single_pass up to (but not including) tables.sort()."""
    name_to_nodes = name_to_nodes_map(ts, suffix_to_strip=suffix_to_strip)
    seq_len = float(ts.sequence_length)

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
        edges = tables.edges
        e_left, e_right, e_parent, e_child = edges.left, edges.right, edges.parent, edges.child
        target = np.fromiter(node_intervals.keys(), dtype=e_child.dtype)
        is_target = np.isin(e_child, target)
        out_left = [e_left[~is_target]]
        out_right = [e_right[~is_target]]
        out_parent = [e_parent[~is_target]]
        out_child = [e_child[~is_target]]
        for node, pairs in node_intervals.items():
            sel = np.flatnonzero(e_child == node)
            if sel.size == 0:
                continue
            merged = merge_intervals(pairs)
            if not merged:
                out_left.append(e_left[sel]); out_right.append(e_right[sel])
                out_parent.append(e_parent[sel]); out_child.append(e_child[sel])
                continue
            rem = np.asarray(merged, dtype=np.float64)
            keep_left = np.concatenate(([0.0], rem[:, 1]))
            keep_right = np.concatenate((rem[:, 0], [seq_len]))
            ne = keep_left < keep_right
            keep_left, keep_right = keep_left[ne], keep_right[ne]
            fL, fR, fP = _intersect_edges_with_keep(e_left[sel], e_right[sel], e_parent[sel], keep_left, keep_right)
            out_left.append(fL); out_right.append(fR); out_parent.append(fP)
            out_child.append(np.full(fL.size, node, dtype=e_child.dtype))
        edges.set_columns(
            left=np.concatenate(out_left), right=np.concatenate(out_right),
            parent=np.concatenate(out_parent).astype(e_parent.dtype),
            child=np.concatenate(out_child).astype(e_child.dtype),
        )
    summary = {"names_removed": names_removed, "intervals_applied": intervals_applied,
               "sample_nodes_removed": sample_nodes_removed}
    return tables, summary


def _concat_ragged(data_list, off_list):
    """Concatenate ragged (data, offset) column pairs into one (data, offset)."""
    data = np.concatenate(data_list) if data_list else np.array([], dtype=np.int8)
    offs = [np.array([0], dtype=np.uint64)]
    total = np.uint64(0)
    for off in off_list:
        offs.append(off[1:].astype(np.uint64) + total)
        total += off[-1]
    return data, np.concatenate(offs)


def trim_samples_chunked(ts, remove_intervals, suffix_to_strip="", chunk_size=1.0e7):
    tables, summary = _surgery_tables(ts, remove_intervals, suffix_to_strip)
    tables.sort()
    tables.build_index()
    tables.compute_mutation_parents()
    tts = tables.tree_sequence()
    seq_len = float(tts.sequence_length)

    eL, eR, eP, eC = [], [], [], []
    spos = []
    sanc_data, sanc_off = [], []
    msite, mnode, mtime = [], [], []
    mder_data, mder_off = [], []
    site_offset = 0

    a = 0.0
    while a < seq_len:
        b = min(a + chunk_size, seq_len)
        ct = tts.keep_intervals([[a, b]], simplify=False).dump_tables()
        ct.simplify(filter_nodes=False)  # cheap on a small region; keeps global node IDs
        e = ct.edges
        eL.append(e.left); eR.append(e.right); eP.append(e.parent); eC.append(e.child)
        s = ct.sites
        spos.append(s.position)
        sanc_data.append(s.ancestral_state); sanc_off.append(s.ancestral_state_offset)
        m = ct.mutations
        msite.append(m.site.astype(np.int64) + site_offset)
        mnode.append(m.node); mtime.append(m.time)
        mder_data.append(m.derived_state); mder_off.append(m.derived_state_offset)
        site_offset += s.num_rows
        a = b

    final = tts.dump_tables()  # full node/individual/population tables (filter_nodes=False kept all)
    final.edges.set_columns(
        left=np.concatenate(eL), right=np.concatenate(eR),
        parent=np.concatenate(eP), child=np.concatenate(eC),
    )
    anc_data, anc_off = _concat_ragged(sanc_data, sanc_off)
    final.sites.set_columns(position=np.concatenate(spos),
                            ancestral_state=anc_data, ancestral_state_offset=anc_off)
    der_data, der_off = _concat_ragged(mder_data, mder_off)
    msite_all = np.concatenate(msite).astype(np.int32)
    final.mutations.set_columns(
        site=msite_all, node=np.concatenate(mnode).astype(np.int32),
        time=np.concatenate(mtime),
        derived_state=der_data, derived_state_offset=der_off,
        parent=np.full(msite_all.size, -1, dtype=np.int32),
    )
    final.sort()
    final.edges.squash()
    final.build_index()
    final.compute_mutation_parents()
    return final.tree_sequence(), summary
