#!/usr/bin/env python
from __future__ import annotations

import math
import pickle
from collections import Counter
from pathlib import Path

import numpy as np
import tskit

try:
    import tszip
except Exception:  # pragma: no cover - optional dependency
    tszip = None

try:
    import msprime
except Exception:  # pragma: no cover - optional dependency
    msprime = None


def _require_msprime():
    if msprime is None:
        raise RuntimeError("msprime is required for accessibility mask and collapse helpers")


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


def get_individual_name(ind, name_substring_to_remove="") -> str:
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
    return nm.replace(name_substring_to_remove, "")


def sample_names(ts: tskit.TreeSequence, name_substring_to_remove=""):
    # Map each sample node to its individual name (or node id if missing).
    names = []
    for u in ts.samples():
        ind_id = ts.node(u).individual
        if ind_id != tskit.NULL:
            nm = get_individual_name(
                ts.individual(ind_id),
                name_substring_to_remove=name_substring_to_remove,
            )
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


def name_to_nodes_map(ts: tskit.TreeSequence, name_substring_to_remove=""):
    # Build lookup from individual name to all node ids for that individual.
    mapping = {}
    for ind_id, ind in enumerate(ts.individuals()):
        nm = get_individual_name(
            ind, name_substring_to_remove=name_substring_to_remove
        )
        nodes = list(ts.individual(ind_id).nodes)
        mapping[nm] = nodes
    return mapping


def audit_individual_contract(ts, name_substring_to_remove=""):
    """Return warnings about the proposed individual/ploidy input contract.

    This is deliberately diagnostic: callers may report the returned strings,
    but pipeline behavior must not depend on the result.
    """
    samples = {int(u) for u in ts.samples()}
    unassigned = sorted(u for u in samples if ts.node(u).individual == tskit.NULL)

    represented_names = []
    represented_ploidies = []
    for ind in ts.individuals():
        sample_nodes = [int(u) for u in ind.nodes if int(u) in samples]
        if not sample_nodes:
            continue
        represented_names.append(
            get_individual_name(
                ind, name_substring_to_remove=name_substring_to_remove
            )
        )
        represented_ploidies.append(len(sample_nodes))

    warnings = []
    if unassigned:
        warnings.append(f"{len(unassigned)} sample node(s) have no individual")

    duplicate_names = sorted(
        name for name, count in Counter(represented_names).items() if count > 1
    )
    if duplicate_names:
        warnings.append(
            "duplicate normalized individual name(s): " + ", ".join(duplicate_names)
        )

    ploidies = sorted(set(represented_ploidies))
    if len(ploidies) > 1:
        warnings.append(
            "mixed ploidy among represented individuals: "
            + ", ".join(map(str, ploidies))
        )

    parent_samples = sorted(samples.intersection(map(int, ts.edges_parent)))
    if parent_samples:
        warnings.append(f"{len(parent_samples)} sample node(s) are edge parents")
    return warnings


def mutational_load(
    ts: tskit.TreeSequence,
    windows: np.ndarray | None = None,
) -> np.ndarray:
    # Compute derived mutation counts per sample, optionally per window.
    genome_windows = np.array([0, ts.sequence_length]) if windows is None else windows
    if genome_windows[0] != 0 or genome_windows[-1] != ts.sequence_length:
        raise ValueError(
            "mutation-load windows must start at 0 and end at the tree-sequence length"
        )
    load = np.zeros((genome_windows.size - 1, ts.num_samples))
    site_windows = None
    for tree in ts.trees(sample_lists=True):
        for s in tree.sites():
            if site_windows is None:
                site_windows = np.digitize(ts.sites_position, genome_windows) - 1
            window = int(site_windows[s.id])
            if window < 0 or window >= genome_windows.size - 1:
                continue
            for m in s.mutations:
                if m.edge != tskit.NULL:
                    load[window, list(tree.samples(m.node))] += 1.0
    return load.squeeze(0) if windows is None else load


def build_bp_windows(ts_or_length, window_size: float) -> np.ndarray:
    """Build a full-sequence grid of fixed-width base-pair windows."""
    if window_size <= 0:
        raise ValueError("window size must be > 0")
    sequence_length = float(
        getattr(ts_or_length, "sequence_length", ts_or_length)
    )
    windows = np.arange(0, sequence_length + window_size, window_size, dtype=float)
    windows[-1] = sequence_length
    return windows


def build_snp_windows(ts, snp_window: int) -> np.ndarray:
    """Build a full-sequence grid with an edge after every Nth site."""
    if snp_window <= 0:
        raise ValueError("SNP window size must be > 0")
    positions = np.asarray(ts.sites_position, dtype=float)
    sequence_length = float(ts.sequence_length)
    edges = positions[snp_window::snp_window]
    return np.concatenate(([0.0], edges, [sequence_length]))


def simulate_expected_load(ts, windows, sample_names, mutation_rate, seed):
    """Return one seeded, simulation-based expected load per individual."""
    _require_msprime()
    sim_ts = msprime.sim_mutations(
        ts, rate=mutation_rate, keep=False, random_seed=seed
    )
    expected = mutational_load(sim_ts, windows=windows)
    expected, _ = aggregate_by_individual(expected, sample_names)
    return expected


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


def merge_intervals(intervals):
    # Merge overlapping or adjacent half-open intervals [left, right).
    if len(intervals) == 0:
        return []
    intervals = np.asarray(intervals, dtype=float)
    intervals = intervals[np.argsort(intervals[:, 0])]
    merged = [intervals[0].tolist()]
    for left, right in intervals[1:]:
        if left <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], right)
        else:
            merged.append([left, right])
    return merged


def load_mask_intervals(path: Path, sequence_length: float):
    intervals = []
    with open(path, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 3:
                raise ValueError(f"Invalid BED line in {path}: {line}")
            start = max(0.0, float(parts[1]))
            end = min(float(sequence_length), float(parts[2]))
            if end > start:
                intervals.append([start, end])
    return merge_intervals(intervals)


def complement_intervals(intervals, sequence_length: float):
    keep = []
    last = 0.0
    for left, right in intervals:
        if left > last:
            keep.append([last, left])
        last = max(last, right)
    if last < float(sequence_length):
        keep.append([last, float(sequence_length)])
    return keep


def chrom_offsets_from_metadata(ts: tskit.TreeSequence):
    """Return the recorded [{chrom, offset, length}, ...] table for a merged TS, or None.

    Populated by merge_treefiles_by_replicate.py. Older merged files predating that
    field return None.
    """
    md = ts.metadata if isinstance(ts.metadata, dict) else {}
    return md.get("chrom_offsets")


def genome_position(ts: tskit.TreeSequence, chrom, position: float) -> float:
    """Map (chromosome, within-chromosome position) to a coordinate in a merged TS.

    Uses the ``chrom_offsets`` metadata written at merge time, so no reference-genome
    lengths or manual offset arithmetic are needed. Raises KeyError if the metadata is
    absent (re-run the merge to record it) or the chromosome is unknown, and ValueError
    if the position falls outside that chromosome.
    """
    entries = chrom_offsets_from_metadata(ts)
    if not entries:
        raise KeyError(
            "Tree sequence has no 'chrom_offsets' metadata; re-run "
            "merge_treefiles_by_replicate.py on the inputs to record it."
        )
    chrom = str(chrom)
    for e in entries:
        if str(e["chrom"]) == chrom:
            length = float(e["length"])
            pos = float(position)
            if not math.isfinite(pos) or pos < 0 or pos >= length:
                raise ValueError(
                    f"position {position} out of range [0, {length}) for chromosome {chrom}"
                )
            return float(e["offset"]) + pos
    known = ", ".join(str(e["chrom"]) for e in entries)
    raise KeyError(f"chromosome {chrom!r} not found in merged tree sequence; known: {known}")


def chrom_position_from_genome(ts: tskit.TreeSequence, genome_pos: float):
    """Inverse of genome_position: map a merged coordinate to (chrom, within-chrom pos)."""
    entries = chrom_offsets_from_metadata(ts)
    if not entries:
        raise KeyError(
            "Tree sequence has no 'chrom_offsets' metadata; re-run "
            "merge_treefiles_by_replicate.py on the inputs to record it."
        )
    g = float(genome_pos)
    for e in entries:
        off = float(e["offset"])
        length = float(e["length"])
        if off <= g < off + length:
            return str(e["chrom"]), g - off
    raise ValueError(f"genome position {genome_pos} is outside the merged sequence")


def tree_at_chrom_position(ts: tskit.TreeSequence, chrom, position: float):
    """Return the tree covering (chromosome, within-chromosome position) in a merged TS.

    If the position lies in a masked/trimmed region the returned tree has no topology
    (``tree.num_edges == 0``); check that if a real local tree is required.
    """
    return ts.at(genome_position(ts, chrom, position))


def ratemap_to_metadata(mu) -> dict:
    """Serialize a RateMap to a JSON-safe dict for storage in tree-sequence metadata."""
    return {
        "mu_position": list(float(x) for x in mu.position),
        "mu_rate": list(float(x) for x in mu.rate),
    }


def ratemap_from_metadata(md: dict):
    """Reconstruct a RateMap from metadata produced by ratemap_to_metadata, or None."""
    if not md or "mu_position" not in md or "mu_rate" not in md:
        return None
    _require_msprime()
    return msprime.RateMap(position=md["mu_position"], rate=md["mu_rate"])


def resolve_mu_rate(ts: tskit.TreeSequence, ts_path: Path, scalar_fallback: float | None = None):
    """Resolve a mutation rate for `msprime.sim_mutations` from (in order):
    ts metadata ratemap, *.mut_rate.p sibling file, or a scalar fallback.
    Returns an msprime.RateMap or float; raises if nothing is available.
    """
    _require_msprime()
    md_rate = ratemap_from_metadata((ts.metadata or {}) if ts.metadata is not None else {})
    if md_rate is not None:
        return md_rate
    try:
        mu_path = infer_mu_path(ts_path)
    except FileNotFoundError:
        mu_path = None
    if mu_path is not None:
        with open(mu_path, "rb") as fh:
            obj = pickle.load(fh)
        if isinstance(obj, msprime.RateMap):
            return obj
        raise RuntimeError(f"Unrecognized mutation-rate object in {mu_path}: {type(obj)}")
    if scalar_fallback is not None:
        return float(scalar_fallback)
    raise FileNotFoundError(
        f"No mutation rate available for {ts_path}: ts.metadata has no ratemap, "
        f"no *.mut_rate.p file found, and no --mutation-rate fallback was given."
    )


def merge_ratemaps(ratemaps: list, offsets: list):
    """Concatenate per-chromosome RateMaps into one by shifting each by its offset."""
    _require_msprime()
    all_pos: list[float] = []
    all_rate: list[float] = []
    for mu, offset in zip(ratemaps, offsets):
        pos = [float(x) + offset for x in mu.position]
        rate = [float(x) for x in mu.rate]
        if all_pos:
            all_pos.extend(pos[1:])  # drop duplicate junction point
        else:
            all_pos.extend(pos)
        all_rate.extend(rate)
    return msprime.RateMap(position=all_pos, rate=all_rate)


def accessible_intervals_from_mu(mu):
    # Convert mutation-rate map into explicit accessible intervals.
    pos = np.asarray(mu.position, dtype=float)
    rate = np.asarray(mu.rate)
    keep = rate > 0
    lefts = pos[:-1][keep]
    rights = pos[1:][keep]
    return np.column_stack([lefts, rights])


def tree_covered_accessible_bp(ts: tskit.TreeSequence, acc_intervals=None) -> float:
    # Total accessible bp that is also covered by a non-empty tree (num_edges > 0).
    # Matches singer-snakemake's extract_accessible_ratemap approach.
    # acc_intervals: array-like of [left, right] pairs, or None to use full sequence.
    if acc_intervals is None:
        return float(sum(
            t.interval.right - t.interval.left
            for t in ts.trees() if t.num_edges > 0
        ))
    acc = np.asarray(acc_intervals, dtype=float)
    total = 0.0
    ai = 0
    n = len(acc)
    for tree in ts.trees():
        if tree.num_edges == 0:
            continue
        tl, tr = float(tree.interval.left), float(tree.interval.right)
        while ai < n and acc[ai, 1] <= tl:
            ai += 1
        j = ai
        while j < n and acc[j, 0] < tr:
            total += min(acc[j, 1], tr) - max(acc[j, 0], tl)
            j += 1
    return total


def infer_mu_path(ts_path: Path) -> Path:
    """Find a mutation map using a small, deterministic set of exact paths."""
    parent = ts_path.parent
    names = [f"{ts_path.stem}.mut_rate.p", f"{parent.name}.mut_rate.p"]
    tried = []
    seen = set()
    for directory in (parent, parent.parent):
        for name in names:
            candidate = directory / name
            if candidate in seen:
                continue
            seen.add(candidate)
            tried.append(candidate)
            if candidate.is_file():
                return candidate
    raise FileNotFoundError(
        f"Could not infer mutation map for {ts_path}. Tried exact paths: "
        + ", ".join(str(path) for path in tried)
    )


def overlap_lengths(intervals, windows):
    # Return total overlap length per window for sorted half-open intervals.
    intervals = np.asarray(intervals, dtype=float)
    windows = np.asarray(windows, dtype=float)
    k = len(windows) - 1
    out = np.zeros(k, dtype=float)
    i = 0
    m = len(intervals)
    for w in range(k):
        wl, wr = windows[w], windows[w + 1]
        while i < m and intervals[i, 1] <= wl:
            i += 1
        j = i
        while j < m and intervals[j, 0] < wr:
            ol = max(wl, intervals[j, 0])
            or_ = min(wr, intervals[j, 1])
            if or_ > ol:
                out[w] += or_ - ol
            if intervals[j, 1] >= wr:
                break
            j += 1
    return out
