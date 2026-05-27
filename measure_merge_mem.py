#!/usr/bin/env python
"""Measure peak RSS of the merge_replicates step on one replicate.

Faithfully mirrors merge_treefiles_by_replicate.merge_group: load all chroms
for a replicate, sequentially concatenate, then run the metadata-rebuild path
(which calls dump_tables -> another full copy). Reports RSS at each stage plus
the process peak (ru_maxrss). Writes nothing to disk.

Usage:
    python measure_merge_mem.py [REP] [GLOB]
        REP   replicate id to test            (default: 0)
        GLOB  glob with a single {rep} field  (default: admix layout below)
"""
import glob
import resource
import sys
from pathlib import Path

sys.path.insert(0, "scripts")
from argtest_common import (  # noqa: E402
    load_ts,
    merge_ratemaps,
    ratemap_from_metadata,
    ratemap_to_metadata,
)
import tskit  # noqa: E402


def rss_gb() -> float:
    with open("/proc/self/status") as f:
        for line in f:
            if line.startswith("VmRSS"):
                return int(line.split()[1]) / 1e6  # kB -> GB
    return float("nan")


def peak_gb() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6  # kB -> GB


def main() -> None:
    rep = sys.argv[1] if len(sys.argv) > 1 else "0"
    pattern = sys.argv[2] if len(sys.argv) > 2 else f"admix/combined.*/trees/combined.*.{rep}.tsz"

    def chrom_num(p: str) -> int:
        return int(Path(p).parts[1].split(".")[1])

    paths = sorted(glob.glob(pattern), key=chrom_num)
    if not paths:
        sys.exit(f"No files matched: {pattern}")
    print(f"{len(paths)} chroms for replicate {rep}\n")

    tseqs = []
    tot_nbytes = 0.0
    print(f"{'chrom':<14}{'Mbp':>8}{'samples':>9}{'edges':>14}{'nbytes(GB)':>12}")
    for p in paths:
        ts = load_ts(Path(p))
        tseqs.append(ts)
        tot_nbytes += ts.nbytes
        print(f"{Path(p).parts[1]:<14}{ts.sequence_length/1e6:>8.1f}"
              f"{ts.num_samples:>9}{ts.num_edges:>14,}{ts.nbytes/1e9:>12.2f}")
    print(f"\nwhole-genome table size (sum nbytes): {tot_nbytes/1e9:.2f} GB")
    print(f"RSS after loading all {len(paths)} chroms:      {rss_gb():.2f} GB")

    merged = tseqs[0]
    for ts in tseqs[1:]:
        merged = merged.concatenate(ts)
    print(f"RSS after sequential concatenate:     {rss_gb():.2f} GB"
          f"   (merged nbytes={merged.nbytes/1e9:.2f} GB)")

    offsets, cum = [], 0.0
    for ts in tseqs:
        offsets.append(cum)
        cum += float(ts.sequence_length)
    extra: dict = {}
    ratemaps = [ratemap_from_metadata(ts.metadata or {}) for ts in tseqs]
    if all(m is not None for m in ratemaps):
        extra.update(ratemap_to_metadata(merge_ratemaps(ratemaps, offsets)))
    kept = [(ts.metadata or {}).get("kept_intervals") for ts in tseqs]
    if all(k is not None for k in kept):
        mk: list = []
        for off, iv in zip(offsets, kept):
            mk.extend([[float(l) + off, float(r) + off] for l, r in iv])
        extra["kept_intervals"] = mk
    print(f"metadata rebuild triggered: {bool(extra)}  (keys={list(extra)})")
    if extra:
        tables = merged.dump_tables()
        existing = merged.metadata if isinstance(merged.metadata, dict) else {}
        tables.metadata_schema = tskit.MetadataSchema({"codec": "json"})
        tables.metadata = {**existing, **extra}
        merged = tables.tree_sequence()
        print(f"RSS after dump_tables rebuild:        {rss_gb():.2f} GB")

    print(f"\n>>> PEAK RSS for the merge: {peak_gb():.2f} GB <<<")
    print(f"merged genome: {merged.sequence_length/1e6:.0f} Mbp, "
          f"{merged.num_samples} samples, {merged.num_edges:,} edges")


if __name__ == "__main__":
    main()
