#!/usr/bin/env python
"""Measure peak RSS of the merge_replicates step on one replicate.

Mirrors merge_treefiles_by_replicate.merge_group: load all chroms for a
replicate, concatenate, then run the metadata-rebuild path (which calls
dump_tables -> another full copy). Reports RSS at each stage plus the process
peak (ru_maxrss). Writes nothing to disk.

Two concatenation strategies can be measured:

    batched      first.concatenate(*rest) in one call, then drop every input
                 reference before dump_tables. This is what v2.0 ships.
    incremental  merged = merged.concatenate(ts) in a loop. This is v1.9.

v2.0 adopted `batched` on the theory that it avoids materialising a
progressively larger intermediate, but that was never measured on real data.
Batching can plausibly be WORSE, since it holds every input tree sequence alive
simultaneously. Run both and keep whichever wins:

    python measure_merge_mem.py --mode incremental --rep 0 --glob '...'
    python measure_merge_mem.py --mode batched     --rep 0 --glob '...'

Each mode must run as its OWN process — ru_maxrss is a high-water mark and never
decreases, so measuring both in one process reports only the larger.

On the cluster, per AGENTS.md, submit rather than running on the head node:
    HPC_MEM=192G ~/.claude/bin/hpc_run 'python measure_merge_mem.py --mode batched'

Legacy positional form (REP, GLOB) is still accepted.
"""
import argparse
import glob
import re
import resource
import sys
import time
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


def natural_key(path: str):
    """Order chromosome paths numerically without assuming a fixed layout depth."""
    parts = re.split(r"(\d+)", str(path))
    return [int(p) if p.isdigit() else p.lower() for p in parts]


def parse_args():
    # Legacy positional form: measure_merge_mem.py [REP] [GLOB]
    argv = sys.argv[1:]
    if argv and not argv[0].startswith("-"):
        rep = argv[0]
        pattern = argv[1] if len(argv) > 1 else None
        return argparse.Namespace(rep=rep, glob=pattern, mode="batched")
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--rep", default="0", help="Replicate id to test (default: 0).")
    p.add_argument("--glob", default=None,
                   help="Glob for one replicate's per-chromosome files.")
    p.add_argument("--mode", choices=["batched", "incremental"], default="batched",
                   help="Concatenation strategy (default: batched, as shipped in v2.0).")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    rep = args.rep
    pattern = args.glob or f"admix/combined.*/trees/combined.*.{rep}.tsz"

    paths = sorted(glob.glob(pattern), key=natural_key)
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
        print(f"{Path(p).parent.name[:13]:<14}{ts.sequence_length/1e6:>8.1f}"
              f"{ts.num_samples:>9}{ts.num_edges:>14,}{ts.nbytes/1e9:>12.2f}")
    print(f"\nwhole-genome table size (sum nbytes): {tot_nbytes/1e9:.2f} GB")
    print(f"RSS after loading all {len(paths)} chroms:      {rss_gb():.2f} GB")

    # Collect the small per-chromosome metadata BEFORE concatenating, so the
    # batched mode can drop every input reference immediately afterwards.
    offsets, cum = [], 0.0
    for ts in tseqs:
        offsets.append(cum)
        cum += float(ts.sequence_length)
    ratemaps = [ratemap_from_metadata(ts.metadata or {}) for ts in tseqs]
    kept = [(ts.metadata or {}).get("kept_intervals") for ts in tseqs]

    t0 = time.perf_counter()
    if args.mode == "batched":
        first, *rest = tseqs
        merged = first.concatenate(*rest) if rest else first
        # Starred unpacking binds independent references to every element, so all
        # four names must go for the inputs to actually be freed.
        del tseqs, rest, first, ts
    else:
        merged = tseqs[0]
        for ts in tseqs[1:]:
            merged = merged.concatenate(ts)
        del tseqs, ts
    concat_s = time.perf_counter() - t0
    print(f"RSS after {args.mode} concatenate ({concat_s:6.1f}s): {rss_gb():.2f} GB"
          f"   (merged nbytes={merged.nbytes/1e9:.2f} GB)")
    # Uses the ratemaps/kept collected before concatenation — the input tree
    # sequences are gone by now, which is the whole point of the batched mode.
    extra: dict = {}
    if all(m is not None for m in ratemaps):
        extra.update(ratemap_to_metadata(merge_ratemaps(ratemaps, offsets)))
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
