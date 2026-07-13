#!/usr/bin/env python
"""Locate the local tree at a (chromosome, position) in a merged tree sequence.

The merged tree sequence lays chromosomes end-to-end; this maps a within-chromosome
coordinate to the concatenated coordinate using the ``chrom_offsets`` metadata written
by merge_treefiles_by_replicate.py, then reports the tree covering that site.

Example:
    python scripts/locate_tree.py --ts combined/amaranth.combined.1.trees --chrom 8 --position 1234
"""
from __future__ import annotations

import argparse
from pathlib import Path

from argtest_common import genome_position, load_ts


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--ts", required=True, type=Path, help="Merged tree sequence file.")
    p.add_argument("--chrom", required=True, help="Chromosome name as used at merge time (e.g. 8).")
    p.add_argument("--position", required=True, type=float, help="Within-chromosome position.")
    return p.parse_args()


def main():
    args = parse_args()
    ts = load_ts(args.ts)
    gpos = genome_position(ts, args.chrom, args.position)
    tree = ts.at(gpos)
    masked = tree.num_edges == 0

    print(f"chromosome:        {args.chrom}")
    print(f"chrom position:    {args.position}")
    print(f"genome position:   {gpos}")
    print(f"tree index:        {tree.index}")
    print(f"tree interval:     [{tree.interval.left}, {tree.interval.right})")
    print(f"tree num_edges:    {tree.num_edges}")
    print(f"tree num_roots:    {tree.num_roots}")
    if masked:
        print("NOTE: this position is in a masked/trimmed region (no local tree topology).")


if __name__ == "__main__":
    main()
