#!/usr/bin/env python
"""Dry-run the v1.10 mutation-map discovery over a dataset. Loads no ARG data.

v1.10 narrowed ``infer_mu_path`` to exact candidates and made step 4 fail when no
rate resolves, so a layout that worked under v1.9's fuzzy matching can now break.
This reports, per (chromosome, replicate), whether a map would be found — using
only ``stat`` calls, so it is safe to run on a head node and needs no ARG data.

By default it checks ONLY the filesystem candidates. ``resolve_mu_rate`` also
accepts an embedded ts.metadata ratemap first, and falls back to a scalar
``--mutation-rate`` last; pass ``--check-metadata`` to also open each tree
sequence (much slower, real I/O) and ``--mutation-rate R`` to model the config
fallback.

Usage:
    python scripts/check_mu_paths.py --root-dir amaranth
    python scripts/check_mu_paths.py --root-dir amaranth --tree-subdir trees
    python scripts/check_mu_paths.py --root-dir amaranth --mutation-rate 3e-8

Exit status is 1 if any pair would fail step 4, so it can gate a release.
"""
from __future__ import annotations

import argparse
import fnmatch
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from argtest_common import infer_mu_path  # noqa: E402

VALID_SUFFIXES = {".ts", ".trees", ".tsz"}


def natural_key(text: str):
    return [int(p) if p.isdigit() else p.lower() for p in re.split(r"(\d+)", text)]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--root-dir", required=True, type=Path,
                   help="Dataset root: one subdirectory per chromosome (config root_dir).")
    p.add_argument("--tree-subdir", default="",
                   help="Optional subdirectory inside each chromosome dir (config tree_subdir).")
    p.add_argument("--tree-pattern", default="*",
                   help="Filename glob (config tree_pattern; default '*').")
    p.add_argument("--mutation-rate", type=float, default=None,
                   help="Scalar fallback from config; if set, step 4 cannot fail on discovery.")
    p.add_argument("--check-metadata", action="store_true",
                   help="Also open each tree sequence to look for an embedded ratemap (slow).")
    p.add_argument("--per-chrom", action="store_true",
                   help="Report only the first replicate per chromosome (fast smoke test).")
    return p.parse_args()


def main():
    args = parse_args()
    root = args.root_dir.expanduser().resolve()
    if not root.is_dir():
        sys.exit(f"root_dir is not a directory: {root}")

    ratemap_from_metadata = load_ts = None
    if args.check_metadata:
        from argtest_common import load_ts, ratemap_from_metadata  # noqa: F811

    rows, failures = [], []
    chrom_dirs = sorted((d for d in root.iterdir() if d.is_dir()),
                        key=lambda d: natural_key(d.name))
    for chrom_dir in chrom_dirs:
        search = chrom_dir / args.tree_subdir if args.tree_subdir else chrom_dir
        if not search.is_dir():
            continue
        trees = sorted((p for p in search.iterdir()
                        if p.is_file() and p.suffix in VALID_SUFFIXES
                        and fnmatch.fnmatch(p.name, args.tree_pattern)),
                       key=lambda p: natural_key(p.name))
        if args.per_chrom:
            trees = trees[:1]
        for ts_path in trees:
            source, detail = None, ""
            if args.check_metadata:
                try:
                    if ratemap_from_metadata(load_ts(ts_path).metadata or {}) is not None:
                        source, detail = "metadata", "embedded ratemap"
                except Exception as exc:  # noqa: BLE001 - diagnostic tool
                    detail = f"metadata read failed: {exc}"
            if source is None:
                try:
                    source, detail = "sibling", str(infer_mu_path(ts_path))
                except FileNotFoundError as exc:
                    if args.mutation_rate is not None:
                        source, detail = "scalar", f"--mutation-rate {args.mutation_rate}"
                    else:
                        source = "MISSING"
                        detail = str(exc).split("Tried exact paths: ")[-1]
                        failures.append((ts_path, detail))
            rows.append((chrom_dir.name, ts_path.name, source, detail))

    if not rows:
        sys.exit(f"No tree files found under {root} (pattern={args.tree_pattern!r}, "
                 f"tree_subdir={args.tree_subdir!r})")

    width = max(len(f"{c}/{n}") for c, n, _, _ in rows)
    for chrom, name, source, detail in rows:
        mark = "FAIL" if source == "MISSING" else "ok"
        print(f"[{mark:>4}] {chrom + '/' + name:<{width}}  {source:<9} {detail}")

    counts: dict[str, int] = {}
    for _, _, source, _ in rows:
        counts[source] = counts.get(source, 0) + 1
    print(f"\n{len(rows)} tree file(s): " + ", ".join(f"{k}={v}" for k, v in sorted(counts.items())))

    if failures:
        print(f"\n{len(failures)} file(s) would FAIL step 4 under v1.10.")
        print("Fix by adding a map at one of the exact paths listed above, or by "
              "setting `mutation_rate` in the config.")
        sys.exit(1)
    print("\nAll files resolve a mutation rate; step 4 will not fail on discovery.")


if __name__ == "__main__":
    main()
