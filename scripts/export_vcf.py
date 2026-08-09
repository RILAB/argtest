#!/usr/bin/env python3
"""Export a VCF from a (filtered) tree sequence.

Writes one VCF record per site present in the tree sequence. The pipeline does
not synthesise invariant sites, so the emitted records are the variable sites
carried on the ARG (some may be monomorphic among the *retained* samples after
trimming — those are still emitted).

Missing genotypes: ``trim_samples.py`` removes a sample's ancestry over an
interval, which leaves that sample **isolated** in the local trees there. tskit's
``write_vcf(..., isolated_as_missing=True)`` renders an isolated sample as a
missing genotype (``.``), so a sample pruned from particular trees is
automatically reported as missing at any site inside those intervals — without
dropping the site or the sample globally.

Ploidy follows each individual's own sample-node count: a haploid individual
(one sample node, e.g. the admix data) is written as a single allele ``0``/``1``
(not ``0/.``); a diploid individual as ``0/1``; etc.
"""

import argparse
import gzip
import sys
from pathlib import Path

import tskit

# When run as `python scripts/export_vcf.py`, the script's own directory is on
# sys.path[0], so this resolves the shared helpers (matches the other scripts).
from argtest_common import get_individual_name, load_ts


def parse_args():
    p = argparse.ArgumentParser(
        description="Export a VCF (variable sites; pruned samples as missing) from a tree sequence.",
    )
    p.add_argument("--ts", required=True, type=Path, help="Input tree sequence (.ts/.trees/.tsz).")
    p.add_argument("--out", required=True, type=Path, help="Output VCF path (.vcf or .vcf.gz).")
    p.add_argument(
        "--chrom",
        required=True,
        help="Contig name written to the VCF CHROM column and ##contig header.",
    )
    p.add_argument(
        "--name-substring-to-remove",
        default="",
        help='Substring removed globally from individual names before use as VCF sample names (default: "").',
    )
    p.add_argument("--log", type=Path, default=None, help="Optional summary log path.")
    return p.parse_args()


def resolve_samples(ts, name_substring_to_remove=""):
    """Decide how to group sample nodes into VCF columns.

    Returns ``(individuals, names, ploidy)``: a list of individual ids (or
    ``None``), the matching sample names, and the per-individual ploidy. When
    every sample node belongs to an individual and all such individuals carry
    the same number of sample nodes, group by individual (so haploid -> single
    allele, diploid -> a pair). Otherwise fall back to one haploid column per
    sample node, which always yields ``0``/``1``/``.`` genotypes.
    """
    sample_set = set(int(s) for s in ts.samples())
    individuals = []
    names = []
    ploidies = set()
    covered = 0
    for ind in ts.individuals():
        nodes = [n for n in ind.nodes if int(n) in sample_set]
        if not nodes:
            continue
        individuals.append(ind.id)
        names.append(
            get_individual_name(
                ind, name_substring_to_remove=name_substring_to_remove
            )
        )
        ploidies.add(len(nodes))
        covered += len(nodes)

    if individuals and covered == len(sample_set) and len(ploidies) == 1:
        return individuals, names, ploidies.pop()
    # Fallback: per-node haploid columns (covers no-individual or mixed-ploidy ts).
    return None, None, 1


def main():
    args = parse_args()
    ts = load_ts(args.ts)
    args.out.parent.mkdir(parents=True, exist_ok=True)

    individuals, names, ploidy = resolve_samples(
        ts, name_substring_to_remove=args.name_substring_to_remove
    )

    opener = gzip.open if str(args.out).endswith(".gz") else open
    with opener(args.out, "wt") as fh:
        if individuals is not None:
            ts.write_vcf(
                fh,
                contig_id=args.chrom,
                individuals=individuals,
                individual_names=names,
                isolated_as_missing=True,
            )
            n_cols = len(individuals)
        else:
            ts.write_vcf(
                fh,
                contig_id=args.chrom,
                ploidy=1,
                isolated_as_missing=True,
            )
            n_cols = ts.num_samples

    summary = [
        f"Wrote VCF: {args.out}",
        f"chrom={args.chrom}",
        f"sites={ts.num_sites}",
        f"sample_columns={n_cols} (ploidy={ploidy})",
        f"sample_nodes={ts.num_samples}",
    ]
    for ln in summary:
        print(ln)
        print(ln, file=sys.stderr)

    if args.log is not None:
        args.log.parent.mkdir(parents=True, exist_ok=True)
        with open(args.log, "w") as fh:
            fh.write("# export_vcf summary\n")
            fh.write(f"ts={args.ts}\n")
            for ln in summary:
                fh.write(ln + "\n")


if __name__ == "__main__":
    main()
