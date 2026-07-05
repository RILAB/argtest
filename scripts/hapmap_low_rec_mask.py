#!/usr/bin/env python
from __future__ import annotations

import argparse
import csv
import math
import sys
from collections import defaultdict
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Extract the bottom fraction of recombination-rate intervals from a HapMap file, "
            "separately for each chromosome."
        )
    )
    p.add_argument("--hapmap", required=True, type=Path, help="Input HapMap file.")
    p.add_argument("--fai", required=True, type=Path, help="FASTA index (.fai) with chromosome lengths.")
    p.add_argument(
        "--rec-fraction",
        required=True,
        type=float,
        help="Fraction of intervals to keep from the low end of Rate(cM/Mb), per chromosome.",
    )
    p.add_argument(
        "--out-dir",
        type=Path,
        default=Path("."),
        help="Directory for <chr>.low_rec.mask.bed outputs (default: current directory).",
    )
    p.add_argument(
        "--chrom",
        default=None,
        help="If provided, process only this chromosome.",
    )
    p.add_argument(
        "--log",
        type=Path,
        default=None,
        help="Optional log file path (default: <out-dir>/logs/hapmap_low_rec_mask.log).",
    )
    args = p.parse_args()
    if not 0 <= args.rec_fraction <= 1:
        raise ValueError("--rec-fraction must be in the interval [0, 1].")
    return args


def load_hapmap(path: Path, chrom_filter: str | None = None, output_chrom: str | None = None):
    by_chr = defaultdict(list)
    with open(path, "r", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        for row in reader:
            if not row or row[0].strip() == "Chromosome":
                continue
            chrom = row[0].strip()
            if chrom_filter is not None and chrom != chrom_filter:
                continue
            pos = int(float(row[1]))
            rate = float(row[2])
            if output_chrom is not None:
                chrom = output_chrom
            by_chr[chrom].append((pos, rate))
    return by_chr


def hapmap_chroms(path: Path) -> set[str]:
    chroms = set()
    with open(path, "r", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        for row in reader:
            if not row or row[0].strip() == "Chromosome":
                continue
            chroms.add(row[0].strip())
    return chroms


def load_fai(path: Path):
    lengths = {}
    with open(path, "r", newline="") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                raise ValueError(f"Invalid FAI line in {path}: {line}")
            lengths[parts[0]] = int(parts[1])
    return lengths


def build_intervals(rows, chrom_length):
    rows = sorted(rows)
    intervals = []
    if not rows:
        return intervals

    first_pos, first_rate = rows[0]
    if first_pos > 0:
        intervals.append((0, first_pos, first_rate))

    for i in range(len(rows)):
        start, rate = rows[i]
        if i + 1 < len(rows):
            end = rows[i + 1][0]
        else:
            end = chrom_length
        if end > start:
            intervals.append((start, end, rate))
    return intervals


def _strip_prefix(chrom):
    """Strip base-name prefix from pipeline chrom names (e.g. 'amaranth.1' -> '1')."""
    return chrom.split(".", 1)[1] if "." in chrom else chrom


def _resolve_chrom(chrom, available):
    """Resolve a pipeline chrom name against available names (a dict or set).

    Tries the full name, the suffix after the first dot, then common 'chr'-
    prefixed variants of that suffix. Returns the matching key or None. Used for
    both the hapmap and the FAI so a 'chr'-prefixed map (e.g. 'chr1') still
    matches pipeline chrom names like 'combined.1'.
    """
    if chrom in available:
        return chrom
    bare = _strip_prefix(chrom)
    for candidate in (bare, f"chr_{bare}", f"chr{bare}"):
        if candidate in available:
            return candidate
    return None


def main():
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    fai = load_fai(args.fai)

    if args.chrom is not None:
        hapmap_key = _resolve_chrom(args.chrom, hapmap_chroms(args.hapmap))
        if hapmap_key is None:
            raise KeyError(f"Chromosome {args.chrom!r} not found in {args.hapmap}")
        hapmap = load_hapmap(args.hapmap, chrom_filter=hapmap_key, output_chrom=args.chrom)
    else:
        hapmap = load_hapmap(args.hapmap)

    total_written = 0
    per_chrom = {}
    for chrom, rows in sorted(hapmap.items()):
        fai_chrom = _resolve_chrom(chrom, fai)
        if fai_chrom is None:
            raise KeyError(f"Chromosome {chrom!r} not found in {args.fai}")
        intervals = build_intervals(rows, fai[fai_chrom])
        out_path = args.out_dir / f"{chrom}.low_rec.mask.bed"
        if not intervals:
            out_path.write_text("")
            per_chrom[chrom] = 0
            continue

        n_keep = math.ceil(len(intervals) * args.rec_fraction)
        keep = set(sorted(range(len(intervals)), key=lambda i: intervals[i][2])[:n_keep])
        lines = [
            f"{chrom}\t{start}\t{end}\t{rate:.6g}"
            for i, (start, end, rate) in enumerate(intervals)
            if i in keep
        ]
        out_path.write_text("\n".join(lines) + ("\n" if lines else ""))
        per_chrom[chrom] = len(lines)
        total_written += len(lines)

    # Summary: print to stdout and stderr and write to log
    summary_lines = [f"Wrote low-rec BEDs to: {args.out_dir}", f"Chromosomes processed: {len(per_chrom)}", f"Total intervals written: {total_written}"]
    for chrom, n in sorted(per_chrom.items()):
        summary_lines.append(f"  {chrom}: {n} intervals")
    for ln in summary_lines:
        print(ln)
        print(ln, file=sys.stderr)

    log_path = args.log or (args.out_dir / "logs" / "hapmap_low_rec_mask.log")
    try:
        log_path.parent.mkdir(parents=True, exist_ok=True)
        with open(log_path, "w") as fh:
            fh.write("# hapmap_low_rec_mask summary\n")
            for ln in summary_lines:
                fh.write(ln + "\n")
    except Exception:
        print(f"WARNING: failed to write log to {log_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
