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
        "--log",
        type=Path,
        default=None,
        help="Optional log file path (default: <out-dir>/hapmap_low_rec_mask.log).",
    )
    args = p.parse_args()
    if not 0 < args.rec_fraction <= 1:
        raise ValueError("--rec-fraction must be in the interval (0, 1].")
    return args


def load_hapmap(path: Path):
    by_chr = defaultdict(list)
    with open(path, "r", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        required = {"Chromosome", "Position(bp)", "Rate(cM/Mb)"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Missing required HapMap columns: {sorted(missing)}")
        for row in reader:
            chrom = row["Chromosome"].strip()
            pos = int(float(row["Position(bp)"]))
            rate = float(row["Rate(cM/Mb)"])
            by_chr[chrom].append((pos, rate))
    return by_chr


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


def main():
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    hapmap = load_hapmap(args.hapmap)
    fai = load_fai(args.fai)
    total_written = 0
    per_chrom = {}
    for chrom, rows in sorted(hapmap.items()):
        if chrom not in fai:
            raise KeyError(f"Chromosome {chrom} not found in {args.fai}")
        intervals = build_intervals(rows, fai[chrom])
        out_path = args.out_dir / f"{chrom}.low_rec.mask.bed"
        if not intervals:
            out_path.write_text("")
            per_chrom[chrom] = 0
            continue

        n_keep = max(1, math.ceil(len(intervals) * args.rec_fraction))
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

    log_path = args.log or (args.out_dir / "hapmap_low_rec_mask.log")
    try:
        with open(log_path, "w") as fh:
            fh.write("# hapmap_low_rec_mask summary\n")
            for ln in summary_lines:
                fh.write(ln + "\n")
    except Exception:
        print(f"WARNING: failed to write log to {log_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
