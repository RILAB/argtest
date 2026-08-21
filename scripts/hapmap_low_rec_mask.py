#!/usr/bin/env python
from __future__ import annotations

import argparse
import csv
import math
import re
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


def load_hapmap(path: Path):
    by_chr = defaultdict(list)
    with open(path, "r", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        for row in reader:
            if not row or row[0].strip() == "Chromosome":
                continue
            chrom = row[0].strip()
            pos = int(float(row[1]))
            rate = float(row[2])
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


def _strip_prefix(chrom):
    """Strip base-name prefix from pipeline chrom names (e.g. 'amaranth.1' -> '1')."""
    return chrom.split(".", 1)[1] if "." in chrom else chrom


# 'chr', 'chr_' and 'chr-' are all in use as chromosome-name prefixes.
_CHR_PREFIX_RE = re.compile(r"^chr[_-]?", re.IGNORECASE)


def _canonical_chrom(name):
    """Drop a 'chr'/'chr_'/'chr-' prefix so the two naming conventions compare equal.

    Deliberately does NOT strip a base-name prefix: that is applied to the
    QUERY only (see _resolve_chrom). Applying it to the available names as well
    would reduce an unrelated contig like 'chrUn_random.1' to '1' and collide it
    with 'chr1'.
    """
    return _CHR_PREFIX_RE.sub("", name).casefold()


def _resolve_chrom(chrom, available, source=None):
    """Resolve a pipeline chrom name against available names (a dict or set).

    Matching is symmetric in the 'chr' prefix: '1', 'chr1', 'chr_1' and 'Chr-1'
    all resolve to each other in either direction, so a bare-numeric map matches
    'chr'-prefixed pipeline chroms just as a 'chr'-prefixed map matches bare
    ones. Pipeline names may additionally carry a base-name prefix
    ('amaranth.1', 'combined.1'), which is stripped from the query.

    An exact match always wins, so a file that uses one convention consistently
    is never reinterpreted. Returns the matching key, or None if nothing matches.

    Raises ValueError if two different available names resolve to the same
    chromosome (e.g. a map holding both '1' and 'chr1'): picking one silently
    would attach the wrong recombination map or the wrong sequence length, and
    which one you got would depend on iteration order. Pass `source` (a path or
    description) to name the offending file in that error.
    """
    if chrom in available:
        return chrom
    for query in (chrom, _strip_prefix(chrom)):
        target = _canonical_chrom(query)
        matches = sorted({key for key in available if _canonical_chrom(key) == target})
        if len(matches) > 1:
            where = f" in {source}" if source else ""
            raise ValueError(
                f"Chromosome {chrom!r} matches more than one entry{where} "
                f"({', '.join(repr(m) for m in matches)}). Rename them so only one "
                f"spelling of each chromosome is present."
            )
        if matches:
            return matches[0]
    return None


def main():
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    fai = load_fai(args.fai)

    # Load once: chromosome resolution and interval construction use the same
    # parsed map instead of scanning a large HapMap file twice per rule.
    hapmap = load_hapmap(args.hapmap)
    if args.chrom is not None:
        hapmap_key = _resolve_chrom(args.chrom, hapmap, source=args.hapmap)
        if hapmap_key is None:
            raise KeyError(f"Chromosome {args.chrom!r} not found in {args.hapmap}")
        hapmap = {args.chrom: hapmap[hapmap_key]}

    total_written = 0
    per_chrom = {}
    for chrom, rows in sorted(hapmap.items()):
        fai_chrom = _resolve_chrom(chrom, fai, source=args.fai)
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
