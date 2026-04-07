#!/usr/bin/env python
from __future__ import annotations

import argparse
import math
import pickle
from collections import defaultdict
from pathlib import Path

import numpy as np

from argtest_common import (
    aggregate_by_individual,
    collapse_masked_intervals,
    dump_ts,
    load_remove_intervals,
    load_ts,
    merge_intervals,
    mutational_load,
    name_to_nodes_map,
    ratemap_from_keep_intervals,
    sample_names,
    validate_trimmed_ts,
)
from find_low_access_regions import infer_mu_path
from hapmap_low_rec_mask import build_intervals, load_fai, load_hapmap
from merge_treefiles_by_replicate import VALID_SUFFIXES, natural_key
from trim_regions import complement_intervals, load_mask_intervals
from trim_samples import remove_ancestry


def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Run steps 1-5 across chromosome subdirectories and concatenate per replicate. "
            "Expected input layout: <root>/<chromosome>/<replicate><suffix>."
        )
    )
    p.add_argument(
        "--root",
        required=True,
        type=Path,
        help="Root directory containing one subdirectory per chromosome.",
    )
    p.add_argument("--hapmap", required=True, type=Path, help="HapMap file for step 1.")
    p.add_argument("--fai", required=True, type=Path, help="FAI file for step 1.")
    p.add_argument(
        "--rec-fraction",
        required=True,
        type=float,
        help="Step-1 low-recombination fraction in (0, 1].",
    )
    p.add_argument(
        "--window-size",
        required=True,
        type=float,
        help="Step-2 window size (bp) for low-access regions.",
    )
    p.add_argument(
        "--cutoff-bp",
        required=True,
        type=float,
        help="Step-2 accessibility cutoff in bp.",
    )
    mutload_group = p.add_mutually_exclusive_group(required=True)
    mutload_group.add_argument(
        "--mutload-window-size",
        type=float,
        help="Step-3 mutational-load window size in bp.",
    )
    mutload_group.add_argument(
        "--mutload-snp-window",
        type=int,
        help="Step-3 mutational-load window size in variants.",
    )
    p.add_argument(
        "--mutload-cutoff",
        type=float,
        default=0.25,
        help="Step-3 outlier cutoff as a fraction of the window median (default: 0.25).",
    )
    p.add_argument(
        "--mutload-fraction",
        type=float,
        default=None,
        help=(
            "Optional step-3 threshold: write mutation-masked BED rows where outlier fraction "
            "is greater than this value."
        ),
    )
    p.add_argument(
        "--pattern",
        default="*",
        help="Glob pattern for tree files inside each chromosome directory (default: '*').",
    )
    p.add_argument(
        "--suffix-to-strip",
        default="_anchorwave",
        help="Suffix stripped from individual IDs when parsing names (default: _anchorwave).",
    )
    p.add_argument(
        "--out-dir",
        type=Path,
        default=None,
        help="Output directory (default: <root>/batch_steps1_5).",
    )
    p.add_argument(
        "--merged-dir",
        type=Path,
        default=None,
        help="Merged per-replicate output directory (default: <out-dir>/combined).",
    )
    p.add_argument(
        "--base-name",
        default=None,
        help="Base name used for merged outputs (default: root directory name).",
    )
    p.add_argument(
        "--out-suffix",
        choices=sorted(VALID_SUFFIXES),
        default=None,
        help="Optional merged output suffix (default: suffix of the first chromosome file).",
    )
    p.add_argument(
        "--allow-missing-replicates",
        action="store_true",
        help="Allow concatenation when a replicate is missing in some chromosome directories.",
    )
    args = p.parse_args()
    if not 0 < args.rec_fraction <= 1:
        raise ValueError("--rec-fraction must be in the interval (0, 1].")
    if args.window_size <= 0:
        raise ValueError("--window-size must be > 0.")
    if args.cutoff_bp < 0:
        raise ValueError("--cutoff-bp must be >= 0.")
    if args.mutload_window_size is not None and args.mutload_window_size <= 0:
        raise ValueError("--mutload-window-size must be > 0.")
    if args.mutload_snp_window is not None and args.mutload_snp_window <= 0:
        raise ValueError("--mutload-snp-window must be > 0.")
    if args.mutload_fraction is not None and not 0 <= args.mutload_fraction <= 1:
        raise ValueError("--mutload-fraction must be between 0 and 1.")
    return args


def discover_chromosome_dirs(root: Path, pattern: str) -> dict[str, list[Path]]:
    if not root.exists():
        raise FileNotFoundError(f"Root directory does not exist: {root}")
    if not root.is_dir():
        raise NotADirectoryError(f"Expected a directory for --root: {root}")
    chrom_to_files = {}
    for d in sorted([p for p in root.iterdir() if p.is_dir()], key=lambda p: natural_key(p.name)):
        files = sorted(
            [p for p in d.glob(pattern) if p.is_file() and p.suffix in VALID_SUFFIXES],
            key=lambda p: natural_key(p.name),
        )
        if files:
            chrom_to_files[d.name] = files
    if not chrom_to_files:
        raise RuntimeError(
            f"No chromosome subdirectories with tree files found in {root} using pattern '{pattern}'."
        )
    return chrom_to_files


def write_step1_low_rec(args, chrom_names: set[str], out_dir: Path) -> dict[str, Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    hapmap = load_hapmap(args.hapmap)
    fai = load_fai(args.fai)
    out = {}
    for chrom in sorted(chrom_names, key=natural_key):
        if chrom not in hapmap:
            raise KeyError(f"Chromosome '{chrom}' not found in HapMap file: {args.hapmap}")
        if chrom not in fai:
            raise KeyError(f"Chromosome '{chrom}' not found in FAI file: {args.fai}")
        intervals = build_intervals(hapmap[chrom], fai[chrom])
        out_path = out_dir / f"{chrom}.low_rec.mask.bed"
        if intervals:
            n_keep = max(1, math.ceil(len(intervals) * args.rec_fraction))
            keep = set(sorted(range(len(intervals)), key=lambda i: intervals[i][2])[:n_keep])
            lines = [
                f"{chrom}\t{start}\t{end}\t{rate:.6g}"
                for i, (start, end, rate) in enumerate(intervals)
                if i in keep
            ]
            out_path.write_text("\n".join(lines) + "\n")
        else:
            out_path.write_text("")
        out[chrom] = out_path
    return out


def write_step2_low_access(chrom: str, chrom_files: list[Path], window_size: float, cutoff_bp: float, out_path: Path):
    out_path.parent.mkdir(parents=True, exist_ok=True)
    first_ts_path = chrom_files[0]
    ts = load_ts(first_ts_path)
    mu_path = infer_mu_path(first_ts_path)
    with open(mu_path, "rb") as fh:
        mu = pickle.load(fh)
    from argtest_common import accessible_intervals_from_mu, overlap_lengths

    sequence_length = float(ts.sequence_length)
    windows = np.arange(0, sequence_length + window_size, window_size, dtype=float)
    if windows[-1] > sequence_length:
        windows[-1] = sequence_length
    acc_bp = overlap_lengths(accessible_intervals_from_mu(mu), windows)
    lines = []
    for i in range(len(windows) - 1):
        if acc_bp[i] < cutoff_bp:
            lines.append(f"{chrom}\t{int(windows[i])}\t{int(windows[i + 1])}\t{float(acc_bp[i]):.3f}")
    out_path.write_text("\n".join(lines) + ("\n" if lines else ""))


def build_mutload_windows(ts, window_size: float | None, snp_window: int | None):
    if window_size is not None:
        windows = np.arange(0, float(ts.sequence_length) + window_size, window_size, dtype=float)
        if windows[-1] > float(ts.sequence_length):
            windows[-1] = float(ts.sequence_length)
        return windows
    positions = np.asarray(ts.sites_position, dtype=float)
    if positions.size == 0:
        return np.array([0.0, float(ts.sequence_length)], dtype=float)
    edges = positions[snp_window::snp_window]
    return np.concatenate((np.array([0.0]), edges, np.array([float(ts.sequence_length)])))


def write_step3_mutload_masks(
    ts_path: Path,
    chrom: str,
    outlier_path: Path,
    masked_path: Path,
    window_size: float | None,
    snp_window: int | None,
    cutoff: float,
    fraction: float | None,
    suffix_to_strip: str,
):
    outlier_path.parent.mkdir(parents=True, exist_ok=True)
    if masked_path is not None:
        masked_path.parent.mkdir(parents=True, exist_ok=True)

    ts = load_ts(ts_path)
    windows = build_mutload_windows(ts, window_size=window_size, snp_window=snp_window)
    load = mutational_load(ts, windows=windows)
    names = sample_names(ts, suffix_to_strip=suffix_to_strip)
    load, unique_names = aggregate_by_individual(load, names)

    window_medians = np.median(load, axis=1)
    valid = window_medians > 0
    high = (1 + cutoff) * window_medians
    low = (1 - cutoff) * window_medians
    mask = ((load > high[:, None]) | (load < low[:, None])) & valid[:, None]

    masked_window_mask = np.zeros(load.shape[0], dtype=bool)
    if fraction is not None and masked_path is not None:
        outlier_fractions = mask.sum(axis=1) / load.shape[1]
        masked_window_mask = valid & (outlier_fractions > fraction)
        masked_lines = []
        for w in range(load.shape[0]):
            if not masked_window_mask[w]:
                continue
            start = int(windows[w])
            end = int(windows[w + 1])
            masked_lines.append(
                f"{chrom}\t{start}\t{end}\t{outlier_fractions[w]:.3f}\t{int(mask[w].sum())}\t{load.shape[1]}"
            )
        masked_path.write_text("\n".join(masked_lines) + ("\n" if masked_lines else ""))

    outlier_lines = []
    for w in range(load.shape[0]):
        if not valid[w] or masked_window_mask[w]:
            continue
        row_mask = mask[w]
        if not row_mask.any():
            continue
        outlier_names = [unique_names[i] for i in np.where(row_mask)[0]]
        outlier_vals = [f"{load[w, i]:.3f}" for i in np.where(row_mask)[0]]
        start = int(windows[w])
        end = int(windows[w + 1])
        outlier_lines.append(
            f"{chrom}\t{start}\t{end}\t{','.join(outlier_names)}\t{','.join(outlier_vals)}\t{window_medians[w]:.3f}"
        )
    outlier_path.write_text("\n".join(outlier_lines) + ("\n" if outlier_lines else ""))


def read_bed_intervals(path: Path):
    intervals = []
    if not path.exists():
        return intervals
    with open(path, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            start = float(parts[1])
            end = float(parts[2])
            if end > start:
                intervals.append([start, end])
    return intervals


def write_combined_remove_bed(chrom: str, source_paths: list[Path], out_path: Path):
    all_intervals = []
    for path in source_paths:
        all_intervals.extend(read_bed_intervals(path))
    merged = merge_intervals(all_intervals)
    lines = [f"{chrom}\t{int(left)}\t{int(right)}" for left, right in merged]
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(lines) + ("\n" if lines else ""))


def run_step4_trim_regions(ts_path: Path, remove_bed: Path, out_path: Path):
    ts = load_ts(ts_path)
    masked = load_mask_intervals(remove_bed, ts.sequence_length)
    keep = complement_intervals(masked, ts.sequence_length)
    accessible = ratemap_from_keep_intervals(keep, ts.sequence_length)
    trimmed = collapse_masked_intervals(ts, accessible)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    dump_ts(trimmed, out_path)


def bed_has_records(path: Path) -> bool:
    if not path.exists():
        return False
    with open(path, "r") as fh:
        for line in fh:
            line = line.strip()
            if line and not line.startswith("#"):
                return True
    return False


def run_step5_trim_samples(ts_path: Path, remove_bed: Path, out_path: Path, suffix_to_strip: str):
    ts = load_ts(ts_path)
    if not bed_has_records(remove_bed):
        out_path.parent.mkdir(parents=True, exist_ok=True)
        dump_ts(ts, out_path)
        return

    remove_intervals = load_remove_intervals([str(remove_bed)])
    if not remove_intervals:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        dump_ts(ts, out_path)
        return

    trimmed_ts = ts
    for name, intervals in remove_intervals.items():
        node_map = name_to_nodes_map(trimmed_ts, suffix_to_strip=suffix_to_strip)
        samples = node_map.get(name, [])
        if not samples:
            continue
        for left, right in zip(intervals["starts"], intervals["ends"]):
            trimmed_ts = remove_ancestry(trimmed_ts, samples, left, right)
    validate_trimmed_ts(trimmed_ts)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    dump_ts(trimmed_ts, out_path)


def concatenate_by_replicate(
    replicate_map: dict[str, list[tuple[str, Path]]],
    out_dir: Path,
    base_name: str,
    out_suffix: str | None,
    n_chromosomes: int,
    allow_missing: bool,
):
    out_dir.mkdir(parents=True, exist_ok=True)
    written = []
    for replicate, chrom_paths in sorted(replicate_map.items(), key=lambda kv: natural_key(kv[0])):
        if (not allow_missing) and len(chrom_paths) != n_chromosomes:
            missing = n_chromosomes - len(chrom_paths)
            raise RuntimeError(
                f"Replicate '{replicate}' is missing {missing} chromosome file(s); "
                "use --allow-missing-replicates to concatenate partial sets."
            )
        ordered = sorted(chrom_paths, key=lambda item: natural_key(item[0]))
        merged = load_ts(ordered[0][1])
        for _, path in ordered[1:]:
            merged = merged.concatenate(load_ts(path))
        suffix = out_suffix if out_suffix is not None else ordered[0][1].suffix
        out_path = out_dir / f"{base_name}.combined.{replicate}{suffix}"
        dump_ts(merged, out_path)
        written.append(out_path)
    return written


def main():
    args = parse_args()
    out_dir = args.out_dir if args.out_dir is not None else (args.root / "batch_steps1_5")
    merged_dir = args.merged_dir if args.merged_dir is not None else (out_dir / "combined")
    base_name = args.base_name or args.root.name

    chrom_to_files = discover_chromosome_dirs(args.root, args.pattern)
    chrom_names = set(chrom_to_files.keys())

    step1_dir = out_dir / "step1_low_rec"
    step2_dir = out_dir / "step2_low_access"
    step3_dir = out_dir / "step3_mutload"
    step4_mask_dir = out_dir / "step4_masks"
    step4_trimmed_dir = out_dir / "step4_trimmed_regions"
    step5_trimmed_dir = out_dir / "step5_trimmed_samples"
    summary_log = out_dir / "logs" / "run_steps1_5_and_concat.log"

    step1_masks = write_step1_low_rec(args, chrom_names=chrom_names, out_dir=step1_dir)
    replicate_outputs = defaultdict(list)

    for chrom in sorted(chrom_to_files.keys(), key=natural_key):
        chrom_files = chrom_to_files[chrom]
        low_access_bed = step2_dir / chrom / f"{chrom}.low_access.bed"
        write_step2_low_access(
            chrom=chrom,
            chrom_files=chrom_files,
            window_size=args.window_size,
            cutoff_bp=args.cutoff_bp,
            out_path=low_access_bed,
        )

        for ts_path in chrom_files:
            replicate = ts_path.stem
            outlier_bed = step3_dir / chrom / f"{replicate}.outliers.bed"
            masked_bed = step3_dir / chrom / f"{replicate}.mutation_masked.bed"
            write_step3_mutload_masks(
                ts_path=ts_path,
                chrom=chrom,
                outlier_path=outlier_bed,
                masked_path=masked_bed,
                window_size=args.mutload_window_size,
                snp_window=args.mutload_snp_window,
                cutoff=args.mutload_cutoff,
                fraction=args.mutload_fraction,
                suffix_to_strip=args.suffix_to_strip,
            )

            combined_bed = step4_mask_dir / chrom / f"{replicate}.remove_regions.bed"
            bed_sources = [step1_masks[chrom], low_access_bed]
            if args.mutload_fraction is not None:
                bed_sources.append(masked_bed)
            write_combined_remove_bed(chrom=chrom, source_paths=bed_sources, out_path=combined_bed)

            step4_out = step4_trimmed_dir / chrom / f"{replicate}{ts_path.suffix}"
            run_step4_trim_regions(ts_path=ts_path, remove_bed=combined_bed, out_path=step4_out)

            step5_out = step5_trimmed_dir / chrom / f"{replicate}{ts_path.suffix}"
            run_step5_trim_samples(
                ts_path=step4_out,
                remove_bed=outlier_bed,
                out_path=step5_out,
                suffix_to_strip=args.suffix_to_strip,
            )
            replicate_outputs[replicate].append((chrom, step5_out))

    merged = concatenate_by_replicate(
        replicate_map=replicate_outputs,
        out_dir=merged_dir,
        base_name=base_name,
        out_suffix=args.out_suffix,
        n_chromosomes=len(chrom_to_files),
        allow_missing=args.allow_missing_replicates,
    )

    summary_log.parent.mkdir(parents=True, exist_ok=True)
    with open(summary_log, "w") as fh:
        fh.write("# run_steps1_5_and_concat summary\n")
        fh.write(f"root={args.root}\n")
        fh.write(f"chromosomes={len(chrom_to_files)}\n")
        fh.write(f"replicates={len(replicate_outputs)}\n")
        fh.write(f"merged_outputs={len(merged)}\n")
        for path in merged:
            fh.write(f"merged={path}\n")

    print(f"Finished. Wrote merged outputs to: {merged_dir}")
    print(f"Summary log: {summary_log}")


if __name__ == "__main__":
    main()
