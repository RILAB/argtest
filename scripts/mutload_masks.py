#!/usr/bin/env python
from __future__ import annotations

import argparse
from pathlib import Path

import msprime
import numpy as np

from argtest_common import (
    aggregate_by_individual,
    load_ts,
    mutational_load,
    resolve_mu_rate,
    sample_names,
)


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Write mutational-load outlier BEDs for one tree sequence. "
            "Outputs are chromosome-scoped to avoid replicate-name collisions."
        )
    )
    parser.add_argument("--ts", required=True, type=Path, help="Input tree sequence file.")
    parser.add_argument("--chrom", required=True, help="Chromosome label written in BED column 1.")
    window_group = parser.add_mutually_exclusive_group(required=True)
    window_group.add_argument("--window-size", type=float, help="Window size in bp.")
    window_group.add_argument("--snp-window", type=int, help="Window size in number of variants.")
    parser.add_argument(
        "--cutoff",
        type=float,
        default=0.5,
        help=(
            "Outlier cutoff as a fraction of the per-individual sim-based expected "
            "load (default: 0.5)."
        ),
    )
    parser.add_argument(
        "--fraction",
        type=float,
        default=None,
        help=(
            "If provided, windows with outlier fraction > this value are written to --masked-bed "
            "and excluded from --outlier-bed."
        ),
    )
    parser.add_argument(
        "--mutation-rate",
        type=float,
        default=None,
        help=(
            "Scalar mutation-rate fallback for the sim-based expected load when neither "
            "ts.metadata nor a sibling *.mut_rate.p file provides a ratemap."
        ),
    )
    parser.add_argument(
        "--random-seed",
        type=int,
        default=1,
        help="Seed for msprime.sim_mutations when computing expected load (default: 1).",
    )
    parser.add_argument(
        "--suffix-to-strip",
        default="",
        help="Suffix removed from sample names before grouping to individuals.",
    )
    parser.add_argument("--outlier-bed", required=True, type=Path, help="Output outlier BED path.")
    parser.add_argument(
        "--masked-bed",
        required=True,
        type=Path,
        help=(
            "Output mutation-masked BED path. "
            "If --fraction is omitted, this file is created as empty."
        ),
    )
    parser.add_argument(
        "--log",
        type=Path,
        default=None,
        help=(
            "Optional log file path (default: <outlier-bed.parent>/logs/"
            "<chrom>.<ts_stem>.mutload.log)."
        ),
    )
    args = parser.parse_args()
    if args.window_size is not None and args.window_size <= 0:
        raise ValueError("--window-size must be > 0")
    if args.snp_window is not None and args.snp_window <= 0:
        raise ValueError("--snp-window must be > 0")
    if args.fraction is not None and not 0 <= args.fraction <= 1:
        raise ValueError("--fraction must be between 0 and 1")
    return args


def build_bp_windows(ts, window_size: float) -> np.ndarray:
    sequence_length = float(ts.sequence_length)
    windows = np.arange(0, sequence_length + window_size, window_size, dtype=float)
    if windows[-1] > sequence_length:
        windows[-1] = sequence_length
    return windows


def build_snp_windows(ts, snp_window: int) -> np.ndarray:
    positions = np.asarray(ts.sites_position, dtype=float)
    sequence_length = float(ts.sequence_length)
    if positions.size == 0:
        return np.array([0.0, sequence_length], dtype=float)
    edges = positions[snp_window::snp_window]
    return np.concatenate((np.array([0.0]), edges, np.array([sequence_length])))


def simulate_expected_load(ts, windows, names, mu, seed):
    # Single-sim expected per-individual load matrix (windows x individuals).
    sim_ts = msprime.sim_mutations(ts, rate=mu, keep=False, random_seed=seed)
    expected = mutational_load(sim_ts, windows=windows)
    expected, unique = aggregate_by_individual(expected, names)
    return expected, unique


def main():
    args = parse_args()
    args.outlier_bed.parent.mkdir(parents=True, exist_ok=True)
    args.masked_bed.parent.mkdir(parents=True, exist_ok=True)
    log_path = args.log or (
        args.outlier_bed.parent / "logs" / f"{args.chrom}.{args.ts.stem}.mutload.log"
    )
    log_path.parent.mkdir(parents=True, exist_ok=True)

    ts = load_ts(args.ts)
    if args.window_size is not None:
        windows = build_bp_windows(ts, args.window_size)
    else:
        windows = build_snp_windows(ts, args.snp_window)

    load = mutational_load(ts, windows=windows)
    names = sample_names(ts, suffix_to_strip=args.suffix_to_strip)
    load, unique_names = aggregate_by_individual(load, names)

    mu = resolve_mu_rate(ts, args.ts, scalar_fallback=args.mutation_rate)
    expected, expected_names = simulate_expected_load(
        ts, windows, names, mu=mu, seed=args.random_seed
    )
    if expected_names != unique_names:
        raise RuntimeError("Sim and observed individual orderings disagree")

    high = (1 + args.cutoff) * expected
    low = (1 - args.cutoff) * expected
    outlier_mask = (load > high) | (load < low)

    masked_window_mask = np.zeros(load.shape[0], dtype=bool)
    masked_lines = []
    if args.fraction is not None:
        # A window is "masked" if too many individuals look like outliers.
        outlier_fractions = outlier_mask.sum(axis=1) / load.shape[1]
        masked_window_mask = outlier_fractions > args.fraction
        for w in range(load.shape[0]):
            if not masked_window_mask[w]:
                continue
            start = int(windows[w])
            end = int(windows[w + 1])
            masked_lines.append(
                f"{args.chrom}\t{start}\t{end}\t{outlier_fractions[w]:.3f}\t{int(outlier_mask[w].sum())}\t{load.shape[1]}"
            )
    args.masked_bed.write_text("\n".join(masked_lines) + ("\n" if masked_lines else ""))

    outlier_lines = []
    for w in range(load.shape[0]):
        if masked_window_mask[w]:
            continue
        row_mask = outlier_mask[w]
        if not row_mask.any():
            continue
        idx = np.where(row_mask)[0]
        outlier_names = [unique_names[i] for i in idx]
        outlier_vals = [f"{load[w, i]:.3f}" for i in idx]
        outlier_exp = [f"{expected[w, i]:.3f}" for i in idx]
        start = int(windows[w])
        end = int(windows[w + 1])
        outlier_lines.append(
            f"{args.chrom}\t{start}\t{end}\t{','.join(outlier_names)}\t{','.join(outlier_vals)}\t{','.join(outlier_exp)}"
        )
    args.outlier_bed.write_text("\n".join(outlier_lines) + ("\n" if outlier_lines else ""))

    with open(log_path, "w") as fh:
        fh.write("# mutload_masks summary\n")
        fh.write(f"ts={args.ts}\n")
        fh.write(f"chrom={args.chrom}\n")
        fh.write(f"outlier_bed={args.outlier_bed}\n")
        fh.write(f"masked_bed={args.masked_bed}\n")
        fh.write(f"windows={len(windows) - 1}\n")
        fh.write(f"random_seed={args.random_seed}\n")
        fh.write(f"cutoff={args.cutoff}\n")
        fh.write(f"outliers_written={len(outlier_lines)}\n")
        fh.write(f"masked_windows={len(masked_lines)}\n")

    print(f"Wrote outlier BED: {args.outlier_bed}")
    print(f"Wrote masked BED: {args.masked_bed}")
    print(f"Wrote log: {log_path}")


if __name__ == "__main__":
    main()
