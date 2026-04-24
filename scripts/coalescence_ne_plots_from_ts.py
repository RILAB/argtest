#!/usr/bin/env python
"""Generate pair-coalescence and Ne plots from tree-sequence replicates.

WARNING: tskit computes coalescence rates under the assumption that every
tree spans the full sample set. On ARGs with partial trees — e.g. after
`trim_regions.py` / `trim_samples.py` have removed intervals or individuals,
or any tree sequence with locally reduced sample membership — the rate is
not corrected for the per-tree effective sample count and the resulting Ne
estimate is biased. Treat Ne output from such tree sequences with caution
until the rate calculation is patched (see project TODO list).
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path
import warnings

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp")

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from argtest_common import load_ts

matplotlib.rcParams["figure.dpi"] = 300


def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Generate pair coalescence and Ne plots from a set of tree "
            "sequence replicates. WARNING: tskit's coalescence-rate "
            "calculation does not correct for partial trees (trees that "
            "cover only a subset of samples after region/sample trimming), "
            "so Ne estimates on trimmed ARGs are biased — see the module "
            "docstring for details."
        ),
    )
    p.add_argument(
        "--ts-dir",
        required=True,
        type=Path,
        help="Directory containing tree sequence files (.tsz, .ts, .trees).",
    )
    p.add_argument(
        "--pattern",
        default="*.tsz",
        help="Glob pattern for input trees (default: *.tsz).",
    )
    p.add_argument(
        "--out-dir",
        type=Path,
        default=Path("results/coalescence_ne_plots"),
        help="Output directory for plots (default: results/coalescence_ne_plots).",
    )
    p.add_argument(
        "--time-bins-file",
        required=True,
        type=Path,
        help=(
            "File of explicit time-bin edges (one per line) defining the "
            "coalescence time grid."
        ),
    )
    p.add_argument(
        "--burnin-frac",
        type=float,
        default=0.5,
        help=(
            "Fraction of initial MCMC trees to discard as burn-in when "
            "computing posterior means (default: 0.5)."
        ),
    )
    p.add_argument(
        "--tail-cutoff",
        type=float,
        default=1e-12,
        help=(
            "Minimum probability-mass threshold for pair-coalescence tail "
            "trimming (default: 1e-12)."
        ),
    )
    p.add_argument(
        "--time-adjust",
        type=float,
        default=1.0,
        help=(
            "Divide plotted time-axis values by this factor (e.g. generations "
            "per year) to convert to calendar time (default: 1.0)."
        ),
    )
    p.add_argument(
        "--log-rates",
        action="store_true",
        help="Use log y-axis for pair-coalescence-rates and Ne plots.",
    )
    p.add_argument(
        "--year",
        type=float,
        default=None,
        help="Optional x-position for a red dashed vertical marker on the Ne plot.",
    )
    p.add_argument(
        "--sim",
        type=int,
        default=0,
        help=(
            "Number of 1 Mb coalescent simulations to run under a Demes model built from "
            "the inferred Ne trajectory (default: 0 = no simulations)."
        ),
    )
    p.add_argument(
        "--mu",
        type=float,
        default=None,
        help=(
            "Mutation rate for optional simulations. Recombination rate is set equal to this "
            "value. Required when --sim > 0."
        ),
    )
    p.add_argument(
        "--sim-outfile",
        type=Path,
        default=None,
        help=(
            "Output TSV path for simulated 50 Kb window statistics. "
            "Default: <out-dir>/<prefix>sim-window-stats.tsv"
        ),
    )
    p.add_argument(
        "--sim-sfs-outfile",
        type=Path,
        default=None,
        help=(
            "Output TSV path for simulated site frequency spectra across simulations. "
            "Default: <out-dir>/<prefix>sim-sfs.tsv"
        ),
    )
    p.add_argument(
        "--sim-length",
        type=float,
        default=1_000_000.0,
        help="Sequence length (bp) for each optional simulation (default: 1000000).",
    )
    p.add_argument(
        "--sim-window-size",
        type=float,
        default=50_000.0,
        help="Window size (bp) for simulated diversity/Tajima's D TSV (default: 50000).",
    )
    p.add_argument(
        "--prefix",
        default="",
        help="Optional prefix prepended to output filenames (default: none).",
    )
    return p.parse_args()


def find_tree_files(ts_dir: Path, pattern: str) -> list[Path]:
    files = sorted(
        [p for p in ts_dir.glob(pattern) if p.is_file() and p.suffix in {".tsz", ".ts", ".trees"}]
    )
    if not files:
        raise RuntimeError(f"No tree files found in {ts_dir} matching pattern '{pattern}'.")
    return files


def load_time_windows(path: Path) -> np.ndarray:
    if not path.exists():
        raise FileNotFoundError(f"Time bins file not found: {path}")
    text = path.read_text().replace(",", " ")
    vals = [float(x) for x in text.split() if x.strip()]
    if len(vals) < 2:
        raise ValueError("Time bins file must contain at least two bin edges.")
    windows = np.array(vals, dtype=float)
    if not np.all(np.diff(windows) > 0):
        raise ValueError("Time bins must be strictly increasing.")
    if windows[0] > 0:
        windows = np.append([0.0], windows)
    if not np.isinf(windows[-1]):
        windows = np.append(windows, np.inf)
    return windows


def plottable_interval_mask(time_windows: np.ndarray) -> np.ndarray:
    left = time_windows[:-1]
    right = time_windows[1:]
    return np.logical_and(np.isfinite(right), left > 0)


def finite_interval_mask(time_windows: np.ndarray) -> np.ndarray:
    return np.isfinite(time_windows[1:])


def compute_pair_coal(ts, time_windows: np.ndarray, tail_cutoff: float):
    pdf = ts.pair_coalescence_counts(time_windows=time_windows, pair_normalise=True)
    rates = ts.pair_coalescence_rates(time_windows=time_windows)
    survival = np.append(1.0, 1.0 - np.cumsum(pdf))
    rates[survival[:-1] <= tail_cutoff] = np.nan
    return pdf, rates


def safe_nanmean(a: np.ndarray, axis=None):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        return np.nanmean(a, axis=axis)


def fill_ne_trajectory(ne: np.ndarray) -> np.ndarray:
    """
    Replace invalid Ne values with nearest-neighbor interpolation on interval index.
    """
    out = np.array(ne, dtype=float)
    valid = np.isfinite(out) & (out > 0)
    if not np.any(valid):
        raise RuntimeError("No valid Ne intervals available to build simulation model.")
    idx = np.arange(out.size)
    out[~valid] = np.interp(idx[~valid], idx[valid], out[valid])
    return out


def build_demes_graph_from_ne(time_windows: np.ndarray, ne_finite: np.ndarray):
    import demes

    finite_mask = finite_interval_mask(time_windows)
    lefts = time_windows[:-1][finite_mask]
    rights = time_windows[1:][finite_mask]
    if len(lefts) != len(ne_finite):
        raise RuntimeError("Ne trajectory length does not match finite time intervals.")
    # Demes epochs are listed oldest -> youngest, each with end_time and constant size.
    epochs = []
    for left, _right, ne in reversed(list(zip(lefts, rights, ne_finite))):
        epochs.append(
            {
                "end_time": float(left),
                "start_size": float(ne),
                "end_size": float(ne),
                "size_function": "constant",
            }
        )
    builder = demes.Builder(time_units="generations")
    builder.add_deme("pop", epochs=epochs)
    return builder.resolve()


def simulate_window_stats_from_ne(
    graph,
    n_samples: int,
    n_sims: int,
    mu: float,
    out_path: Path,
    sfs_out_path: Path,
    sequence_length: float = 1_000_000.0,
    window_size: float = 50_000.0,
):
    import msprime

    demography = msprime.Demography.from_demes(graph)
    windows = np.arange(0, sequence_length + window_size, window_size, dtype=float)
    if windows[-1] > sequence_length:
        windows[-1] = sequence_length

    lines = ["sim\twindow_index\tstart\tend\tdiversity\ttajima_d"]
    sfs_lines = ["sim\tmode\tbin\tvalue"]
    for sim_idx in range(n_sims):
        seed_base = 10_000 + sim_idx * 2
        ts = msprime.sim_ancestry(
            samples={"pop": int(n_samples)},
            demography=demography,
            sequence_length=sequence_length,
            recombination_rate=mu,
            random_seed=seed_base,
        )
        mts = msprime.sim_mutations(ts, rate=mu, random_seed=seed_base + 1, keep=False)
        pi = mts.diversity(mode="site", windows=windows)
        td = mts.Tajimas_D(mode="site", windows=windows)
        afs_pol = mts.allele_frequency_spectrum(mode="site", polarised=True)
        afs_fold = mts.allele_frequency_spectrum(mode="site", polarised=False)
        for w in range(len(windows) - 1):
            lines.append(
                "\t".join(
                    [
                        str(sim_idx),
                        str(w),
                        str(int(windows[w])),
                        str(int(windows[w + 1])),
                        f"{float(pi[w]):.10g}",
                        f"{float(td[w]):.10g}",
                    ]
                )
            )
        for i, val in enumerate(afs_pol):
            sfs_lines.append(f"{sim_idx}\tpolarised\t{i}\t{float(val):.10g}")
        for i, val in enumerate(afs_fold):
            sfs_lines.append(f"{sim_idx}\tfolded\t{i}\t{float(val):.10g}")
    out_path.write_text("\n".join(lines) + "\n")
    sfs_out_path.write_text("\n".join(sfs_lines) + "\n")
    return out_path, sfs_out_path


def main():
    args = parse_args()
    if args.time_adjust <= 0:
        raise ValueError("--time-adjust must be > 0")
    if args.sim < 0:
        raise ValueError("--sim must be >= 0")
    if args.sim > 0 and (args.mu is None or args.mu <= 0):
        raise ValueError("--mu must be > 0 when --sim > 0")
    if args.sim_length <= 0:
        raise ValueError("--sim-length must be > 0")
    if args.sim_window_size <= 0:
        raise ValueError("--sim-window-size must be > 0")
    ts_files = find_tree_files(args.ts_dir, args.pattern)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    time_windows = load_time_windows(args.time_bins_file)
    finite_mask = finite_interval_mask(time_windows)
    keep_mask = plottable_interval_mask(time_windows)
    breaks = time_windows[:-1][keep_mask]
    plot_breaks = breaks / args.time_adjust

    burnin = int(np.floor(len(ts_files) * args.burnin_frac))
    keep_idx = np.arange(len(ts_files))
    keep_post = keep_idx[burnin:] if burnin < len(ts_files) else keep_idx[-1:]

    pdf_vals_full = []
    rate_vals_full = []
    n_samples = None
    seq_lengths = []
    for ts_path in ts_files:
        ts = load_ts(ts_path)
        if n_samples is None:
            n_samples = ts.num_samples
        elif ts.num_samples != n_samples:
            raise RuntimeError(
                f"Sample count mismatch: {ts_path} has {ts.num_samples}, expected {n_samples}."
            )
        pdf, rates = compute_pair_coal(ts, time_windows, args.tail_cutoff)
        pdf_vals_full.append(pdf)
        rate_vals_full.append(rates)
        seq_lengths.append(float(ts.sequence_length))

    pdf_vals_full = np.stack(pdf_vals_full, axis=0)
    rate_vals_full = np.stack(rate_vals_full, axis=0)
    pdf_vals = pdf_vals_full[:, keep_mask]
    rate_vals = rate_vals_full[:, keep_mask]

    # Convert PMF to density on log scale: divide by log bin width so equal
    # area on the log x-axis equals equal probability.
    bin_widths = np.log(time_windows[1:][keep_mask] / time_windows[:-1][keep_mask])
    pdf_density = pdf_vals / bin_widths

    mean_pdf = safe_nanmean(pdf_density[keep_post], axis=0)
    mean_rates = safe_nanmean(rate_vals[keep_post], axis=0)

    reps_kwargs = {"color": "gray", "alpha": 0.15}
    mean_kwargs = {"color": "black", "linewidth": 1.5}

    fig, ax = plt.subplots(1, 1, figsize=(5, 4), constrained_layout=True)
    for val in pdf_density:
        ax.step(plot_breaks, val, **reps_kwargs)
    ax.step(plot_breaks, mean_pdf, **mean_kwargs)
    ax.set_xlabel("Adjusted generations in past")
    ax.set_ylabel("Coalescence density (proportion / log-generation)")
    ax.set_xscale("log")
    pdf_path = args.out_dir / f"{args.prefix}pair-coalescence-pdf.png"
    plt.savefig(pdf_path)
    plt.clf()

    fig, ax = plt.subplots(1, 1, figsize=(5, 4), constrained_layout=True)
    for val in rate_vals:
        ax.step(plot_breaks, val, **reps_kwargs)
    ax.step(plot_breaks, mean_rates, **mean_kwargs)
    ax.set_xlabel("Adjusted generations in past")
    ax.set_ylabel("Pair coalescence rate")
    ax.set_xscale("log")
    if args.log_rates:
        ax.set_yscale("log")
    rate_path = args.out_dir / f"{args.prefix}pair-coalescence-rates.png"
    plt.savefig(rate_path)
    plt.clf()

    ne_vals = np.full_like(rate_vals, np.nan, dtype=float)
    valid_rates = np.isfinite(rate_vals) & (rate_vals > 0)
    ne_vals[valid_rates] = 1.0 / (2.0 * rate_vals[valid_rates])
    mean_ne = safe_nanmean(ne_vals[keep_post], axis=0)
    fig, ax = plt.subplots(1, 1, figsize=(5, 4), constrained_layout=True)
    for val in ne_vals:
        ax.step(plot_breaks, val, **reps_kwargs)
    ax.step(plot_breaks, mean_ne, **mean_kwargs)
    if args.year is not None:
        ax.axvline(args.year, color="red", linestyle="--", linewidth=1.2)
    ax.set_xlabel("Adjusted generations in past")
    ax.set_ylabel("Estimated Ne = 1 / (2 * coal. rate)")
    ax.set_xscale("log")
    if args.log_rates:
        ax.set_yscale("log")
    ne_path = args.out_dir / f"{args.prefix}effective-pop-size.png"
    plt.savefig(ne_path)
    plt.clf()

    sim_out_path = None
    sim_sfs_out_path = None
    sim_graph = None
    if args.sim > 0:
        mean_rates_finite = safe_nanmean(rate_vals_full[keep_post][:, finite_mask], axis=0)
        ne_finite = np.full_like(mean_rates_finite, np.nan, dtype=float)
        valid = np.isfinite(mean_rates_finite) & (mean_rates_finite > 0)
        ne_finite[valid] = 1.0 / (2.0 * mean_rates_finite[valid])
        ne_finite = fill_ne_trajectory(ne_finite)
        sim_graph = build_demes_graph_from_ne(time_windows, ne_finite)
        sim_out_path = args.sim_outfile or (args.out_dir / f"{args.prefix}sim-window-stats.tsv")
        sim_sfs_out_path = args.sim_sfs_outfile or (args.out_dir / f"{args.prefix}sim-sfs.tsv")
        simulate_window_stats_from_ne(
            sim_graph,
            n_samples=n_samples,
            n_sims=args.sim,
            mu=float(args.mu),
            out_path=sim_out_path,
            sfs_out_path=sim_sfs_out_path,
            sequence_length=float(args.sim_length),
            window_size=float(args.sim_window_size),
        )

    summary_path = args.out_dir / f"{args.prefix}summary.txt"
    summary_path.write_text(
        "\n".join(
            [
                f"ts_dir={args.ts_dir}",
                f"pattern={args.pattern}",
                f"n_files={len(ts_files)}",
                f"burnin_frac={args.burnin_frac}",
                f"burnin_index={burnin}",
                f"time_bins_file={args.time_bins_file}",
                f"time_windows={time_windows.tolist()}",
                f"time_adjust={args.time_adjust}",
                f"year_marker={args.year}",
                f"tail_cutoff={args.tail_cutoff}",
                f"simulations={args.sim}",
                f"simulation_mu={args.mu}",
                f"simulation_length_bp={args.sim_length}",
                f"simulation_window_size_bp={args.sim_window_size}",
                f"n_samples={n_samples}",
                f"sequence_length_min={float(np.min(seq_lengths))}",
                f"sequence_length_max={float(np.max(seq_lengths))}",
                f"pair_coalescence_pdf_plot={pdf_path}",
                f"pair_coalescence_rates_plot={rate_path}",
                f"effective_pop_size_plot={ne_path}",
                f"simulation_window_stats_tsv={sim_out_path}",
                f"simulation_sfs_tsv={sim_sfs_out_path}",
            ]
        )
        + "\n"
    )

    print(f"Wrote: {pdf_path}")
    print(f"Wrote: {rate_path}")
    print(f"Wrote: {ne_path}")
    if sim_out_path is not None:
        print(f"Wrote: {sim_out_path}")
    if sim_sfs_out_path is not None:
        print(f"Wrote: {sim_sfs_out_path}")
    print(f"Wrote: {summary_path}")


if __name__ == "__main__":
    main()
