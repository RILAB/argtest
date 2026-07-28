#!/usr/bin/env python
"""Generate pair-coalescence and Ne plots from tree-sequence replicates.

Pair-coalescence quantiles, counts, and rates use tskit's native
partial-missing-data normalization, which conditions each result on the
pair-spans actually at risk of coalescence. Intervals where a sample is
isolated and has no MRCA with another sample therefore do not inflate the
survival denominator or bias Ne upward.
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path
import re
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
            "sequence replicates. Native tskit partial-missing-data "
            "normalization prevents isolated-sample regions (e.g. from "
            "trim_samples.py) from biasing the survival denominator."
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
        type=Path,
        default=None,
        help=(
            "File of explicit time-bin edges (one per line) defining the "
            "coalescence time grid. Mutually exclusive with --num-bins."
        ),
    )
    p.add_argument(
        "--num-quantiles",
        type=int,
        default=None,
        help=(
            "Choose this many equal-coalescence-event bins automatically. "
            "Edges are derived from conditional connected-pair coalescence "
            "quantiles on each post-burnin replicate (quantiles 0, 1/N, ..., "
            "1) and averaged across replicates. This conditions out pair-spans "
            "that cannot coalesce because samples are locally isolated. Zero "
            "and inf are then padded on so tskit's rate calculation accepts "
            "the grid. The two padded intervals are dropped from plotting, "
            "leaving N plotted bins each holding ~1/N of connected-pair "
            "coalescence mass. Mutually exclusive with --time-bins-file and "
            "--num-bins."
        ),
    )
    p.add_argument(
        "--num-bins",
        type=int,
        default=None,
        help=(
            "Choose this many bins spaced uniformly on a log scale across the "
            "coalescence-time range of the tree sequences. The range runs from "
            "the youngest to the oldest node time observed across the "
            "post-burnin replicates; N+1 log-spaced edges are laid over it, then "
            "0 and inf are padded on so tskit's rate calc accepts the grid. The "
            "two padded intervals are dropped from plotting, leaving N plotted "
            "bins of equal width in log-time. Mutually exclusive with "
            "--time-bins-file and --num-quantiles."
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
    files = [
        p
        for p in ts_dir.glob(pattern)
        if p.is_file() and p.suffix in {".tsz", ".ts", ".trees"}
    ]
    files.sort(
        key=lambda p: [
            int(part) if part.isdigit() else part.casefold()
            for part in re.split(r"(\d+)", str(p.relative_to(ts_dir)))
        ]
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


def compute_quantile_time_windows(
    ts_files: list[Path],
    post_burnin_indices: np.ndarray,
    num_quantiles: int,
) -> np.ndarray:
    """Build N+1 equal-connected-event bin edges averaged across replicates."""
    if num_quantiles < 2:
        raise ValueError("--num-quantiles must be >= 2.")
    if len(post_burnin_indices) == 0:
        raise RuntimeError("No post-burnin replicates available for --num-quantiles.")
    quantiles = np.linspace(0.0, 1.0, num_quantiles + 1)
    edge_arrays = []
    for idx in post_burnin_indices:
        ts = load_ts(ts_files[idx])
        edges = np.asarray(
            ts.pair_coalescence_quantiles(quantiles=quantiles),
            dtype=float,
        )
        edge_arrays.append(edges)
    mean_edges = np.mean(np.stack(edge_arrays, axis=0), axis=0)
    if not np.all(np.diff(mean_edges) > 0):
        # The conditional CDF is a step function with one weighted atom per
        # coalescence node. If one atom spans multiple requested probabilities,
        # adjacent quantiles have the same time and cannot define positive-width
        # bins. Nonfinite values here indicate an unexpected inversion failure,
        # since locally unconnected pair mass has already been conditioned out.
        undefined = [
            f"{q:.0%}" for q, e in zip(quantiles, mean_edges) if not np.isfinite(e)
        ]
        detail = (
            f" ({len(undefined)} of {mean_edges.size} quantile edges were "
            f"undefined, at quantiles {', '.join(undefined)})"
            if undefined
            else " (adjacent quantile edges are tied — a single node holds more "
            "than one full bin's worth of pair-coalescence mass)"
        )
        raise RuntimeError(
            "The conditional connected-pair coalescence CDF cannot define "
            f"{num_quantiles} positive-width equal-mass time bins"
            + detail
            + ". Re-run with fixed edges via --time-bins-file, with uniform "
            "log-spaced bins via --num-bins, or reduce --num-quantiles."
        )
    if mean_edges[0] <= 0:
        raise RuntimeError(
            "Averaged minimum coalescence-time edge is non-positive; cannot "
            "log-plot. Reduce --num-quantiles or supply fixed edges via "
            "--time-bins-file."
        )
    # Pad with 0 (tskit requires the grid to start at sample time) and inf
    # (tskit requires the rates grid to end at infinity). The two padded
    # intervals are empty in expectation and dropped by plottable_interval_mask.
    return np.concatenate(([0.0], mean_edges, [np.inf]))


def compute_logspaced_time_windows(
    ts_files: list[Path],
    post_burnin_indices: np.ndarray,
    num_bins: int,
) -> np.ndarray:
    """Build N uniform log-spaced bins across the coalescence-time range.

    The range is the youngest to oldest (non-sample) node time observed across
    the post-burnin replicates; N+1 log-spaced edges span it, and 0 / inf are
    padded on so tskit's rate calc accepts the grid. The two padded intervals
    are dropped by plottable_interval_mask, leaving N plotted bins."""
    if num_bins < 2:
        raise ValueError("--num-bins must be >= 2.")
    if len(post_burnin_indices) == 0:
        raise RuntimeError("No post-burnin replicates available for --num-bins.")
    lo = np.inf
    hi = 0.0
    for idx in post_burnin_indices:
        ts = load_ts(ts_files[idx])
        times = np.asarray(ts.tables.nodes.time, dtype=float)
        positive = times[times > 0]
        if positive.size:
            lo = min(lo, float(positive.min()))
            hi = max(hi, float(positive.max()))
    if not np.isfinite(lo) or hi <= lo:
        raise RuntimeError(
            "Could not determine a positive coalescence-time range for "
            "--num-bins (need at least two distinct positive node times across "
            "the post-burnin replicates). Supply fixed edges via "
            "--time-bins-file instead."
        )
    edges = np.logspace(np.log10(lo), np.log10(hi), num_bins + 1)
    return np.concatenate(([0.0], edges, [np.inf]))


def plottable_interval_mask(time_windows: np.ndarray) -> np.ndarray:
    left = time_windows[:-1]
    right = time_windows[1:]
    return np.logical_and(np.isfinite(right), left > 0)


def finite_interval_mask(time_windows: np.ndarray) -> np.ndarray:
    return np.isfinite(time_windows[1:])


def compute_pair_coal(ts, time_windows: np.ndarray, tail_cutoff: float):
    # The pinned tskit build normalizes over locally non-missing pair-spans.
    pdf = np.asarray(
        ts.pair_coalescence_counts(
            time_windows=time_windows,
            pair_normalise=True,
        ),
        dtype=float,
    )
    rates = np.asarray(
        ts.pair_coalescence_rates(time_windows=time_windows),
        dtype=float,
    )
    survival = np.append(1.0, 1.0 - np.cumsum(pdf))
    rates[survival[:-1] <= tail_cutoff] = np.nan
    return pdf, rates


def safe_nanmean(a: np.ndarray, axis=None):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        return np.nanmean(a, axis=axis)


def plot_postburn_replicates(ax, x, values, postburn_indices, **kwargs):
    """Draw only replicate trajectories that contribute to the mean."""
    for val in values[postburn_indices]:
        ax.step(x, val, **kwargs)


def write_coalescence_estimates(
    out_path: Path,
    ts_files: list[Path],
    postburn_indices: np.ndarray,
    time_windows: np.ndarray,
    keep_mask: np.ndarray,
    time_adjust: float,
    pdf_mass: np.ndarray,
    pdf_density: np.ndarray,
    rates: np.ndarray,
    ne: np.ndarray,
    mean_pdf_mass: np.ndarray,
    mean_pdf_density: np.ndarray,
    mean_rates: np.ndarray,
    mean_ne: np.ndarray,
) -> Path:
    """Save plotted post-burnin trajectories and their posterior means."""
    left = time_windows[:-1][keep_mask]
    right = time_windows[1:][keep_mask]
    header = [
        "series",
        "replicate_index",
        "tree_file",
        "bin_index",
        "time_left",
        "time_right",
        "adjusted_time_left",
        "adjusted_time_right",
        "pair_coalescence_mass",
        "pair_coalescence_log_density",
        "pair_coalescence_rate",
        "effective_population_size",
    ]
    lines = ["\t".join(header)]

    def add_series(
        series,
        replicate_index,
        tree_file,
        mass_values,
        density_values,
        rate_values,
        ne_values,
    ):
        for j in range(left.size):
            row = [
                series,
                str(replicate_index),
                str(tree_file),
                str(j),
                f"{left[j]:.10g}",
                f"{right[j]:.10g}",
                f"{left[j] / time_adjust:.10g}",
                f"{right[j] / time_adjust:.10g}",
                f"{mass_values[j]:.10g}",
                f"{density_values[j]:.10g}",
                f"{rate_values[j]:.10g}",
                f"{ne_values[j]:.10g}",
            ]
            lines.append("\t".join(row))

    for i in postburn_indices:
        add_series(
            "replicate",
            int(i),
            ts_files[i],
            pdf_mass[i],
            pdf_density[i],
            rates[i],
            ne[i],
        )
    add_series(
        "posterior_mean",
        "",
        "",
        mean_pdf_mass,
        mean_pdf_density,
        mean_rates,
        mean_ne,
    )
    out_path.write_text("\n".join(lines) + "\n")
    return out_path


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
    grid_opts = [
        args.time_bins_file is not None,
        args.num_quantiles is not None,
        args.num_bins is not None,
    ]
    if sum(grid_opts) != 1:
        raise ValueError(
            "Provide exactly one of --time-bins-file, --num-quantiles, or "
            "--num-bins."
        )
    if args.num_quantiles is not None and args.num_quantiles < 2:
        raise ValueError("--num-quantiles must be >= 2.")
    if args.num_bins is not None and args.num_bins < 2:
        raise ValueError("--num-bins must be >= 2.")
    ts_files = find_tree_files(args.ts_dir, args.pattern)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    burnin = int(np.floor(len(ts_files) * args.burnin_frac))
    keep_idx = np.arange(len(ts_files))
    keep_post = keep_idx[burnin:] if burnin < len(ts_files) else keep_idx[-1:]

    if args.time_bins_file is not None:
        time_windows = load_time_windows(args.time_bins_file)
    elif args.num_quantiles is not None:
        time_windows = compute_quantile_time_windows(
            ts_files, keep_post, args.num_quantiles
        )
    else:
        time_windows = compute_logspaced_time_windows(
            ts_files, keep_post, args.num_bins
        )
    finite_mask = finite_interval_mask(time_windows)
    keep_mask = plottable_interval_mask(time_windows)
    breaks = time_windows[:-1][keep_mask]
    plot_breaks = breaks / args.time_adjust

    pdf_vals = []
    rate_vals = []
    rate_vals_finite_post = []
    n_samples = None
    seq_lengths = []
    for i, ts_path in enumerate(ts_files):
        ts = load_ts(ts_path)
        if n_samples is None:
            n_samples = ts.num_samples
        elif ts.num_samples != n_samples:
            raise RuntimeError(
                f"Sample count mismatch: {ts_path} has {ts.num_samples}, expected {n_samples}."
            )
        pdf, rates = compute_pair_coal(ts, time_windows, args.tail_cutoff)
        pdf_vals.append(pdf[keep_mask])
        rate_vals.append(rates[keep_mask])
        if i in keep_post:
            rate_vals_finite_post.append(rates[finite_mask])
        seq_lengths.append(float(ts.sequence_length))

    pdf_vals = np.stack(pdf_vals, axis=0)
    rate_vals = np.stack(rate_vals, axis=0)

    # Convert PMF to density on log scale: divide by log bin width so equal
    # area on the log x-axis equals equal probability.
    bin_widths = np.log(time_windows[1:][keep_mask] / time_windows[:-1][keep_mask])
    pdf_density = pdf_vals / bin_widths

    mean_pdf_mass = safe_nanmean(pdf_vals[keep_post], axis=0)
    mean_pdf = safe_nanmean(pdf_density[keep_post], axis=0)
    mean_rates = safe_nanmean(rate_vals[keep_post], axis=0)

    reps_kwargs = {"color": "gray", "alpha": 0.15}
    mean_kwargs = {"color": "black", "linewidth": 1.5}

    fig, ax = plt.subplots(1, 1, figsize=(5, 4), constrained_layout=True)
    plot_postburn_replicates(
        ax, plot_breaks, pdf_density, keep_post, **reps_kwargs
    )
    ax.step(plot_breaks, mean_pdf, **mean_kwargs)
    ax.set_xlabel("Adjusted generations in past")
    ax.set_ylabel("Coalescence density (proportion / log-generation)")
    ax.set_xscale("log")
    pdf_path = args.out_dir / f"{args.prefix}pair-coalescence-pdf.png"
    plt.savefig(pdf_path)
    plt.clf()

    fig, ax = plt.subplots(1, 1, figsize=(5, 4), constrained_layout=True)
    plot_postburn_replicates(ax, plot_breaks, rate_vals, keep_post, **reps_kwargs)
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
    plot_postburn_replicates(ax, plot_breaks, ne_vals, keep_post, **reps_kwargs)
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

    estimates_path = write_coalescence_estimates(
        args.out_dir / f"{args.prefix}coalescence-ne-estimates.tsv",
        ts_files,
        keep_post,
        time_windows,
        keep_mask,
        args.time_adjust,
        pdf_vals,
        pdf_density,
        rate_vals,
        ne_vals,
        mean_pdf_mass,
        mean_pdf,
        mean_rates,
        mean_ne,
    )

    sim_out_path = None
    sim_sfs_out_path = None
    sim_graph = None
    if args.sim > 0:
        mean_rates_finite = safe_nanmean(np.stack(rate_vals_finite_post, axis=0), axis=0)
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
                f"num_quantiles={args.num_quantiles}",
                f"num_bins={args.num_bins}",
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
                f"coalescence_estimates_tsv={estimates_path}",
                f"simulation_window_stats_tsv={sim_out_path}",
                f"simulation_sfs_tsv={sim_sfs_out_path}",
            ]
        )
        + "\n"
    )

    print(f"Wrote: {pdf_path}")
    print(f"Wrote: {rate_path}")
    print(f"Wrote: {ne_path}")
    print(f"Wrote: {estimates_path}")
    if sim_out_path is not None:
        print(f"Wrote: {sim_out_path}")
    if sim_sfs_out_path is not None:
        print(f"Wrote: {sim_sfs_out_path}")
    print(f"Wrote: {summary_path}")


if __name__ == "__main__":
    main()
