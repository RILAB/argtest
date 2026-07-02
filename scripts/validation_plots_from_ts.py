#!/usr/bin/env python
from __future__ import annotations

import argparse
import os
import pickle
import re
from pathlib import Path
import warnings
import csv

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp")

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

try:
    import msprime
except ImportError:
    msprime = None

from argtest_common import (
    accessible_intervals_from_mu,
    infer_mu_path,
    load_ts,
    mutational_load,
    ratemap_from_metadata,
    tree_covered_accessible_bp,
)

matplotlib.rcParams["figure.dpi"] = 300


def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Generate SINGER-style validation/diagnostic plots from a set of tree sequences."
        )
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
        default=Path("results/validation_plots"),
        help="Output directory for plots.",
    )
    p.add_argument(
        "--burnin-frac",
        type=float,
        default=0.5,
        help="Fraction of initial trees to ignore when computing posterior means.",
    )
    p.add_argument(
        "--prefix",
        default="",
        help="Optional prefix for output plot filenames.",
    )
    p.add_argument(
        "--folded",
        action="store_true",
        help="Plot folded SFS (minor-allele frequency) instead of polarised derived-frequency SFS.",
    )
    p.add_argument(
        "--sim",
        type=Path,
        default=None,
        help=(
            "Optional TSV from scripts/coalescence_ne_plots_from_ts.py --sim output "
            "(sim-window-stats.tsv) for observed-vs-simulated pi/Tajima's D density plots."
        ),
    )
    p.add_argument(
        "--sim-sfs",
        type=Path,
        default=None,
        help=(
            "Optional TSV from scripts/coalescence_ne_plots_from_ts.py simulated SFS output "
            "(sim-sfs.tsv) for observed-vs-simulated SFS plot."
        ),
    )
    p.add_argument(
        "--window-size",
        type=float,
        default=1.0e5,
        help=(
            "Window size (bp) for the diversity, Tajima's D, and "
            "segregating-sites validation plots (default: 100000)."
        ),
    )
    p.add_argument(
        "--mutation-rate",
        type=float,
        default=None,
        help=(
            "Fallback uniform mutation rate for --sim-branch when no *.mut_rate.p "
            "ratemap is discoverable. Optional when --sim-branch can use a ratemap; "
            "unused when --sim-branch is not set."
        ),
    )
    p.add_argument(
        "--compare",
        type=Path,
        default=None,
        help=(
            "Optional second tree-sequence directory to overlay on all plots for comparison "
            "(e.g. pre- vs post-pipeline). Uses the same --pattern, --window-size, "
            "--burnin-frac, and --mutation-rate as the primary directory."
        ),
    )
    p.add_argument(
        "--sim-branch",
        action="store_true",
        help=(
            "Use msprime.sim_mutations to simulate site mutations on each ARG replicate, "
            "then compute site-mode statistics on the simulated TS — a posterior predictive "
            "check matching the original nspope/singer-snakemake approach. Requires a "
            "*.mut_rate.p file co-located with the tree sequences; falls back to uniform "
            "--mutation-rate if the file cannot be found."
        ),
    )
    p.add_argument(
        "--mcmc-thin",
        type=int,
        default=1,
        help=(
            "Thinning interval between replicate IDs for trace x-axes. "
            "MCMC iteration is computed as replicate_id * mcmc_thin (default: 1)."
        ),
    )
    return p.parse_args()


def natural_sort_key(text: str):
    parts = re.split(r"(\d+)", text)
    key = []
    for part in parts:
        if part.isdigit():
            key.append((0, int(part)))
        else:
            key.append((1, part.lower()))
    return key


def extract_replicate_id(path: Path):
    m = re.search(r"(\d+)$", path.stem)
    return int(m.group(1)) if m else None


def find_tree_files(ts_dir: Path, pattern: str) -> list[Path]:
    files = sorted(
        [p for p in ts_dir.glob(pattern) if p.is_file() and p.suffix in {".tsz", ".ts", ".trees"}],
        key=lambda p: natural_sort_key(p.stem),
    )
    if not files:
        raise RuntimeError(f"No tree files found in {ts_dir} matching pattern '{pattern}'.")
    return files


def _try_mu_ratemap(ts_files: list[Path]):
    """Return a RateMap from ts metadata (preferred) or nearest *.mut_rate.p file, or None."""
    if not ts_files:
        return None
    try:
        mu = ratemap_from_metadata(load_ts(ts_files[0]).metadata or {})
        if mu is not None:
            return mu
    except Exception:
        pass
    try:
        mu_path = infer_mu_path(ts_files[0])
        with open(mu_path, "rb") as fh:
            return pickle.load(fh)
    except Exception:
        return None


def _try_mu_intervals(ts_files: list[Path]):
    """Return accessible intervals from the nearest *.mut_rate.p file, or None if not found."""
    mu = _try_mu_ratemap(ts_files)
    return accessible_intervals_from_mu(mu) if mu is not None else None


def genome_windows(sequence_length: float, window_size: float) -> np.ndarray:
    windows = np.arange(0, float(sequence_length) + float(window_size), float(window_size), dtype=float)
    if windows[-1] > float(sequence_length):
        windows[-1] = float(sequence_length)
    return windows


def accessible_per_window(stat_windows: np.ndarray, kept_intervals: list) -> np.ndarray:
    """Return accessible bp per window given a list of [[left, right], ...] kept intervals."""
    wl = stat_windows[:-1]
    wr = stat_windows[1:]
    acc = np.zeros(len(wl), dtype=float)
    for left, right in kept_intervals:
        acc += np.maximum(0.0, np.minimum(right, wr) - np.maximum(left, wl))
    return acc


def safe_nanmean(a: np.ndarray, axis=None):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        return np.nanmean(a, axis=axis)


def safe_nanquantile(a: np.ndarray, q, axis=None):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        return np.nanquantile(a, q, axis=axis)


def safe_log_yscale(ax, *arrays, context: str = "plot") -> None:
    """Apply a log y-scale only if some plotted value is positive and finite.

    matplotlib raises "Data cannot be log-scaled because all values are <= 0"
    at savefig time when every value is zero/NaN (e.g. a site-frequency
    spectrum from tree sequences that carry no mutations). In that case leave
    the axis linear and warn instead of crashing the whole pipeline.
    """
    has_positive = any(
        np.any(np.isfinite(a) & (np.asarray(a, dtype=float) > 0))
        for a in arrays
        if a is not None and np.size(a) > 0
    )
    if has_positive:
        ax.set_yscale("log")
    else:
        warnings.warn(
            f"{context}: no positive values to plot (empty site spectrum?); "
            "leaving y-axis linear instead of log-scaled."
        )


def load_sim_window_stats(path: Path) -> tuple[np.ndarray, np.ndarray]:
    if not path.exists():
        raise FileNotFoundError(f"Simulation TSV not found: {path}")
    sim_pi = []
    sim_td = []
    with open(path, "r", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        required = {"diversity", "tajima_d"}
        if reader.fieldnames is None or not required.issubset(set(reader.fieldnames)):
            raise ValueError(
                f"Simulation TSV {path} must contain columns: {sorted(required)}"
            )
        for row in reader:
            try:
                sim_pi.append(float(row["diversity"]))
            except Exception:
                sim_pi.append(np.nan)
            try:
                sim_td.append(float(row["tajima_d"]))
            except Exception:
                sim_td.append(np.nan)
    return np.array(sim_pi, dtype=float), np.array(sim_td, dtype=float)


def load_sim_sfs(path: Path, *, folded: bool) -> np.ndarray:
    if not path.exists():
        raise FileNotFoundError(f"Simulation SFS TSV not found: {path}")
    mode_keep = "folded" if folded else "polarised"
    rows = []
    with open(path, "r", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        required = {"sim", "mode", "bin", "value"}
        if reader.fieldnames is None or not required.issubset(set(reader.fieldnames)):
            raise ValueError(f"Simulation SFS TSV {path} missing required columns {sorted(required)}")
        for row in reader:
            if row["mode"] != mode_keep:
                continue
            rows.append((int(row["sim"]), int(row["bin"]), float(row["value"])))
    if not rows:
        raise ValueError(f"No rows for mode='{mode_keep}' found in {path}")
    n_sim = max(r[0] for r in rows) + 1
    n_bin = max(r[1] for r in rows) + 1
    out = np.full((n_bin, n_sim), np.nan, dtype=float)
    for sim_idx, bin_idx, val in rows:
        out[bin_idx, sim_idx] = val
    return out


def collect_stats(ts_files: list[Path], window_size: float, burnin_frac: float,
                  mutation_rate: float,
                  mcmc_thin: int = 1,
                  sim_branch: bool = False, mu_ratemap=None) -> dict:
    """Load all tree sequences and return per-replicate site statistics.

    When sim_branch=True, additionally simulate site mutations on each replicate with
    msprime.sim_mutations and return paired 'sim' statistics for posterior predictive checks.
    """
    if sim_branch and msprime is None:
        raise RuntimeError("msprime is required for --sim-branch mode")

    _sim_rate = mu_ratemap if mu_ratemap is not None else mutation_rate

    load_vals, seq_lengths = [], []
    sim_load_vals = [] if sim_branch else None
    site_div_vals = []
    site_td_vals = []
    site_s_vals = []
    site_afs_unfolded_vals = []
    site_afs_folded_vals = []
    sim_div_vals = [] if sim_branch else None
    sim_td_vals = [] if sim_branch else None
    sim_s_vals = [] if sim_branch else None
    sim_afs_unfolded_vals = [] if sim_branch else None
    sim_afs_folded_vals = [] if sim_branch else None
    n_samples = None

    if mu_ratemap is not None:
        mu_intervals = accessible_intervals_from_mu(mu_ratemap)
    else:
        mu_intervals = _try_mu_intervals(ts_files)

    parsed_rep_ids = [extract_replicate_id(p) for p in ts_files]
    replicate_ids = np.array(
        [rid if rid is not None else i for i, rid in enumerate(parsed_rep_ids)],
        dtype=int,
    )

    for i, (ts_path, rep_id) in enumerate(zip(ts_files, replicate_ids)):
        ts = load_ts(ts_path)
        # tskit requires windows[-1] == sequence_length exactly, so each replicate
        # gets its own window grid.
        rep_windows = genome_windows(ts.sequence_length, window_size)
        if n_samples is None:
            n_samples = ts.num_samples
        elif ts.num_samples != n_samples:
            raise RuntimeError(
                f"Sample count mismatch: {ts_path} has {ts.num_samples}, expected {n_samples}."
            )

        # Determine accessible intervals for this replicate:
        #   1. kept_intervals metadata (written by trim_regions_single) — exact post-pipeline mask
        #   2. mu_intervals from *.mut_rate.p — pre-pipeline accessibility (rate > 0)
        #   3. fallback: treat entire sequence as accessible
        kept_intervals = ts.metadata.get("kept_intervals") if ts.metadata else None
        if kept_intervals is not None:
            acc_intervals = np.asarray(kept_intervals, dtype=float)
        elif mu_intervals is not None:
            acc_intervals = mu_intervals
        else:
            acc_intervals = None
        total_accessible = tree_covered_accessible_bp(ts, acc_intervals)

        if acc_intervals is not None:
            rep_acc = accessible_per_window(rep_windows, acc_intervals)
            window_spans = rep_windows[1:] - rep_windows[:-1]
            with np.errstate(invalid="ignore", divide="ignore"):
                div_scale = np.where(rep_acc > 0, window_spans / rep_acc, np.nan)
        else:
            rep_acc = None
            div_scale = None

        load = mutational_load(ts) / total_accessible
        load_vals.append(load)
        seq_lengths.append(float(ts.sequence_length))

        site_div = ts.diversity(mode="site", windows=rep_windows)
        site_td = ts.Tajimas_D(mode="site", windows=rep_windows)
        site_s_raw = ts.segregating_sites(
            mode="site", windows=rep_windows, span_normalise=False
        )
        if rep_acc is not None:
            with np.errstate(invalid="ignore", divide="ignore"):
                site_s = np.where(rep_acc > 0, site_s_raw / rep_acc, np.nan)
        else:
            window_spans = rep_windows[1:] - rep_windows[:-1]
            with np.errstate(invalid="ignore", divide="ignore"):
                site_s = site_s_raw / window_spans
        if div_scale is not None:
            site_div = site_div * div_scale
            site_td = np.where(np.isnan(div_scale), np.nan, site_td)
        site_afs_unfolded = (
            ts.allele_frequency_spectrum(mode="site", polarised=True, span_normalise=False)
            / total_accessible
        )
        site_afs_folded = (
            ts.allele_frequency_spectrum(mode="site", polarised=False, span_normalise=False)
            / total_accessible
        )
        site_div_vals.append(site_div)
        site_td_vals.append(site_td)
        site_s_vals.append(site_s)
        site_afs_unfolded_vals.append(site_afs_unfolded)
        site_afs_folded_vals.append(site_afs_folded)

        if sim_branch:
            # Posterior predictive check: simulate mutations on the ARG topology at
            # the observed rate, then compute site-mode stats on the result. Seed is
            # deterministic from replicate id so runs are reproducible.
            sim_ts = msprime.sim_mutations(
                ts, rate=_sim_rate, keep=False, random_seed=1 + int(rep_id) * 1000
            )
            # Per-sample expected load from this simulation draw (matches the
            # observed per-sample load normalization above).
            sim_load_vals.append(mutational_load(sim_ts) / total_accessible)
            sim_div_raw = sim_ts.diversity(
                mode="site", windows=rep_windows, span_normalise=False
            )
            sim_td = sim_ts.Tajimas_D(mode="site", windows=rep_windows)
            sim_s_raw = sim_ts.segregating_sites(
                mode="site", windows=rep_windows, span_normalise=False
            )
            if rep_acc is not None:
                with np.errstate(invalid="ignore", divide="ignore"):
                    sim_div = np.where(rep_acc > 0, sim_div_raw / rep_acc, np.nan)
                    sim_s = np.where(rep_acc > 0, sim_s_raw / rep_acc, np.nan)
                sim_td = np.where(rep_acc <= 0, np.nan, sim_td)
            else:
                window_spans = rep_windows[1:] - rep_windows[:-1]
                with np.errstate(invalid="ignore", divide="ignore"):
                    sim_div = sim_div_raw / window_spans
                    sim_s = sim_s_raw / window_spans
            sim_afs_unfolded = (
                sim_ts.allele_frequency_spectrum(mode="site", polarised=True, span_normalise=False)
                / total_accessible
            )
            sim_afs_folded = (
                sim_ts.allele_frequency_spectrum(mode="site", polarised=False, span_normalise=False)
                / total_accessible
            )
            sim_div_vals.append(sim_div)
            sim_td_vals.append(sim_td)
            sim_s_vals.append(sim_s)
            sim_afs_unfolded_vals.append(sim_afs_unfolded)
            sim_afs_folded_vals.append(sim_afs_folded)

    min_n_windows = min(len(a) for a in site_div_vals)
    site_div_vals = [a[:min_n_windows] for a in site_div_vals]
    site_td_vals = [a[:min_n_windows] for a in site_td_vals]
    site_s_vals = [a[:min_n_windows] for a in site_s_vals]
    if sim_branch:
        sim_div_vals = [a[:min_n_windows] for a in sim_div_vals]
        sim_td_vals = [a[:min_n_windows] for a in sim_td_vals]
        sim_s_vals = [a[:min_n_windows] for a in sim_s_vals]
    stat_windows = genome_windows(min(seq_lengths), window_size)

    n_files = len(ts_files)
    burnin = int(np.floor(n_files * burnin_frac))
    keep_idx = np.arange(n_files)
    mcmc_iterates = replicate_ids * int(mcmc_thin)
    burnin_iterate = int(mcmc_iterates[burnin]) if burnin < n_files else int(mcmc_iterates[-1])
    keep_post = keep_idx[burnin:] if burnin < n_files else keep_idx[-1:]

    load_vals = np.stack(load_vals, axis=-1)
    site_div_vals = np.stack(site_div_vals, axis=-1)
    site_td_vals = np.stack(site_td_vals, axis=-1)
    site_s_vals = np.stack(site_s_vals, axis=-1)
    site_afs_unfolded_vals = np.stack(site_afs_unfolded_vals, axis=-1)
    site_afs_folded_vals = np.stack(site_afs_folded_vals, axis=-1)
    coord = stat_windows[:-1] / 2.0 + stat_windows[1:] / 2.0

    out = {
        "n_files": n_files,
        "n_samples": n_samples,
        "seq_lengths": np.array(seq_lengths),
        "coord": coord,
        "keep_idx": keep_idx,
        "replicate_ids": replicate_ids,
        "mcmc_iterates": mcmc_iterates,
        "burnin": burnin,
        "burnin_iterate": burnin_iterate,
        "keep_post": keep_post,
        "sim_branch": sim_branch,
        "load_vals": load_vals,
        "site_div_vals": site_div_vals,
        "site_td_vals": site_td_vals,
        "site_s_vals": site_s_vals,
        "mean_load": safe_nanmean(load_vals[:, keep_post], axis=-1),
        "q_load": safe_nanquantile(load_vals[:, keep_post], [0.025, 0.975], axis=-1),
        "mean_site_div": safe_nanmean(site_div_vals[:, keep_post], axis=-1),
        "mean_site_td": safe_nanmean(site_td_vals[:, keep_post], axis=-1),
        "mean_site_s": safe_nanmean(site_s_vals[:, keep_post], axis=-1),
        "trace_site_div": safe_nanmean(site_div_vals, axis=0),
        "trace_site_td": safe_nanmean(site_td_vals, axis=0),
        "trace_site_s": safe_nanmean(site_s_vals, axis=0),
        "mean_site_afs_unfolded": safe_nanmean(site_afs_unfolded_vals[:, keep_post], axis=-1),
        "mean_site_afs_folded": safe_nanmean(site_afs_folded_vals[:, keep_post], axis=-1),
    }
    if sim_branch:
        sim_load_vals = np.stack(sim_load_vals, axis=-1)
        sim_div_vals = np.stack(sim_div_vals, axis=-1)
        sim_td_vals = np.stack(sim_td_vals, axis=-1)
        sim_s_vals = np.stack(sim_s_vals, axis=-1)
        sim_afs_unfolded_vals = np.stack(sim_afs_unfolded_vals, axis=-1)
        sim_afs_folded_vals = np.stack(sim_afs_folded_vals, axis=-1)
        out.update({
            "sim_load_vals": sim_load_vals,
            "mean_sim_load": safe_nanmean(sim_load_vals[:, keep_post], axis=-1),
            "q_sim_load": safe_nanquantile(sim_load_vals[:, keep_post], [0.025, 0.975], axis=-1),
            "sim_div_vals": sim_div_vals,
            "sim_td_vals": sim_td_vals,
            "sim_s_vals": sim_s_vals,
            "mean_sim_div": safe_nanmean(sim_div_vals[:, keep_post], axis=-1),
            "q_sim_div": safe_nanquantile(sim_div_vals[:, keep_post], [0.025, 0.975], axis=-1),
            "mean_sim_td": safe_nanmean(sim_td_vals[:, keep_post], axis=-1),
            "q_sim_td": safe_nanquantile(sim_td_vals[:, keep_post], [0.025, 0.975], axis=-1),
            "mean_sim_s": safe_nanmean(sim_s_vals[:, keep_post], axis=-1),
            "q_sim_s": safe_nanquantile(sim_s_vals[:, keep_post], [0.025, 0.975], axis=-1),
            "trace_sim_div": safe_nanmean(sim_div_vals, axis=0),
            "trace_sim_td": safe_nanmean(sim_td_vals, axis=0),
            "trace_sim_s": safe_nanmean(sim_s_vals, axis=0),
            "mean_sim_afs_unfolded": safe_nanmean(sim_afs_unfolded_vals[:, keep_post], axis=-1),
            "q_sim_afs_unfolded": safe_nanquantile(sim_afs_unfolded_vals[:, keep_post], [0.025, 0.975], axis=-1),
            "mean_sim_afs_folded": safe_nanmean(sim_afs_folded_vals[:, keep_post], axis=-1),
            "q_sim_afs_folded": safe_nanquantile(sim_afs_folded_vals[:, keep_post], [0.025, 0.975], axis=-1),
        })
    return out


def main():
    args = parse_args()
    if args.mcmc_thin <= 0:
        raise ValueError("--mcmc-thin must be > 0.")
    ts_files = find_tree_files(args.ts_dir, args.pattern)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    mu_ratemap = None
    if args.sim_branch:
        if msprime is None:
            raise RuntimeError("msprime is required for --sim-branch mode")
        mu_ratemap = _try_mu_ratemap(ts_files)
        if mu_ratemap is None:
            if args.mutation_rate is None:
                raise RuntimeError(
                    "--sim-branch requires either an embedded/sibling mutation-rate "
                    "ratemap or --mutation-rate as a scalar fallback."
                )
            warnings.warn(
                "--sim-branch requested but no *.mut_rate.p file found; "
                "simulating with uniform rate (--mutation-rate)."
            )

    pri = collect_stats(
        ts_files, args.window_size, args.burnin_frac, args.mutation_rate,
        mcmc_thin=args.mcmc_thin,
        sim_branch=args.sim_branch, mu_ratemap=mu_ratemap,
    )
    pri_label = args.ts_dir.name

    cmp = None
    cmp_label = None
    if args.compare is not None:
        cmp_files = find_tree_files(args.compare, args.pattern)
        cmp_mu_ratemap = _try_mu_ratemap(cmp_files) if args.sim_branch else None
        if args.sim_branch and cmp_mu_ratemap is None and args.mutation_rate is None:
            raise RuntimeError(
                "--sim-branch with --compare requires the comparison directory to "
                "have an embedded/sibling mutation-rate ratemap, or --mutation-rate "
                "as a scalar fallback."
            )
        cmp = collect_stats(
            cmp_files, args.window_size, args.burnin_frac, args.mutation_rate,
            mcmc_thin=args.mcmc_thin,
            sim_branch=args.sim_branch, mu_ratemap=cmp_mu_ratemap,
        )
        cmp_label = args.compare.name

    # Color / style scheme:
    #   black = observed site stats, firebrick = sim overlay (only with --sim-branch)
    #   solid = primary dataset, dashed = --compare dataset
    PRI_SITE = "black"
    PRI_SIM = "firebrick"
    CMP_SITE = "dimgray"
    CMP_SIM = "steelblue"
    CMP_LS = "--"

    # ------------------------------------------------------------------ #
    # Mutational load by sample (posterior summary)
    # ------------------------------------------------------------------ #
    n_width = max(pri["mean_load"].size, cmp["mean_load"].size if cmp else 0)
    fig, ax = plt.subplots(1, 1, figsize=(max(5, 0.1 * n_width), 4), constrained_layout=True)
    samples = np.arange(pri["mean_load"].size)
    # Expected load: with --sim-branch, plot the per-sample simulated expectation
    # (mutations dropped on each ARG) so each sample is compared to its own
    # genealogy-based expectation; otherwise fall back to a flat genome-wide mean.
    if "mean_sim_load" in pri:
        ax.plot(samples, pri["mean_sim_load"], "_", color=PRI_SIM, markersize=8,
                markeredgewidth=1.5, label=f"{pri_label} expected (sim)")
        ax.vlines(samples, pri["q_sim_load"][0], pri["q_sim_load"][1], color=PRI_SIM,
                  alpha=0.5, linewidth=3)
    else:
        ax.axhline(y=float(np.nanmean(pri["mean_load"])), color=PRI_SIM, linestyle="dashed")
    ax.plot(samples, pri["mean_load"], "o", color=PRI_SITE, markersize=2, label=pri_label)
    ax.vlines(samples, pri["q_load"][0], pri["q_load"][1], color=PRI_SITE)
    if cmp is not None:
        c_samples = np.arange(cmp["mean_load"].size)
        if "mean_sim_load" in cmp:
            ax.plot(c_samples, cmp["mean_sim_load"], "_", color=CMP_SIM, markersize=8,
                    markeredgewidth=1.5, alpha=0.7, label=f"{cmp_label} expected (sim)")
            ax.vlines(c_samples, cmp["q_sim_load"][0], cmp["q_sim_load"][1], color=CMP_SIM,
                      alpha=0.4, linewidth=3)
        else:
            ax.axhline(y=float(np.nanmean(cmp["mean_load"])), color=CMP_SIM, linestyle="dashed",
                       alpha=0.7)
        ax.plot(c_samples, cmp["mean_load"], "o", color=CMP_SIM, markersize=2, label=cmp_label)
        ax.vlines(c_samples, cmp["q_load"][0], cmp["q_load"][1], color=CMP_SIM, alpha=0.6)
    ax.set_xlabel("Sample ID")
    ax.set_ylabel("Derived mutations / base")
    ax.legend()
    load_path = args.out_dir / f"{args.prefix}mutational-load.png"
    plt.savefig(load_path)
    plt.clf()

    # ------------------------------------------------------------------ #
    # Mutational load trace across replicates
    # ------------------------------------------------------------------ #
    fig, ax = plt.subplots(1, 1, figsize=(5, 4), constrained_layout=True)
    for i in range(pri["load_vals"].shape[0]):
        ax.plot(pri["mcmc_iterates"], pri["load_vals"][i], "-", color=PRI_SITE,
                linewidth=0.5, alpha=0.2)
    ax.axvline(pri["burnin_iterate"], color=PRI_SIM, linestyle="--", linewidth=1,
               label=f"{pri_label} burnin")
    if cmp is not None:
        for i in range(cmp["load_vals"].shape[0]):
            ax.plot(cmp["mcmc_iterates"], cmp["load_vals"][i], "-", color=CMP_SIM,
                    linewidth=0.5, alpha=0.2)
        ax.axvline(cmp["burnin_iterate"], color=CMP_SIM, linestyle="--", linewidth=1,
                   label=f"{cmp_label} burnin")
    ax.plot([], [], "-", color=PRI_SITE, linewidth=1, label=pri_label)
    if cmp is not None:
        ax.plot([], [], "-", color=CMP_SIM, linewidth=1, label=cmp_label)
    ax.set_xlabel("MCMC iteration")
    ax.set_ylabel("Derived mutations / base in each sample")
    ax.legend()
    trace_path = args.out_dir / f"{args.prefix}mutational-load-trace.png"
    plt.savefig(trace_path)
    plt.clf()

    # ------------------------------------------------------------------ #
    # Diversity: site-vs-sim scatter (only meaningful with --sim-branch)
    # ------------------------------------------------------------------ #
    div_scatter_path = None
    if args.sim_branch:
        fig, ax = plt.subplots(1, 1, figsize=(5, 4), constrained_layout=True)
        ax.scatter(pri["mean_site_div"], pri["mean_sim_div"], c=PRI_SIM, s=8, label=pri_label)
        if cmp is not None:
            ax.scatter(cmp["mean_site_div"], cmp["mean_sim_div"], c=CMP_SIM, s=8,
                       marker="^", label=cmp_label)
        x = float(np.nanmean(pri["mean_site_div"]))
        ax.axline((x, x), slope=1.0, color="black", linestyle="dashed")
        ax.set_xlabel("Observed diversity per window (site)")
        ax.set_ylabel("Simulated diversity per window (sim)")
        ax.legend()
        div_scatter_path = args.out_dir / f"{args.prefix}diversity-scatter.png"
        plt.savefig(div_scatter_path)
        plt.clf()

    # ------------------------------------------------------------------ #
    # Diversity skyline
    # ------------------------------------------------------------------ #
    fig, ax = plt.subplots(1, 1, figsize=(8, 4), constrained_layout=True)
    if args.sim_branch:
        ax.fill_between(pri["coord"], pri["q_sim_div"][0], pri["q_sim_div"][1],
                        color=PRI_SIM, alpha=0.1)
        ax.scatter(pri["coord"], pri["mean_sim_div"], c=PRI_SIM, s=8,
                   label=f"{pri_label} sim")
    ax.scatter(pri["coord"], pri["mean_site_div"], c=PRI_SITE, s=8,
               label=f"{pri_label} site")
    if cmp is not None:
        if args.sim_branch:
            ax.fill_between(cmp["coord"], cmp["q_sim_div"][0], cmp["q_sim_div"][1],
                            color=CMP_SIM, alpha=0.1)
            ax.scatter(cmp["coord"], cmp["mean_sim_div"], c=CMP_SIM, s=8, marker="^",
                       label=f"{cmp_label} sim")
        ax.scatter(cmp["coord"], cmp["mean_site_div"], c=CMP_SITE, s=8, marker="^",
                   label=f"{cmp_label} site")
    ax.set_xlabel("Position on chromosome")
    ax.set_ylabel("Diversity")
    ax.legend(fontsize=7)
    div_skyline_path = args.out_dir / f"{args.prefix}diversity-skyline.png"
    plt.savefig(div_skyline_path)
    plt.clf()

    # ------------------------------------------------------------------ #
    # Diversity trace across replicates
    # ------------------------------------------------------------------ #
    fig, ax = plt.subplots(1, 1, figsize=(5, 4), constrained_layout=True)
    ax.plot(pri["mcmc_iterates"], pri["trace_site_div"], "-", c=PRI_SITE,
            label=f"{pri_label} site")
    if args.sim_branch:
        ax.plot(pri["mcmc_iterates"], pri["trace_sim_div"], "-", c=PRI_SIM,
                label=f"{pri_label} sim")
    if cmp is not None:
        ax.plot(cmp["mcmc_iterates"], cmp["trace_site_div"], CMP_LS, c=CMP_SITE,
                label=f"{cmp_label} site")
        if args.sim_branch:
            ax.plot(cmp["mcmc_iterates"], cmp["trace_sim_div"], CMP_LS, c=CMP_SIM,
                    label=f"{cmp_label} sim")
    ax.set_xlabel("MCMC iteration")
    ax.set_ylabel("Genome-wide diversity")
    ax.legend()
    div_trace_path = args.out_dir / f"{args.prefix}diversity-trace.png"
    plt.savefig(div_trace_path)
    plt.clf()

    # ------------------------------------------------------------------ #
    # Segregating sites: site-vs-sim scatter (only with --sim-branch)
    # ------------------------------------------------------------------ #
    s_scatter_path = None
    if args.sim_branch:
        fig, ax = plt.subplots(1, 1, figsize=(5, 4), constrained_layout=True)
        ax.scatter(pri["mean_site_s"], pri["mean_sim_s"], c=PRI_SIM, s=8, label=pri_label)
        if cmp is not None:
            ax.scatter(cmp["mean_site_s"], cmp["mean_sim_s"], c=CMP_SIM, s=8,
                       marker="^", label=cmp_label)
        x = float(np.nanmean(pri["mean_site_s"]))
        ax.axline((x, x), slope=1.0, color="black", linestyle="dashed")
        ax.set_xlabel("Observed segregating sites / base (site)")
        ax.set_ylabel("Simulated segregating sites / base (sim)")
        ax.legend()
        s_scatter_path = args.out_dir / f"{args.prefix}segsites-scatter.png"
        plt.savefig(s_scatter_path)
        plt.clf()

    # ------------------------------------------------------------------ #
    # Segregating sites skyline
    # ------------------------------------------------------------------ #
    fig, ax = plt.subplots(1, 1, figsize=(8, 4), constrained_layout=True)
    if args.sim_branch:
        ax.fill_between(pri["coord"], pri["q_sim_s"][0], pri["q_sim_s"][1],
                        color=PRI_SIM, alpha=0.1)
        ax.scatter(pri["coord"], pri["mean_sim_s"], c=PRI_SIM, s=8,
                   label=f"{pri_label} sim")
    ax.scatter(pri["coord"], pri["mean_site_s"], c=PRI_SITE, s=8,
               label=f"{pri_label} site")
    if cmp is not None:
        if args.sim_branch:
            ax.fill_between(cmp["coord"], cmp["q_sim_s"][0], cmp["q_sim_s"][1],
                            color=CMP_SIM, alpha=0.1)
            ax.scatter(cmp["coord"], cmp["mean_sim_s"], c=CMP_SIM, s=8, marker="^",
                       label=f"{cmp_label} sim")
        ax.scatter(cmp["coord"], cmp["mean_site_s"], c=CMP_SITE, s=8, marker="^",
                   label=f"{cmp_label} site")
    ax.set_xlabel("Position on chromosome")
    ax.set_ylabel("Segregating sites / base")
    ax.legend(fontsize=7)
    s_skyline_path = args.out_dir / f"{args.prefix}segsites-skyline.png"
    plt.savefig(s_skyline_path)
    plt.clf()

    # ------------------------------------------------------------------ #
    # Segregating sites trace
    # ------------------------------------------------------------------ #
    fig, ax = plt.subplots(1, 1, figsize=(5, 4), constrained_layout=True)
    ax.plot(pri["mcmc_iterates"], pri["trace_site_s"], "-", c=PRI_SITE,
            label=f"{pri_label} site")
    if args.sim_branch:
        ax.plot(pri["mcmc_iterates"], pri["trace_sim_s"], "-", c=PRI_SIM,
                label=f"{pri_label} sim")
    if cmp is not None:
        ax.plot(cmp["mcmc_iterates"], cmp["trace_site_s"], CMP_LS, c=CMP_SITE,
                label=f"{cmp_label} site")
        if args.sim_branch:
            ax.plot(cmp["mcmc_iterates"], cmp["trace_sim_s"], CMP_LS, c=CMP_SIM,
                    label=f"{cmp_label} sim")
    ax.set_xlabel("MCMC iteration")
    ax.set_ylabel("Genome-wide segregating sites / base")
    ax.legend()
    s_trace_path = args.out_dir / f"{args.prefix}segsites-trace.png"
    plt.savefig(s_trace_path)
    plt.clf()

    # ------------------------------------------------------------------ #
    # Tajima's D: site-vs-sim scatter (only with --sim-branch)
    # ------------------------------------------------------------------ #
    td_scatter_path = None
    if args.sim_branch:
        fig, ax = plt.subplots(1, 1, figsize=(5, 4), constrained_layout=True)
        ax.scatter(pri["mean_site_td"], pri["mean_sim_td"], c=PRI_SIM, s=8, label=pri_label)
        if cmp is not None:
            ax.scatter(cmp["mean_site_td"], cmp["mean_sim_td"], c=CMP_SIM, s=8,
                       marker="^", label=cmp_label)
        x = float(np.nanmean(pri["mean_site_td"]))
        ax.axline((x, x), slope=1.0, color="black", linestyle="dashed")
        ax.set_xlabel("Observed Tajima's D per window (site)")
        ax.set_ylabel("Simulated Tajima's D per window (sim)")
        ax.legend()
        td_scatter_path = args.out_dir / f"{args.prefix}tajima-d-scatter.png"
        plt.savefig(td_scatter_path)
        plt.clf()

    # ------------------------------------------------------------------ #
    # Tajima's D skyline
    # ------------------------------------------------------------------ #
    fig, ax = plt.subplots(1, 1, figsize=(8, 4), constrained_layout=True)
    if args.sim_branch:
        ax.fill_between(pri["coord"], pri["q_sim_td"][0], pri["q_sim_td"][1],
                        color=PRI_SIM, alpha=0.1)
        ax.scatter(pri["coord"], pri["mean_sim_td"], c=PRI_SIM, s=8,
                   label=f"{pri_label} sim")
    ax.scatter(pri["coord"], pri["mean_site_td"], c=PRI_SITE, s=8,
               label=f"{pri_label} site")
    if cmp is not None:
        if args.sim_branch:
            ax.fill_between(cmp["coord"], cmp["q_sim_td"][0], cmp["q_sim_td"][1],
                            color=CMP_SIM, alpha=0.1)
            ax.scatter(cmp["coord"], cmp["mean_sim_td"], c=CMP_SIM, s=8, marker="^",
                       label=f"{cmp_label} sim")
        ax.scatter(cmp["coord"], cmp["mean_site_td"], c=CMP_SITE, s=8, marker="^",
                   label=f"{cmp_label} site")
    ax.set_xlabel("Position on chromosome")
    ax.set_ylabel("Tajima's D")
    ax.legend(fontsize=7)
    td_skyline_path = args.out_dir / f"{args.prefix}tajima-d-skyline.png"
    plt.savefig(td_skyline_path)
    plt.clf()

    # ------------------------------------------------------------------ #
    # Tajima's D trace
    # ------------------------------------------------------------------ #
    fig, ax = plt.subplots(1, 1, figsize=(5, 4), constrained_layout=True)
    ax.plot(pri["mcmc_iterates"], pri["trace_site_td"], "-", c=PRI_SITE,
            label=f"{pri_label} site")
    if args.sim_branch:
        ax.plot(pri["mcmc_iterates"], pri["trace_sim_td"], "-", c=PRI_SIM,
                label=f"{pri_label} sim")
    if cmp is not None:
        ax.plot(cmp["mcmc_iterates"], cmp["trace_site_td"], CMP_LS, c=CMP_SITE,
                label=f"{cmp_label} site")
        if args.sim_branch:
            ax.plot(cmp["mcmc_iterates"], cmp["trace_sim_td"], CMP_LS, c=CMP_SIM,
                    label=f"{cmp_label} sim")
    ax.set_xlabel("MCMC iteration")
    ax.set_ylabel("Genome-wide Tajima's D")
    ax.legend()
    td_trace_path = args.out_dir / f"{args.prefix}tajima-d-trace.png"
    plt.savefig(td_trace_path)
    plt.clf()

    # ------------------------------------------------------------------ #
    # Site SFS (with optional sim overlay) — unfolded and folded
    # ------------------------------------------------------------------ #
    sfs_unfolded_path = None
    sfs_folded_path = None
    for _fold, _xlabel, _suffix in [
        (False, "Derived allele frequency", "unfolded"),
        (True,  "Minor allele frequency",   "folded"),
    ]:
        _site_key = f"mean_site_afs_{'folded' if _fold else 'unfolded'}"
        freq = np.arange(1, pri[_site_key].size)
        fig, ax = plt.subplots(1, 1, figsize=(8, 4), constrained_layout=True)
        if args.sim_branch:
            _sim_key = f"mean_sim_afs_{'folded' if _fold else 'unfolded'}"
            _q_key = f"q_sim_afs_{'folded' if _fold else 'unfolded'}"
            ax.fill_between(freq, pri[_q_key][0, 1:], pri[_q_key][1, 1:],
                            color=PRI_SIM, alpha=0.1)
            ax.scatter(freq, pri[_sim_key][1:], c=PRI_SIM, s=8,
                       label=f"{pri_label} sim")
        ax.scatter(freq, pri[_site_key][1:], c=PRI_SITE, s=8,
                   label=f"{pri_label} site")
        if cmp is not None:
            c_freq = np.arange(1, cmp[_site_key].size)
            if args.sim_branch:
                ax.fill_between(c_freq, cmp[_q_key][0, 1:], cmp[_q_key][1, 1:],
                                color=CMP_SIM, alpha=0.1)
                ax.scatter(c_freq, cmp[_sim_key][1:], c=CMP_SIM, s=8, marker="^",
                           label=f"{cmp_label} sim")
            ax.scatter(c_freq, cmp[_site_key][1:], c=CMP_SITE, s=8, marker="^",
                       label=f"{cmp_label} site")
        ax.set_xlabel(_xlabel)
        ax.set_ylabel("# of variants / base")
        _log_series = [pri[_site_key][1:]]
        if args.sim_branch:
            _log_series.append(pri[_sim_key][1:])
        if cmp is not None:
            _log_series.append(cmp[_site_key][1:])
        safe_log_yscale(ax, *_log_series,
                        context=f"frequency-spectrum-{_suffix}")
        ax.legend(fontsize=7)
        _plot_path = args.out_dir / f"{args.prefix}frequency-spectrum-{_suffix}.png"
        plt.savefig(_plot_path)
        plt.clf()
        if _suffix == "unfolded":
            sfs_unfolded_path = _plot_path
        else:
            sfs_folded_path = _plot_path

    # ------------------------------------------------------------------ #
    # Optional: observed vs simulated density plots (primary dataset only)
    # ------------------------------------------------------------------ #
    sim_pi_density_path = None
    sim_td_density_path = None
    if args.sim is not None:
        sim_pi, sim_td = load_sim_window_stats(args.sim)
        obs_pi = pri["site_div_vals"][:, pri["keep_post"]].reshape(-1)
        obs_td = pri["site_td_vals"][:, pri["keep_post"]].reshape(-1)
        obs_pi = obs_pi[np.isfinite(obs_pi)]
        obs_td = obs_td[np.isfinite(obs_td)]
        sim_pi = sim_pi[np.isfinite(sim_pi)]
        sim_td = sim_td[np.isfinite(sim_td)]

        if obs_pi.size and sim_pi.size:
            fig, ax = plt.subplots(1, 1, figsize=(6, 4), constrained_layout=True)
            bins = 50
            ax.hist(obs_pi, bins=bins, density=True, alpha=0.4, color="black",
                    label=f"observed — {pri_label}")
            ax.hist(sim_pi, bins=bins, density=True, alpha=0.4, color=PRI_SIM,
                    label="simulated")
            ax.set_xlabel("Nucleotide diversity (pi) per window")
            ax.set_ylabel("Density")
            ax.legend()
            sim_pi_density_path = args.out_dir / f"{args.prefix}diversity-density-vs-sim.png"
            plt.savefig(sim_pi_density_path)
            plt.clf()

        if obs_td.size and sim_td.size:
            fig, ax = plt.subplots(1, 1, figsize=(6, 4), constrained_layout=True)
            bins = 50
            ax.hist(obs_td, bins=bins, density=True, alpha=0.4, color="black",
                    label=f"observed — {pri_label}")
            ax.hist(sim_td, bins=bins, density=True, alpha=0.4, color=PRI_SIM,
                    label="simulated")
            ax.set_xlabel("Tajima's D per window")
            ax.set_ylabel("Density")
            ax.legend()
            sim_td_density_path = args.out_dir / f"{args.prefix}tajima-d-density-vs-sim.png"
            plt.savefig(sim_td_density_path)
            plt.clf()

    sim_sfs_unfolded_plot_path = None
    sim_sfs_folded_plot_path = None
    if args.sim_sfs is not None:
        for _fold, _xlabel, _suffix in [
            (False, "Derived allele frequency", "unfolded"),
            (True,  "Minor allele frequency",   "folded"),
        ]:
            _obs_key = f"mean_site_afs_{'folded' if _fold else 'unfolded'}"
            sim_sfs_vals = load_sim_sfs(args.sim_sfs, folded=_fold)
            sim_mean_sfs = safe_nanmean(sim_sfs_vals, axis=-1)
            sim_q_sfs = safe_nanquantile(sim_sfs_vals, [0.025, 0.975], axis=-1)
            n = min(pri[_obs_key].size, sim_mean_sfs.size)
            freq = np.arange(1, n)
            fig, ax = plt.subplots(1, 1, figsize=(8, 4), constrained_layout=True)
            ax.fill_between(freq, sim_q_sfs[0, 1:n], sim_q_sfs[1, 1:n], color="steelblue", alpha=0.15)
            ax.scatter(freq, sim_mean_sfs[1:n], c="steelblue", label="simulated", s=8)
            ax.scatter(freq, pri[_obs_key][1:n], c="black",
                       label=f"observed — {pri_label}", s=8)
            ax.set_xlabel(_xlabel)
            ax.set_ylabel("# of variants / base")
            safe_log_yscale(ax, sim_mean_sfs[1:n], pri[_obs_key][1:n],
                            context=f"frequency-spectrum-vs-sim-{_suffix}")
            ax.legend()
            _plot_path = args.out_dir / f"{args.prefix}frequency-spectrum-vs-sim-{_suffix}.png"
            plt.savefig(_plot_path)
            plt.clf()
            if _suffix == "unfolded":
                sim_sfs_unfolded_plot_path = _plot_path
            else:
                sim_sfs_folded_plot_path = _plot_path

    # ------------------------------------------------------------------ #
    # Summary
    # ------------------------------------------------------------------ #
    summary_path = args.out_dir / f"{args.prefix}summary.txt"
    summary_path.write_text(
        "\n".join([
            f"ts_dir={args.ts_dir}",
            f"pattern={args.pattern}",
            f"n_files={pri['n_files']}",
            f"burnin_frac={args.burnin_frac}",
            f"burnin_index={pri['burnin']}",
            f"mcmc_thin={args.mcmc_thin}",
            f"mcmc_iter_min={int(pri['mcmc_iterates'][0])}",
            f"mcmc_iter_max={int(pri['mcmc_iterates'][-1])}",
            f"window_size={args.window_size}",
            f"mutation_rate={args.mutation_rate}",
            f"sim_branch={args.sim_branch}",
            f"compare_dir={args.compare}",
            f"n_files_compare={cmp['n_files'] if cmp else 'n/a'}",
            f"sim_tsv={args.sim}",
            f"sim_sfs_tsv={args.sim_sfs}",
            f"n_samples={pri['n_samples']}",
            f"sequence_length_min={float(np.min(pri['seq_lengths']))}",
            f"sequence_length_max={float(np.max(pri['seq_lengths']))}",
            f"mutational_load_plot={load_path}",
            f"mutational_load_trace_plot={trace_path}",
            f"diversity_scatter_plot={div_scatter_path}",
            f"diversity_skyline_plot={div_skyline_path}",
            f"diversity_trace_plot={div_trace_path}",
            f"tajima_d_scatter_plot={td_scatter_path}",
            f"tajima_d_skyline_plot={td_skyline_path}",
            f"tajima_d_trace_plot={td_trace_path}",
            f"segsites_scatter_plot={s_scatter_path}",
            f"segsites_skyline_plot={s_skyline_path}",
            f"segsites_trace_plot={s_trace_path}",
            f"frequency_spectrum_unfolded_plot={sfs_unfolded_path}",
            f"frequency_spectrum_folded_plot={sfs_folded_path}",
            f"frequency_spectrum_vs_sim_unfolded_plot={sim_sfs_unfolded_plot_path}",
            f"frequency_spectrum_vs_sim_folded_plot={sim_sfs_folded_plot_path}",
            f"diversity_density_vs_sim_plot={sim_pi_density_path}",
            f"tajima_d_density_vs_sim_plot={sim_td_density_path}",
        ]) + "\n"
    )

    for p in [load_path, trace_path, div_skyline_path, div_trace_path,
              td_skyline_path, td_trace_path,
              s_skyline_path, s_trace_path,
              sfs_unfolded_path, sfs_folded_path]:
        print(f"Wrote: {p}")
    for p in [div_scatter_path, td_scatter_path, s_scatter_path,
              sim_sfs_unfolded_plot_path, sim_sfs_folded_plot_path,
              sim_pi_density_path, sim_td_density_path]:
        if p is not None:
            print(f"Wrote: {p}")
    print(f"Wrote: {summary_path}")


if __name__ == "__main__":
    main()
