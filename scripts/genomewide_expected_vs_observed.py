#!/usr/bin/env python
"""Genome-wide expected-vs-observed scatter plots, coloured by recombination rate.

The per-chromosome scatters written by validation_plots_from_ts.py compare each
window's observed (site-mode) statistic to its expectation (mutations simulated
on the inferred ARG), one chromosome at a time. This script pools every window
from every chromosome into one panel per statistic and colours each point by the
recombination rate of that window, so a rate-dependent bias -- the signature of
ARG inference struggling in low-recombination regions -- shows up as colour
structure around the 1:1 line rather than as a difference between chromosomes.

It reads the windows.tsv dumps rather than the tree sequences, so it is cheap
and never needs to reload an ARG.
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
os.environ.setdefault("XDG_CACHE_HOME", "/tmp")

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm, Normalize

from hapmap_low_rec_mask import (
    _resolve_chrom,
    build_intervals,
    load_fai,
    load_hapmap,
)

matplotlib.rcParams["figure.dpi"] = 300

# (file stem, windows.tsv column prefix, axis label)
STATS = [
    ("diversity", "pi", "Nucleotide diversity (pi) per window"),
    ("tajima-d", "tajimas_d", "Tajima's D per window"),
]


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description=(
            "Pool per-window statistics across chromosomes and plot expected vs "
            "observed, coloured by the recombination rate of each window."
        )
    )
    p.add_argument(
        "--windows",
        type=Path,
        nargs="+",
        required=True,
        help=(
            "Per-chromosome windows.tsv files written by validation_plots_from_ts.py. "
            "The chromosome name is taken from the grandparent directory, matching "
            "the step-6 layout <step6_dir>/<chrom>/<variant>/windows.tsv."
        ),
    )
    p.add_argument(
        "--hapmap",
        type=Path,
        required=True,
        help="HapMap-format recombination map (same file step 1 uses).",
    )
    p.add_argument(
        "--fai",
        type=Path,
        required=True,
        help="FASTA index supplying chromosome lengths, to close the last map interval.",
    )
    p.add_argument(
        "--out-dir",
        type=Path,
        required=True,
        help="Directory for the pooled TSV and the PNGs.",
    )
    p.add_argument(
        "--prefix",
        default="",
        help="Optional filename prefix for the outputs.",
    )
    return p.parse_args(argv)


def read_windows_tsv(path: Path) -> tuple[str, list[dict]]:
    """Read one windows.tsv, keeping only its first dataset.

    Returns (dataset_label, rows). A windows.tsv can hold a second dataset when
    validation_plots_from_ts.py ran with --compare; pooling two datasets into one
    scatter would silently mix them, so only the first is used and the rest are
    reported to stderr.
    """
    lines = path.read_text().splitlines()
    if not lines:
        raise ValueError(f"{path} is empty")
    header = lines[0].split("\t")
    rows = []
    label = None
    skipped = set()
    for line in lines[1:]:
        if not line.strip():
            continue
        row = dict(zip(header, line.split("\t")))
        if label is None:
            label = row["dataset"]
        if row["dataset"] != label:
            skipped.add(row["dataset"])
            continue
        rows.append(row)
    if skipped:
        print(
            f"NOTE: {path} also holds dataset(s) {sorted(skipped)}; pooling only "
            f"{label!r}.",
            file=sys.stderr,
        )
    return label, rows


def _cell(row: dict, key: str) -> float:
    """A windows.tsv cell as a float; blank cells are gaps, i.e. NaN."""
    text = row.get(key, "")
    return float(text) if text else np.nan


def mean_rate_per_window(
    starts: np.ndarray, ends: np.ndarray, intervals: list[tuple[int, int, float]]
) -> np.ndarray:
    """Length-weighted mean map rate over each [start, end) window.

    Windows with no overlapping map interval (or zero total overlap) come back as
    NaN rather than 0, so "no map here" is not plotted as "no recombination here".
    """
    if not intervals:
        return np.full(starts.shape, np.nan)
    iv = np.asarray(intervals, dtype=float)
    order = np.argsort(iv[:, 0])
    iv = iv[order]
    ivs, ive, ivr = iv[:, 0], iv[:, 1], iv[:, 2]

    out = np.full(starts.shape, np.nan)
    for i, (ws, we) in enumerate(zip(starts, ends)):
        if not (np.isfinite(ws) and np.isfinite(we)) or we <= ws:
            continue
        overlap = np.clip(np.minimum(ive, we) - np.maximum(ivs, ws), 0.0, None)
        total = overlap.sum()
        if total > 0:
            out[i] = float((overlap * ivr).sum() / total)
    return out


def collect(windows_paths, hapmap, fai) -> list[dict]:
    """Pool every window from every chromosome, annotated with its map rate."""
    pooled = []
    for path in sorted(windows_paths, key=lambda p: p.parent.parent.name):
        chrom = path.parent.parent.name
        label, rows = read_windows_tsv(path)
        if not rows:
            print(f"NOTE: {path} has no data rows; skipping.", file=sys.stderr)
            continue

        hap_key = _resolve_chrom(chrom, hapmap)
        fai_key = _resolve_chrom(chrom, fai)
        if hap_key is None or fai_key is None:
            print(
                f"WARNING: no recombination map ({hap_key is None}) or FAI length "
                f"({fai_key is None}) for {chrom}; its windows are dropped.",
                file=sys.stderr,
            )
            continue
        intervals = build_intervals(hapmap[hap_key], fai[fai_key])

        starts = np.array([_cell(r, "window_start") for r in rows])
        ends = np.array([_cell(r, "window_end") for r in rows])
        rates = mean_rate_per_window(starts, ends, intervals)
        for row, start, end, rate in zip(rows, starts, ends, rates):
            entry = {
                "chrom": chrom,
                "dataset": label,
                "window_start": start,
                "window_end": end,
                "rec_rate_cm_per_mb": rate,
            }
            for _stem, col, _lbl in STATS:
                entry[f"obs_{col}"] = _cell(row, f"obs_{col}")
                entry[f"exp_{col}"] = _cell(row, f"exp_{col}")
            pooled.append(entry)
    return pooled


def write_pooled_tsv(path: Path, pooled: list[dict]) -> None:
    """The data behind the plots below, in the same layout as windows.tsv."""
    cols = ["chrom", "dataset", "window_start", "window_end", "rec_rate_cm_per_mb"]
    for _stem, col, _lbl in STATS:
        cols += [f"obs_{col}", f"exp_{col}"]
    lines = ["\t".join(cols)]
    for entry in pooled:
        cells = []
        for col in cols:
            val = entry[col]
            if isinstance(val, str):
                cells.append(val)
            else:
                cells.append("" if not np.isfinite(val) else repr(float(val)))
        lines.append("\t".join(cells))
    path.write_text("\n".join(lines) + "\n")


def rate_norm(rates: np.ndarray):
    """Colour normalisation for map rates, plus the rates as coloured.

    Recombination rates usually span orders of magnitude, so a log scale is used
    when the positive rates cover more than a decade. Zero-rate windows have no
    place on a log axis; they are clipped to the smallest positive rate and the
    colourbar label says so, which is honest about the floor without dropping
    the windows that matter most for a low-recombination bias.
    """
    positive = rates[np.isfinite(rates) & (rates > 0)]
    if positive.size and positive.max() / positive.min() > 10:
        floor = float(positive.min())
        clipped = np.where(np.isfinite(rates) & (rates <= 0), floor, rates)
        label = "Recombination rate (cM/Mb, log scale)"
        if np.any(np.isfinite(rates) & (rates <= 0)):
            label += f"\nzero-rate windows clipped to {floor:.3g}"
        return clipped, LogNorm(vmin=floor, vmax=float(positive.max())), label
    finite = rates[np.isfinite(rates)]
    vmin = float(finite.min()) if finite.size else 0.0
    vmax = float(finite.max()) if finite.size else 1.0
    if vmax <= vmin:
        vmax = vmin + 1.0
    return rates, Normalize(vmin=vmin, vmax=vmax), "Recombination rate (cM/Mb)"


def scatter_plot(path: Path, obs, exp, rates, axis_label: str, n_chroms: int) -> bool:
    """One expected-vs-observed panel. Returns False if there was nothing to draw."""
    keep = np.isfinite(obs) & np.isfinite(exp)
    if not keep.any():
        return False
    obs, exp, rates = obs[keep], exp[keep], rates[keep]
    colors, norm, cbar_label = rate_norm(rates)

    fig, ax = plt.subplots(1, 1, figsize=(5.5, 4.5), constrained_layout=True)
    # Windows with no map coverage still belong on the plot; draw them as
    # uncoloured background so the point count matches the pooled TSV.
    grey = ~np.isfinite(colors)
    if grey.any():
        ax.scatter(obs[grey], exp[grey], c="lightgray", s=6, linewidths=0,
                   label="no map coverage")
    sc = ax.scatter(obs[~grey], exp[~grey], c=colors[~grey], norm=norm,
                    cmap="viridis", s=6, linewidths=0)
    lo = float(min(np.min(obs), np.min(exp)))
    hi = float(max(np.max(obs), np.max(exp)))
    ax.plot([lo, hi], [lo, hi], color="black", linestyle="dashed", linewidth=1,
            label="1:1")
    ax.set_xlabel(f"Observed — {axis_label}")
    ax.set_ylabel(f"Expected (sim) — {axis_label}")
    ax.set_title(f"{obs.size} windows pooled across {n_chroms} chromosomes",
                 fontsize=9)
    ax.legend(fontsize=7)
    fig.colorbar(sc, ax=ax, label=cbar_label)
    fig.savefig(path)
    plt.close(fig)
    return True


def main():
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    hapmap = load_hapmap(args.hapmap)
    fai = load_fai(args.fai)
    pooled = collect(args.windows, hapmap, fai)
    if not pooled:
        raise RuntimeError(
            "No windows could be pooled: none of the windows.tsv inputs had rows "
            "with a matching recombination map and FAI entry."
        )

    tsv_path = args.out_dir / f"{args.prefix}genomewide-windows.tsv"
    write_pooled_tsv(tsv_path, pooled)
    print(f"Wrote: {tsv_path}")

    n_chroms = len({e["chrom"] for e in pooled})
    rates = np.array([e["rec_rate_cm_per_mb"] for e in pooled])
    for stem, col, axis_label in STATS:
        obs = np.array([e[f"obs_{col}"] for e in pooled])
        exp = np.array([e[f"exp_{col}"] for e in pooled])
        png = args.out_dir / f"{args.prefix}{stem}-expected-vs-observed-by-rec.png"
        if scatter_plot(png, obs, exp, rates, axis_label, n_chroms):
            print(f"Wrote: {png}")
        else:
            # Without --sim-branch there is no expectation to plot against, so
            # this is a normal outcome, not a failure.
            print(
                f"NOTE: no window had both an observed and an expected {col}; "
                f"skipped {png.name} (was the pipeline run without sim_branch?)"
            )


if __name__ == "__main__":
    main()
