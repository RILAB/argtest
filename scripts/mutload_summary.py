#!/usr/bin/env python
from __future__ import annotations

import argparse
import base64
from io import BytesIO
import os
import sys
from pathlib import Path

import numpy as np

plt = None

from argtest_common import (
    aggregate_by_individual,
    load_ts,
    mutational_load,
    sample_names,
)


def fig_to_data_url(fig) -> str:
    # Encode a matplotlib figure as an inline SVG data URL so plots stay sharp.
    buf = BytesIO()
    with plt.rc_context({"svg.hashsalt": "argtest"}):
        fig.savefig(
            buf,
            format="svg",
            bbox_inches="tight",
            metadata={"Date": "1970-01-01T00:00:00"},
        )
    plt.close(fig)
    buf.seek(0)
    svg_b64 = base64.b64encode(buf.read()).decode("ascii")
    return f"data:image/svg+xml;base64,{svg_b64}"


def plot_single(load, names, title):
    # Simple per-individual bar chart.
    fig, ax = plt.subplots(figsize=(max(6, 0.35 * len(names)), 4))
    ax.bar(names, load, color="#777777")
    ax.set_ylabel("Derived mutational load")
    ax.set_title(title)
    ax.tick_params(axis="x", labelrotation=90, labelsize=8)
    fig.tight_layout()
    return fig


def lineage_label(name: str) -> str:
    # Use the leading token as the lineage label when names encode subgroups.
    for sep in ("_", "-", "."):
        if sep in name:
            head = name.split(sep, 1)[0].strip()
            if head:
                return head
    return name


def summarize_lineage_outliers(names, outlier_counts):
    # Sum outlier windows across individuals sharing the same lineage token.
    summary = {}
    for name, count in zip(names, outlier_counts):
        lineage = lineage_label(name)
        summary[lineage] = summary.get(lineage, 0) + int(count)
    return sorted(summary.items(), key=lambda item: (-item[1], item[0]))


def lineage_outlier_table_html(names, outlier_counts):
    rows = summarize_lineage_outliers(names, outlier_counts)
    body = "\n".join(
        f"<tr><td>{lineage}</td><td>{count}</td></tr>"
        for lineage, count in rows
    )
    return (
        "<h2>Outlier windows by lineage</h2>\n"
        "<table>\n"
        "<thead><tr><th>Lineage</th><th>Outlier windows</th></tr></thead>\n"
        f"<tbody>\n{body}\n</tbody>\n"
        "</table>\n"
    )


def plot_windows(load, names, windows, outlier_mask=None):
    # Grid of per-window bar charts with shared y-scale.
    nwin = load.shape[0]
    ncols = 4
    nrows = int(np.ceil(nwin / ncols))
    # Use a shared y-axis scale across all window plots
    ymax = float(load.max()) if load.size else 0.0
    ytop = ymax * 1.05 if ymax > 0 else 1.0
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * 4, nrows * 3), squeeze=False)
    for i in range(nwin):
        r, c = divmod(i, ncols)
        ax = axes[r][c]
        colors = ["#777777"] * len(names)
        if outlier_mask is not None:
            colors = ["#d62728" if is_outlier else "#777777" for is_outlier in outlier_mask[i]]
        ax.bar(names, load[i], color=colors)
        ax.set_ylim(0, ytop)
        left = int(windows[i])
        right = int(windows[i + 1])
        ax.set_title(f"{left:,}-{right:,} bp")
        ax.tick_params(axis="x", labelrotation=90, labelsize=6)
    for j in range(nwin, nrows * ncols):
        r, c = divmod(j, ncols)
        axes[r][c].axis("off")
    fig.tight_layout()
    return fig


def plot_outlier_hist(counts):
    # Histogram of how many outlier windows each individual has.
    counts = [int(c) for c in counts]
    max_count = max(counts) if counts else 0
    bins = list(range(0, max_count + 2))
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(counts, bins=bins, color="#777777", edgecolor="#222222")
    ax.set_xlabel("Outlier windows per individual")
    ax.set_ylabel("Count of individuals")
    ax.set_title("Outlier window counts")
    if max_count <= 4:
        ticks = bins
    else:
        # Keep the histogram readable for large window counts.
        ticks = np.unique(np.rint(np.linspace(0, max_count + 1, 5)).astype(int)).tolist()
    ax.set_xticks(ticks)
    fig.tight_layout()
    return fig


def parse_args():
    p = argparse.ArgumentParser(
        description="Summarize mutational load per individual from a tree sequence",
    )
    p.add_argument("ts", help="Tree sequence file (.ts, .trees, or .tsz)")
    group = p.add_mutually_exclusive_group()
    group.add_argument("--window-size", type=float, help="Window size in bp")
    group.add_argument(
        "--snp-window",
        type=int,
        help="Window size as a fixed number of variants",
    )
    p.add_argument(
        "--cutoff",
        type=float,
        default=0.25,
        help="Outlier cutoff as a fraction of the window median (default: 0.25)",
    )
    p.add_argument(
        "--fraction",
        type=float,
        default=None,
        help="Write windows where the fraction of outlier individuals is greater than this value",
    )
    p.add_argument("--out", default="mutational_load_summary.html")
    p.add_argument("--suffix-to-strip", default="_anchorwave")
    args = p.parse_args()
    if args.fraction is not None and not 0 <= args.fraction <= 1:
        raise ValueError("--fraction must be between 0 and 1")
    return args


def build_bp_windows(ts, window_size: float) -> np.ndarray:
    if window_size <= 0:
        raise ValueError("--window-size must be > 0")
    L = float(ts.sequence_length)
    windows = np.arange(0, L + window_size, window_size, dtype=float)
    if windows[-1] > L:
        windows[-1] = L
    return windows


def build_snp_windows(ts, snp_window: int) -> np.ndarray:
    if snp_window <= 0:
        raise ValueError("--snp-window must be > 0")
    positions = np.asarray(ts.sites_position, dtype=float)
    L = float(ts.sequence_length)
    if positions.size == 0:
        return np.array([0.0, L], dtype=float)
    # Use the position of each (N+1)th variant as the next window edge so each
    # preceding window contains exactly N variants when site positions are unique.
    edges = positions[snp_window::snp_window]
    return np.concatenate((np.array([0.0]), edges, np.array([L])))


class Tee:
    def __init__(self, *streams):
        self.streams = streams

    def write(self, data):
        for s in self.streams:
            s.write(data)
            s.flush()

    def flush(self):
        for s in self.streams:
            s.flush()


def main():
    args = parse_args()
    ts_path = Path(args.ts)
    repo_root = Path(__file__).resolve().parent.parent
    results_dir = repo_root / "results"
    logs_dir = repo_root / "logs"
    results_dir.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)

    out = results_dir / Path(args.out).name

    log_path = logs_dir / f"{out.stem}.log"
    old_stdout = sys.stdout
    old_stderr = sys.stderr
    with open(log_path, "w") as log_fh:
        try:
            sys.stdout = Tee(old_stdout, log_fh)
            sys.stderr = Tee(old_stderr, log_fh)
            # Ensure matplotlib cache is writable and avoid stderr warnings
            os.environ.setdefault("MPLCONFIGDIR", str(logs_dir / ".matplotlib"))
            global plt  # pylint: disable=global-statement
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt  # noqa: E402

            ts = load_ts(ts_path)

            windows = None
            if args.window_size is not None:
                windows = build_bp_windows(ts, args.window_size)
            elif args.snp_window is not None:
                windows = build_snp_windows(ts, args.snp_window)

            load = mutational_load(ts, windows=windows)
            names = sample_names(ts, suffix_to_strip=args.suffix_to_strip)
            load, unique_names = aggregate_by_individual(load, names)

            outlier_counts = None
            mask = None
            if windows is None:
                fig = plot_single(load, unique_names, "Mutational load")
            else:
                window_medians = np.median(load, axis=1)
                valid = window_medians > 0
                high = (1 + args.cutoff) * window_medians
                low = (1 - args.cutoff) * window_medians
                # Identify per-individual outliers in each window.
                mask = (load > high[:, None]) | (load < low[:, None])
                mask &= valid[:, None]
                outlier_counts = mask.sum(axis=0).astype(int).tolist()
                fig = plot_windows(load, unique_names, windows, outlier_mask=mask)

            img_url = fig_to_data_url(fig)
            hist_html = ""
            lineage_html = ""
            if outlier_counts is not None:
                outlier_window_count = int(mask.any(axis=1).sum())
                total_window_count = int(load.shape[0])
                hist_fig = plot_outlier_hist(outlier_counts)
                hist_url = fig_to_data_url(hist_fig)
                hist_html = (
                    "<h2>Outlier window counts</h2>\n"
                    f"<div class=\"meta\">{outlier_window_count} total windows with at least one outlier out of {total_window_count} total windows</div>\n"
                    f"<img src=\"{hist_url}\" alt=\"Outlier window counts histogram\">\n"
                )
                lineage_html = lineage_outlier_table_html(unique_names, outlier_counts)
            html = f"""<!doctype html>
<html lang="en">
<meta charset="utf-8">
<title>Mutational load summary</title>
<style>
body {{ font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif; margin: 24px; }}
h1 {{ font-size: 20px; }}
.meta {{ color: #444; font-size: 13px; margin-bottom: 16px; }}
img {{ max-width: 100%; height: auto; }}
table {{ border-collapse: collapse; margin-bottom: 20px; }}
th, td {{ border: 1px solid #cccccc; padding: 6px 10px; text-align: left; }}
</style>
<h1>Mutational load summary</h1>
{hist_html}
{lineage_html}
<div class="meta">Input: {ts_path.name} | Samples: {ts.num_samples} | Individuals: {len(unique_names)}</div>
<img src="{img_url}" alt="Mutational load plot">
"""

            if windows is not None and args.window_size is not None:
                html += f"<div class=\"meta\">Window size: {int(args.window_size)} bp</div>"
                html += f"<div class=\"meta\">Outlier cutoff: {args.cutoff:.3f} of window median</div>"
            elif windows is not None and args.snp_window is not None:
                html += f"<div class=\"meta\">Window size: {int(args.snp_window)} variants</div>"
                html += f"<div class=\"meta\">Outlier cutoff: {args.cutoff:.3f} of window median</div>"

            html += "</html>\n"
            out.write_text(html)

            if windows is not None:
                # Write one BED listing outliers per window.
                out_path = results_dir / f"{ts_path.stem}_outliers.bed"
                lines = []
                for w in range(load.shape[0]):
                    if not valid[w]:
                        continue
                    row_mask = mask[w]
                    if not row_mask.any():
                        continue
                    outlier_names = [unique_names[i] for i in np.where(row_mask)[0]]
                    outlier_vals = [f"{load[w, i]:.3f}" for i in np.where(row_mask)[0]]
                    start = int(windows[w])
                    end = int(windows[w + 1])
                    lines.append(
                        f"{ts_path.stem}\t{start}\t{end}\t{','.join(outlier_names)}\t{','.join(outlier_vals)}\t{window_medians[w]:.3f}"
                    )
                out_path.write_text("\n".join(lines) + ("\n" if lines else ""))

                if args.fraction is not None:
                    masked_path = results_dir / f"{ts_path.stem}_mutation_masked.bed"
                    masked_lines = []
                    outlier_fractions = mask.sum(axis=1) / load.shape[1]
                    for w in range(load.shape[0]):
                        if not valid[w] or outlier_fractions[w] <= args.fraction:
                            continue
                        start = int(windows[w])
                        end = int(windows[w + 1])
                        masked_lines.append(
                            f"{ts_path.stem}\t{start}\t{end}\t{outlier_fractions[w]:.3f}\t{int(mask[w].sum())}\t{load.shape[1]}"
                        )
                    masked_path.write_text("\n".join(masked_lines) + ("\n" if masked_lines else ""))
        finally:
            sys.stdout = old_stdout
            sys.stderr = old_stderr


if __name__ == "__main__":
    main()
