#!/usr/bin/env python
from __future__ import annotations

import argparse
import html
import sys
from pathlib import Path

import numpy as np

import msprime

from argtest_common import (
    aggregate_by_individual,
    load_ts,
    mutational_load,
    resolve_mu_rate,
    sample_names,
)


def ascii_bar_html(value, max_value, width=40, color="#777777"):
    n = max(0, int(round(value / max_value * width))) if max_value > 0 else 0
    return f'<span style="color:{color}">{"█" * n}</span>'


def load_chart_html(load_1d, names, title, outlier_mask=None):
    max_val = float(np.max(load_1d)) if load_1d.size else 1.0
    label_width = max((len(n) for n in names), default=1)
    rows = []
    for i, (name, val) in enumerate(zip(names, load_1d)):
        is_out = outlier_mask is not None and bool(outlier_mask[i])
        color = "#d62728" if is_out else "#444444"
        bar = ascii_bar_html(val, max_val, color=color)
        label = html.escape(name.ljust(label_width))
        rows.append(
            f'<span style="color:{color}">{label}</span>  {bar}  '
            f'<span style="color:{color}">{val:.4f}</span>'
        )
    inner = "\n".join(rows)
    return (
        f"<h2>{html.escape(title)}</h2>\n"
        f'<pre style="line-height:1.4;font-size:12px">{inner}</pre>\n'
    )


def outlier_hist_html(counts):
    if not counts:
        return ""
    max_count = max(counts)
    bins: dict[int, int] = {}
    for c in counts:
        bins[c] = bins.get(c, 0) + 1
    max_freq = max(bins.values()) if bins else 1
    rows = []
    for k in range(0, max_count + 1):
        freq = bins.get(k, 0)
        bar_len = int(round(freq / max_freq * 40)) if max_freq > 0 else 0
        rows.append(f"{k:>4}  {'█' * bar_len:<42} {freq}")
    inner = "\n".join(rows)
    return (
        "<h2>Outlier window counts per individual</h2>\n"
        f'<pre style="line-height:1.4;font-size:12px">{inner}</pre>\n'
    )


def lineage_label(name: str) -> str:
    for sep in ("_", "-", "."):
        if sep in name:
            head = name.split(sep, 1)[0].strip()
            if head:
                return head
    return name


def summarize_lineage_outliers(names, outlier_counts):
    summary = {}
    for name, count in zip(names, outlier_counts):
        lineage = lineage_label(name)
        summary[lineage] = summary.get(lineage, 0) + int(count)
    return sorted(summary.items(), key=lambda item: (-item[1], item[0]))


def lineage_outlier_table_html(names, outlier_counts):
    rows = summarize_lineage_outliers(names, outlier_counts)
    body = "\n".join(
        f"<tr><td>{html.escape(lineage)}</td><td>{count}</td></tr>"
        for lineage, count in rows
    )
    return (
        "<h2>Outlier windows by lineage</h2>\n"
        "<table>\n"
        "<thead><tr><th>Lineage</th><th>Outlier windows</th></tr></thead>\n"
        f"<tbody>\n{body}\n</tbody>\n"
        "</table>\n"
    )


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
        default=0.5,
        help=(
            "Outlier cutoff as a fraction of the per-individual sim-based expected "
            "load (default: 0.5)."
        ),
    )
    p.add_argument(
        "--mutation-rate",
        type=float,
        default=None,
        help=(
            "Scalar mutation-rate fallback when neither ts.metadata nor a sibling "
            "*.mut_rate.p file provides a ratemap."
        ),
    )
    p.add_argument(
        "--random-seed",
        type=int,
        default=1,
        help="Seed for msprime.sim_mutations when computing expected load (default: 1).",
    )
    p.add_argument(
        "--out",
        default="mutational_load_summary.html",
        help=(
            "Output filename (default: mutational_load_summary.html). Only the "
            "filename part is used; the file is always written to "
            "<repo-root>/results/<filename>."
        ),
    )
    p.add_argument(
        "--suffix-to-strip",
        default="",
        help='Suffix stripped from sample names before display (default: "").',
    )
    return p.parse_args()


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
    edges = positions[snp_window::snp_window]
    return np.concatenate((np.array([0.0]), edges, np.array([L])))


def simulate_expected_load(ts, ts_path, windows, names, scalar_rate, seed):
    # Single-sim per-individual expected load matrix (windows x individuals).
    mu = resolve_mu_rate(ts, ts_path, scalar_fallback=scalar_rate)
    sim_ts = msprime.sim_mutations(ts, rate=mu, keep=False, random_seed=seed)
    expected = mutational_load(sim_ts, windows=windows)
    expected, _ = aggregate_by_individual(expected, names)
    return expected


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

            ts = load_ts(ts_path)

            windows = None
            if args.window_size is not None:
                windows = build_bp_windows(ts, args.window_size)
            elif args.snp_window is not None:
                windows = build_snp_windows(ts, args.snp_window)

            load = mutational_load(ts, windows=windows)
            names = sample_names(ts, suffix_to_strip=args.suffix_to_strip)
            load, unique_names = aggregate_by_individual(load, names)

            body_parts = []
            meta_lines = [
                f"Input: {html.escape(ts_path.name)} | "
                f"Samples: {ts.num_samples} | "
                f"Individuals: {len(unique_names)}"
            ]

            if windows is None:
                body_parts.append(load_chart_html(load, unique_names, "Mutational load"))
            else:
                expected = simulate_expected_load(
                    ts, ts_path, windows, names,
                    scalar_rate=args.mutation_rate,
                    seed=args.random_seed,
                )
                valid = expected > 0
                high = (1 + args.cutoff) * expected
                low = (1 - args.cutoff) * expected
                mask = ((load > high) | (load < low)) & valid
                outlier_counts = mask.sum(axis=0).astype(int).tolist()

                mean_load = load.mean(axis=0)
                any_outlier = mask.any(axis=0)
                body_parts.append(
                    load_chart_html(
                        mean_load, unique_names,
                        "Mean mutational load per individual",
                        outlier_mask=any_outlier,
                    )
                )

                outlier_window_count = int(mask.any(axis=1).sum())
                total_window_count = int(load.shape[0])
                meta_lines.append(
                    f"{outlier_window_count} of {total_window_count} windows "
                    f"have at least one outlier"
                )

                body_parts.append(outlier_hist_html(outlier_counts))
                body_parts.append(lineage_outlier_table_html(unique_names, outlier_counts))

                if args.window_size is not None:
                    meta_lines.append(
                        f"Window size: {int(args.window_size)} bp | "
                        f"Outlier cutoff: {args.cutoff:.3f} of sim expectation"
                    )
                elif args.snp_window is not None:
                    meta_lines.append(
                        f"Window size: {int(args.snp_window)} variants | "
                        f"Outlier cutoff: {args.cutoff:.3f} of sim expectation"
                    )

            meta_html = "\n".join(
                f'<div class="meta">{m}</div>' for m in meta_lines
            )

            html_doc = f"""<!doctype html>
<html lang="en">
<meta charset="utf-8">
<title>Mutational load summary</title>
<style>
body {{ font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif; margin: 24px; }}
h1 {{ font-size: 20px; }}
h2 {{ font-size: 16px; margin-top: 24px; }}
.meta {{ color: #444; font-size: 13px; margin-bottom: 8px; }}
pre {{ background: #f8f8f8; padding: 12px; border-radius: 4px; overflow-x: auto; }}
table {{ border-collapse: collapse; margin-bottom: 20px; }}
th, td {{ border: 1px solid #cccccc; padding: 6px 10px; text-align: left; }}
</style>
<h1>Mutational load summary</h1>
{meta_html}
{"".join(body_parts)}</html>
"""
            out.write_text(html_doc)

        finally:
            sys.stdout = old_stdout
            sys.stderr = old_stderr


if __name__ == "__main__":
    main()
