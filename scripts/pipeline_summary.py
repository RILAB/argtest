#!/usr/bin/env python
"""Generate a self-contained HTML pipeline summary from Snakemake outputs."""
from __future__ import annotations

import argparse
import base64
import math
import re
import subprocess
from collections import defaultdict
from pathlib import Path

from argtest_common import (
    describe_mu_source,
    accessible_intervals_from_mu,
    load_ts,
    ratemap_from_metadata,
    tree_covered_accessible_bp,
)


# ---------------------------------------------------------------------------
# Version stamp
# ---------------------------------------------------------------------------

def git_version() -> str:
    """Return the argtest 'tag (shorthash[-dirty])' stamp, or 'unknown'.

    Derived from the git repo this script lives in, so the summary records
    exactly which pipeline version produced it. Fails soft to 'unknown' when
    git is unavailable or the tree is not a checkout (e.g. exported source).
    """
    repo = Path(__file__).resolve().parent.parent

    def _git(*a: str) -> str:
        return subprocess.run(
            ["git", "-C", str(repo), *a],
            capture_output=True, text=True, check=True,
        ).stdout.strip()

    try:
        commit = _git("rev-parse", "--short", "HEAD")
    except (subprocess.CalledProcessError, OSError):
        return "unknown"

    try:
        tag = _git("describe", "--tags", "--abbrev=0")
    except subprocess.CalledProcessError:
        tag = ""

    dirty = ""
    try:
        if _git("status", "--porcelain"):
            dirty = "-dirty"
    except subprocess.CalledProcessError:
        pass

    return f"{tag} ({commit}{dirty})" if tag else f"{commit}{dirty}"


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Generate a self-contained HTML pipeline summary (run info, "
            "per-chromosome genome retention table, per-individual outlier "
            "counts, and embedded validation plots)."
        ),
    )
    p.add_argument(
        "--out-dir",
        required=True,
        type=Path,
        help=(
            "Snakemake output root containing step1..step6 subdirectories. "
            "Used to locate mask BEDs, trimmed treefiles, and validation PNGs."
        ),
    )
    p.add_argument(
        "--fai",
        required=True,
        type=Path,
        help="FASTA index (.fai) supplying per-chromosome lengths.",
    )
    p.add_argument(
        "--chroms",
        required=True,
        nargs="+",
        help="Chromosome names to include in the report (space-separated).",
    )
    p.add_argument(
        "--replicates",
        required=True,
        nargs="+",
        help="Replicate IDs to aggregate across (space-separated).",
    )
    p.add_argument(
        "--out",
        required=True,
        type=Path,
        help="Output HTML path (typically <out-dir>/pipeline_summary.html).",
    )
    p.add_argument(
        "--filtered-ts",
        required=True,
        nargs="+",
        type=Path,
        help=(
            "Final filtered per-chromosome tree sequences. Paths must use the "
            "<chrom>/<rep>.<suffix> layout produced by steps 5/5b."
        ),
    )
    p.add_argument(
        "--rec-fraction",
        type=float,
        default=None,
        help="Step 1 --rec-fraction value, echoed into the Run-info block.",
    )
    p.add_argument(
        "--low-access-window",
        type=int,
        default=None,
        help="Step 2 --window-size value, echoed into the Run-info block.",
    )
    p.add_argument(
        "--low-access-cutoff",
        type=int,
        default=None,
        help="Step 2 --cutoff-bp value, echoed into the Run-info block.",
    )
    p.add_argument(
        "--mutload-cutoff",
        type=float,
        default=None,
        help="Step 3 --cutoff value, echoed into the Run-info block.",
    )
    p.add_argument(
        "--mutation-rate",
        default=None,
        help="Validation-plot mutation rate, echoed into the Run-info block.",
    )
    p.add_argument(
        "--sim-branch",
        default="false",
        help=(
            'Whether step 6 used --sim-branch. Accepts "true"/"false" as a '
            'string (default: "false"); echoed into the Run-info block.'
        ),
    )
    return p.parse_args()


# ---------------------------------------------------------------------------
# BED / FAI helpers
# ---------------------------------------------------------------------------

def read_fai(path: Path) -> dict[str, int]:
    lengths = {}
    for line in path.read_text().splitlines():
        if not line:
            continue
        parts = line.split("\t")
        lengths[parts[0]] = int(parts[1])
    return lengths


def _fai_num_index(chrom_lengths: dict[str, int]) -> dict[str, int]:
    """Map trailing-digit suffix → length, skipping ambiguous entries."""
    by_num: dict[str, int | None] = {}
    for name, length in chrom_lengths.items():
        m = re.search(r"(\d+)$", name)
        if m:
            num = m.group(1)
            by_num[num] = None if num in by_num else length
    return {k: v for k, v in by_num.items() if v is not None}


def lookup_chrom_len(chrom: str, chrom_lengths: dict[str, int], by_num: dict[str, int]) -> int:
    """Look up chrom length; fall back to numeric-suffix match if name not in FAI."""
    v = chrom_lengths.get(chrom)
    if v is not None:
        return v
    m = re.search(r"(\d+)$", chrom)
    if m:
        return by_num.get(m.group(1), 0)
    return 0


def bed_bp(path: Path) -> int:
    if not path.exists():
        return 0
    total = 0
    for line in path.read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        total += int(parts[2]) - int(parts[1])
    return total


def parse_outlier_bed(path: Path) -> list[tuple[int, int, list[str]]]:
    """Return list of (start, end, [individual_ids]) from a step-3 outlier BED."""
    rows = []
    if not path.exists():
        return rows
    for line in path.read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 4:
            continue
        start, end = int(parts[1]), int(parts[2])
        individuals = [x.strip() for x in parts[3].split(",") if x.strip()]
        rows.append((start, end, individuals))
    return rows


# ---------------------------------------------------------------------------
# Image embedding
# ---------------------------------------------------------------------------

def img_b64(path: Path) -> str | None:
    if not path.exists():
        return None
    data = base64.b64encode(path.read_bytes()).decode("ascii")
    return f"data:image/png;base64,{data}"


# ---------------------------------------------------------------------------
# Stats helpers
# ---------------------------------------------------------------------------

def _mean_sd(values: list[float]) -> tuple[float, float]:
    if not values:
        return 0.0, 0.0
    n = len(values)
    m = sum(values) / n
    if n < 2:
        return m, 0.0
    var = sum((x - m) ** 2 for x in values) / (n - 1)
    return m, math.sqrt(var)


def fmt_bp(bp: float) -> str:
    if bp >= 1_000_000:
        return f"{bp / 1_000_000:.1f} Mb"
    if bp >= 1_000:
        return f"{bp / 1_000:.0f} kb"
    return f"{int(bp)} bp"


def fmt_meansd(values: list[float], fmt=fmt_bp) -> str:
    if not values:
        return "—"
    m, sd = _mean_sd(values)
    if len(values) < 2 or sd < 0.5:
        return fmt(m)
    return f"{fmt(m)} ± {fmt(sd)}"


def fmt_pct_meansd(values: list[float]) -> str:
    if not values:
        return "—"
    mean, sd = _mean_sd(values)
    if len(values) < 2 or sd < 0.05:
        return f"{mean:.1f}%"
    return f"{mean:.1f}% ± {sd:.1f}%"


def percentages_of_length(values: list[float], length: float) -> list[float]:
    if length <= 0:
        return []
    return [100.0 * value / length for value in values]


def pct_color(pct: float) -> str:
    if pct >= 85:
        return "#27ae60"
    if pct >= 70:
        return "#e67e22"
    return "#c0392b"


def bar_html(pct: float) -> str:
    w = max(0, min(120, int(pct * 1.2)))
    color = pct_color(pct)
    return (
        f'<div class="bar-wrap">'
        f'<div class="bar-fill" style="width:{w}px;background:{color}"></div>'
        f"</div>"
    )


# ---------------------------------------------------------------------------
# Data collection
# ---------------------------------------------------------------------------

def index_filtered_ts(paths: list[Path]) -> dict[tuple[str, str], Path]:
    """Index final tree paths by their pipeline <chrom>/<rep> layout."""
    indexed = {}
    for path in paths:
        key = (path.parent.name, path.stem)
        if key in indexed:
            raise ValueError(
                f"Duplicate final tree sequence for chromosome/replicate {key}: "
                f"{indexed[key]} and {path}"
            )
        indexed[key] = path
    return indexed


def retained_bp_from_final_ts(path: Path) -> float:
    """Accessible bp with a non-empty genealogy in one final filtered ARG."""
    ts = load_ts(path)
    metadata = ts.metadata or {}
    accessible = metadata.get("kept_intervals")
    if accessible is None:
        mu = ratemap_from_metadata(metadata)
        accessible = accessible_intervals_from_mu(mu) if mu is not None else None
    return tree_covered_accessible_bp(ts, accessible)


def collect_mu_sources(filtered_ts: list[Path]) -> list[dict]:
    """Report which mutation-rate source step 4 used for each final ARG.

    Reads the ``mu_source`` stamp written by trim_regions_single.py. A "scalar"
    entry means no mutation map was found and a flat rate was substituted, which
    removes the local-rate correction the step-3 outlier test depends on — so it
    is surfaced rather than left silent. ARGs produced before this stamp existed
    report "unrecorded".
    """
    rows = []
    for path in sorted(filtered_ts, key=lambda p: (p.parent.name, p.stem)):
        try:
            metadata = load_ts(path).metadata or {}
        except Exception as exc:  # noqa: BLE001 - reporting must not abort the summary
            rows.append({"chrom": path.parent.name, "rep": path.stem,
                         "kind": "unreadable", "detail": str(exc)})
            continue
        source = metadata.get("mu_source")
        if not source:
            rows.append({"chrom": path.parent.name, "rep": path.stem,
                         "kind": "unrecorded",
                         "detail": "no mu_source stamp (ARG predates this field)"})
            continue
        rows.append({"chrom": path.parent.name, "rep": path.stem,
                     "kind": source.get("kind", "unknown"),
                     "detail": describe_mu_source(source)})
    return rows


def mu_source_section(mu_sources: list[dict]) -> str:
    """Render the mutation-rate provenance block, flagging scalar fallbacks."""
    if not mu_sources:
        return ""
    counts: dict[str, int] = {}
    for row in mu_sources:
        counts[row["kind"]] = counts.get(row["kind"], 0) + 1
    total = len(mu_sources)

    order = ["sibling", "metadata", "scalar", "unrecorded", "unreadable"]
    labels = {
        "sibling": "mutation map file (<code>*.mut_rate.p</code>)",
        "metadata": "ratemap embedded in the input ARG",
        "scalar": "flat scalar <code>mutation_rate</code> fallback",
        "unrecorded": "not recorded",
        "unreadable": "could not be read",
    }
    items = "".join(
        f'<div class="config-item"><span>{labels.get(k, k)}:</span> '
        f'{counts[k]} / {total}</div>'
        for k in order if k in counts
    )

    warn = ""
    fallback = [r for r in mu_sources if r["kind"] == "scalar"]
    if fallback:
        listed = ", ".join(f'{r["chrom"]}/{r["rep"]}' for r in fallback[:12])
        if len(fallback) > 12:
            listed += f", … (+{len(fallback) - 12} more)"
        warn = (
            '<p class="note" style="border-left:4px solid #c33;padding-left:10px">'
            f'<b>{len(fallback)} of {total}</b> filtered ARGs fell back to the flat '
            'scalar <code>mutation_rate</code> because no mutation map was found at '
            'an exact sibling path. A flat rate removes the local mutation-rate '
            'correction that the step-3 outlier test is designed to apply, so those '
            'outlier calls are made against a uniform-rate expectation. '
            f'Affected: {listed}.</p>'
        )
    return (
        "<h2>Mutation-rate source</h2>"
        '<p class="note">Resolution order is: ratemap embedded in the input ARG, '
        'then an exact sibling <code>*.mut_rate.p</code>, then the scalar '
        '<code>mutation_rate</code> config value. Recorded per ARG at step 4.</p>'
        f'<div class="config-grid">{items}</div>{warn}'
    )


def collect_retention(chroms, replicates, out_dir, chrom_lengths, filtered_ts):
    step1_dir = out_dir / "step1_low_rec"
    step2_dir = out_dir / "step2_low_access"
    step3_dir = out_dir / "step3_mutload"
    step4_mask_dir = out_dir / "step4_masks"

    by_num = _fai_num_index(chrom_lengths)
    final_paths = index_filtered_ts(filtered_ts)
    rows = []
    for chrom in chroms:
        seq_len = lookup_chrom_len(chrom, chrom_lengths, by_num)
        s1 = bed_bp(step1_dir / f"{chrom}.low_rec.mask.bed")
        s2 = bed_bp(step2_dir / chrom / f"{chrom}.low_access.bed")

        s3_vals = [
            bed_bp(step3_dir / chrom / f"{rep}.mutation_masked.bed")
            for rep in replicates
        ]
        combined_vals = [
            bed_bp(step4_mask_dir / chrom / f"{rep}.remove_regions.bed")
            for rep in replicates
        ]

        retained_by_rep = {
            rep: retained_bp_from_final_ts(final_paths[(chrom, rep)])
            for rep in replicates
            if (chrom, rep) in final_paths
        }
        retained_vals = list(retained_by_rep.values())
        pct_vals = [100.0 * r / seq_len for r in retained_vals if seq_len > 0]
        mean_pct, _ = _mean_sd(pct_vals)

        rows.append({
            "chrom": chrom,
            "seq_len": seq_len,
            "s1_bp": s1,
            "s1_pct": 100.0 * s1 / seq_len if seq_len else 0,
            "s2_bp": s2,
            "s2_pct": 100.0 * s2 / seq_len if seq_len else 0,
            "s3_vals": s3_vals,
            "combined_vals": combined_vals,
            "retained_by_rep": retained_by_rep,
            "retained_vals": retained_vals,
            "mean_pct": mean_pct,
        })
    return rows


def collect_outliers(chroms, replicates, out_dir):
    """Per-individual: outlier windows and bp removed, aggregated across chroms per rep."""
    step3_dir = out_dir / "step3_mutload"

    # ind -> rep -> (window_count, bp_removed)
    by_ind_rep: dict[str, dict[str, tuple[int, int]]] = defaultdict(lambda: defaultdict(lambda: (0, 0)))
    ind_chroms: dict[str, set] = defaultdict(set)

    for chrom in chroms:
        for rep in replicates:
            rows = parse_outlier_bed(step3_dir / chrom / f"{rep}.outliers.bed")
            for start, end, individuals in rows:
                span = end - start
                for ind in individuals:
                    wc, bp = by_ind_rep[ind][rep]
                    by_ind_rep[ind][rep] = (wc + 1, bp + span)
                    ind_chroms[ind].add(chrom)

    summary = []
    for ind, rep_data in by_ind_rep.items():
        window_counts = [rep_data[rep][0] for rep in replicates]
        bp_vals = [rep_data[rep][1] for rep in replicates]
        mean_w, sd_w = _mean_sd(window_counts)
        mean_bp, _ = _mean_sd(bp_vals)
        summary.append({
            "ind": ind,
            "mean_windows": mean_w,
            "sd_windows": sd_w,
            "mean_bp": mean_bp,
            "chroms": sorted(ind_chroms[ind]),
        })

    summary.sort(key=lambda r: -r["mean_windows"])
    return summary


def weighted_retained_pct(retention, replicates) -> float:
    retained_by_rep: dict[str, float] = defaultdict(float)
    length_by_rep: dict[str, float] = defaultdict(float)
    for row in retention:
        chrom_len = float(row["seq_len"])
        for rep, retained in row["retained_by_rep"].items():
            retained_by_rep[rep] += float(retained)
            length_by_rep[rep] += chrom_len

    pct_vals = [
        100.0 * retained_by_rep[rep] / length_by_rep[rep]
        for rep in replicates
        if length_by_rep[rep]
    ]
    mean_pct, _ = _mean_sd(pct_vals)
    return mean_pct


def totals_by_replicate(retention, replicates, key) -> list[float]:
    """Sum a per-replicate row metric across chromosomes before summarising."""
    totals = []
    for rep_i, rep in enumerate(replicates):
        total = 0.0
        found = False
        for row in retention:
            if key == "retained_vals":
                if rep not in row["retained_by_rep"]:
                    continue
                value = row["retained_by_rep"][rep]
            else:
                values = row[key]
                if rep_i >= len(values):
                    continue
                value = values[rep_i]
            total += float(value)
            found = True
        if found:
            totals.append(total)
    return totals


# ---------------------------------------------------------------------------
# HTML building blocks
# ---------------------------------------------------------------------------

CSS = """
body{font-family:sans-serif;font-size:14px;margin:0;background:#f8f8f8;color:#222}
h1{background:#2c3e50;color:white;margin:0;padding:16px 24px;font-size:20px}
h2{font-size:15px;color:#2c3e50;border-bottom:2px solid #2c3e50;padding-bottom:4px;margin-top:28px}
h3{font-size:13px;color:#555;margin:12px 0 6px 0}
.container{max-width:1200px;margin:0 auto;padding:20px 24px}
.meta{font-size:12px;color:#888;margin-top:4px}
.config-grid{display:grid;grid-template-columns:repeat(3,1fr);gap:8px 24px;margin-top:10px}
.config-item{font-size:12px}.config-item span{font-weight:bold}
table{border-collapse:collapse;width:100%;font-size:13px;margin-top:8px}
th{background:#2c3e50;color:white;padding:6px 10px;text-align:right;font-weight:normal}
th:first-child{text-align:left}
td{padding:5px 10px;border-bottom:1px solid #e0e0e0;text-align:right}
td:first-child{text-align:left;font-weight:bold}
tr:hover{background:#f0f4f8}
.total-row td{font-weight:bold;background:#eef2f7}
.bar-wrap{background:#ddd;border-radius:3px;height:12px;width:120px;display:inline-block;vertical-align:middle}
.bar-fill{height:12px;border-radius:3px}
.plot-grid{display:grid;grid-template-columns:1fr 1fr;gap:16px;margin-top:8px}
.plot-col h3{text-align:center;color:#888;font-size:11px;text-transform:uppercase;letter-spacing:.05em;margin:0 0 8px 0}
.plot-label{font-size:11px;color:#aaa;margin:10px 0 2px 0}
.plot-placeholder{background:#e8e8e8;border:1px solid #ccc;border-radius:4px;height:120px;
  display:flex;align-items:center;justify-content:center;color:#aaa;font-size:12px;font-style:italic}
img.plot{width:100%;border:1px solid #ddd;border-radius:4px;margin-bottom:4px}
.chrom-block{margin-bottom:32px;border-top:1px solid #ddd;padding-top:12px}
.note{font-size:12px;color:#666;margin:0 0 6px 0}
.warn{color:#c0392b;font-weight:bold}
.footer{max-width:1200px;margin:24px auto 0;padding:12px 24px;border-top:1px solid #ddd;
  font-size:11px;color:#999}
.footer code{font-family:monospace;color:#666}
"""

PLOTS = [
    ("mutational-load.png",              "Mutational load by sample"),
    ("diversity-scatter.png",            "Diversity: site vs branch/sim"),
    ("diversity-skyline.png",            "Diversity skyline"),
    ("tajima-d-scatter.png",             "Tajima's D: site vs branch/sim"),
    ("tajima-d-skyline.png",             "Tajima's D skyline"),
    ("frequency-spectrum-unfolded.png",  "SFS — unfolded"),
    ("frequency-spectrum-folded.png",    "SFS — folded"),
]


def plot_img(uri: str | None, label: str) -> str:
    if uri:
        return f'<div class="plot-label">{label}</div><img class="plot" src="{uri}">'
    return (
        f'<div class="plot-label">{label}</div>'
        f'<div class="plot-placeholder">{label} not found</div>'
    )


def validation_section(chroms, out_dir, step6_dir) -> str:
    if not step6_dir.exists():
        return ""
    parts = ['<h2>Validation plots</h2>',
             '<p class="note">Original (pre-pipeline) left · Cleaned (post-pipeline) right. '
             'Posterior means across replicates; shaded bands = 95% CI.</p>']

    for chrom in chroms:
        orig_dir = step6_dir / chrom / "original"
        cln_dir  = step6_dir / chrom / "cleaned"
        if not orig_dir.exists() and not cln_dir.exists():
            continue

        parts.append(f'<div class="chrom-block"><h3>{chrom}</h3><div class="plot-grid">')
        for label, col_dir in [("Original", orig_dir), ("Cleaned", cln_dir)]:
            parts.append(f'<div class="plot-col"><h3>{label}</h3>')
            for fname, plabel in PLOTS:
                parts.append(plot_img(img_b64(col_dir / fname), plabel))
            parts.append("</div>")
        parts.append("</div></div>")

    return "\n".join(parts)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()
    out_dir = args.out_dir
    chroms = args.chroms
    replicates = args.replicates

    chrom_lengths = read_fai(args.fai)
    retention = collect_retention(
        chroms,
        replicates,
        out_dir,
        chrom_lengths,
        args.filtered_ts,
    )
    outliers  = collect_outliers(chroms, replicates, out_dir)

    mu_sources = collect_mu_sources(args.filtered_ts)
    mu_source_html = mu_source_section(mu_sources)

    step6_dir = out_dir / "step6_validation"
    val_html  = validation_section(chroms, out_dir, step6_dir)

    mu_str = args.mutation_rate if args.mutation_rate not in (None, "null") else "—"
    sim_str = "yes" if args.sim_branch.lower() == "true" else "no"
    version = git_version()

    # ---- Config grid -------------------------------------------------------
    def cfg(label, val):
        return f'<div class="config-item"><span>{label}:</span> {val}</div>'

    config_html = "\n".join([
        cfg("chromosomes", len(chroms)),
        cfg("replicates", len(replicates)),
        cfg("rec_fraction", args.rec_fraction if args.rec_fraction is not None else "—"),
        cfg("low_access_window", fmt_bp(args.low_access_window) if args.low_access_window else "—"),
        cfg("low_access_cutoff", fmt_bp(args.low_access_cutoff) if args.low_access_cutoff else "—"),
        cfg("mutload_cutoff", args.mutload_cutoff if args.mutload_cutoff is not None else "—"),
        cfg("validation_mu", mu_str),
        cfg("sim_branch", sim_str),
    ])

    # ---- Retention table ---------------------------------------------------
    ret_rows = []
    total_len = sum(r["seq_len"] for r in retention)
    total_s1  = sum(r["s1_bp"] for r in retention)
    total_s2  = sum(r["s2_bp"] for r in retention)
    all_s3 = totals_by_replicate(retention, replicates, "s3_vals")
    all_comb = totals_by_replicate(retention, replicates, "combined_vals")
    all_retained = totals_by_replicate(retention, replicates, "retained_vals")
    mean_total_pct = weighted_retained_pct(retention, replicates)

    for r in retention:
        pct = r["mean_pct"]
        s3_pct_vals = percentages_of_length(r["s3_vals"], r["seq_len"])
        combined_pct_vals = percentages_of_length(r["combined_vals"], r["seq_len"])
        color = pct_color(pct)
        warn = " <span class='warn'>⚠</span>" if pct < 70 else ""
        ret_rows.append(
            f"<tr>"
            f"<td>{r['chrom']}</td>"
            f"<td>{fmt_bp(r['seq_len'])}</td>"
            f"<td>{fmt_bp(r['s1_bp'])} ({r['s1_pct']:.1f}%)</td>"
            f"<td>{fmt_bp(r['s2_bp'])} ({r['s2_pct']:.1f}%)</td>"
            f"<td>{fmt_meansd(r['s3_vals'])} ({fmt_pct_meansd(s3_pct_vals)})</td>"
            f"<td>{fmt_meansd(r['combined_vals'])} "
            f"({fmt_pct_meansd(combined_pct_vals)})</td>"
            f"<td style='color:{color}'>{fmt_meansd(r['retained_vals'])} "
            f"({pct:.1f}%){warn}</td>"
            f"<td>{bar_html(pct)}</td>"
            f"</tr>"
        )

    total_s1_pct = 100.0 * total_s1 / total_len if total_len else 0
    total_s2_pct = 100.0 * total_s2 / total_len if total_len else 0
    all_s3_pct = percentages_of_length(all_s3, total_len)
    all_comb_pct = percentages_of_length(all_comb, total_len)
    ret_rows.append(
        f"<tr class='total-row'>"
        f"<td>All</td>"
        f"<td>{fmt_bp(total_len)}</td>"
        f"<td>{fmt_bp(total_s1)} ({total_s1_pct:.1f}%)</td>"
        f"<td>{fmt_bp(total_s2)} ({total_s2_pct:.1f}%)</td>"
        f"<td>{fmt_meansd(all_s3)} ({fmt_pct_meansd(all_s3_pct)})</td>"
        f"<td>{fmt_meansd(all_comb)} ({fmt_pct_meansd(all_comb_pct)})</td>"
        f"<td>{fmt_meansd(all_retained)} ({mean_total_pct:.1f}%)</td>"
        f"<td>{bar_html(mean_total_pct)}</td>"
        f"</tr>"
    )

    # ---- Outlier table -----------------------------------------------------
    MAX_OUTLIER_ROWS = 30
    outlier_rows = []
    for o in outliers[:MAX_OUTLIER_ROWS]:
        chrom_str = ", ".join(o["chroms"]) if len(o["chroms"]) <= 4 else \
            f"{', '.join(o['chroms'][:4])} +{len(o['chroms'])-4}"
        sd_str = f" ± {o['sd_windows']:.1f}" if o["sd_windows"] >= 0.5 else ""
        outlier_rows.append(
            f"<tr>"
            f"<td>{o['ind']}</td>"
            f"<td>{o['mean_windows']:.1f}{sd_str}</td>"
            f"<td>{chrom_str}</td>"
            f"<td>{fmt_bp(o['mean_bp'])}</td>"
            f"</tr>"
        )
    if len(outliers) > MAX_OUTLIER_ROWS:
        n_extra = len(outliers) - MAX_OUTLIER_ROWS
        outlier_rows.append(
            f"<tr><td colspan='4' style='color:#aaa;font-style:italic'>"
            f"… {n_extra} additional individuals …</td></tr>"
        )
    if not outlier_rows:
        outlier_rows.append(
            "<tr><td colspan='4' style='color:#aaa'>No outlier windows detected.</td></tr>"
        )

    # ---- Assemble HTML -----------------------------------------------------
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>Pipeline Summary</title>
<style>{CSS}</style>
</head>
<body>
<h1>ARG Pipeline Summary</h1>
<div class="container">

<h2>Run info</h2>
<div class="config-grid">
{config_html}
</div>

<h2>Genome retention by chromosome</h2>
<p class="note">Steps 1–2 masks are replicate-independent.
Step-3, combined-mask, and retained bp are shown as mean ± SD across
{len(replicates)} replicates. Percentages use the full reference length as
their denominator: each chromosome's full .fai length in chromosome rows, and
the sum of those lengths in the All row. Retained bp are measured directly
from the final filtered ARG and exclude both initially inaccessible/empty
regions and pipeline masks.</p>
<table>
<thead><tr>
  <th>Chrom</th><th>Length</th>
  <th>Step 1 masked Mb (% of chromosome)</th>
  <th>Step 2 masked Mb (% of chromosome)</th>
  <th>Step 3 masked Mb (% of chromosome)</th>
  <th>Pipeline removed Mb (% of chromosome; union)</th>
  <th>Final retained Mb (% of reference)</th><th></th>
</tr></thead>
<tbody>
{"".join(ret_rows)}
</tbody>
</table>

<h2>Sample trimming (step 5)</h2>
<p class="note">Individuals with the most outlier windows.
Mean ± SD across replicates; bp removed averaged over all chromosomes and replicates.</p>
<table style="width:auto">
<thead><tr>
  <th style="text-align:left">Individual</th>
  <th>Outlier windows (mean ± SD)</th>
  <th>Chromosomes</th>
  <th>Bp removed (mean)</th>
</tr></thead>
<tbody>
{"".join(outlier_rows)}
</tbody>
</table>

{mu_source_html}

{val_html}

</div>
<div class="footer">Generated by argtest <code>{version}</code></div>
</body>
</html>
"""

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(html)
    print(f"Wrote: {args.out}")


if __name__ == "__main__":
    main()
