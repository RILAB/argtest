#!/usr/bin/env python
"""Generate a self-contained HTML pipeline summary from Snakemake outputs."""
from __future__ import annotations

import argparse
import base64
import math
import re
from collections import defaultdict
from pathlib import Path


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(description="Generate pipeline summary HTML.")
    p.add_argument("--out-dir", required=True, type=Path)
    p.add_argument("--fai", required=True, type=Path)
    p.add_argument("--chroms", required=True, nargs="+")
    p.add_argument("--replicates", required=True, nargs="+")
    p.add_argument("--out", required=True, type=Path)
    p.add_argument("--rec-fraction", type=float, default=None)
    p.add_argument("--low-access-window", type=int, default=None)
    p.add_argument("--low-access-cutoff", type=int, default=None)
    p.add_argument("--mutload-cutoff", type=float, default=None)
    p.add_argument("--mutation-rate", default=None)
    p.add_argument("--sim-branch", default="false")
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

def collect_retention(chroms, replicates, out_dir, chrom_lengths):
    step1_dir = out_dir / "step1_low_rec"
    step2_dir = out_dir / "step2_low_access"
    step3_dir = out_dir / "step3_mutload"
    step4_mask_dir = out_dir / "step4_masks"

    rows = []
    for chrom in chroms:
        seq_len = chrom_lengths.get(chrom, 0)
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

        retained_vals = [seq_len - c for c in combined_vals if seq_len > 0]
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
            "retained_vals": retained_vals,
            "pct_vals": pct_vals,
            "mean_pct": mean_pct,
        })
    return rows


def collect_outliers(chroms, replicates, out_dir):
    """Per-individual: windows flagged and bp removed, aggregated across chroms per rep."""
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
        window_counts = [wc for wc, _ in rep_data.values()]
        bp_vals = [bp for _, bp in rep_data.values()]
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
    retention = collect_retention(chroms, replicates, out_dir, chrom_lengths)
    outliers  = collect_outliers(chroms, replicates, out_dir)

    step6_dir = out_dir / "step6_validation"
    val_html  = validation_section(chroms, out_dir, step6_dir)

    mu_str = args.mutation_rate if args.mutation_rate not in (None, "null") else "—"
    sim_str = "yes" if args.sim_branch.lower() == "true" else "no"

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
    all_s3    = [v for r in retention for v in r["s3_vals"]]
    all_comb  = [v for r in retention for v in r["combined_vals"]]
    all_ret   = [v for r in retention for v in r["retained_vals"]]
    all_pct   = [v for r in retention for v in r["pct_vals"]]
    mean_total_pct, _ = _mean_sd(all_pct)

    for r in retention:
        pct = r["mean_pct"]
        color = pct_color(pct)
        warn = " <span class='warn'>⚠</span>" if pct < 70 else ""
        ret_rows.append(
            f"<tr>"
            f"<td>{r['chrom']}</td>"
            f"<td>{fmt_bp(r['seq_len'])}</td>"
            f"<td>{fmt_bp(r['s1_bp'])} ({r['s1_pct']:.1f}%)</td>"
            f"<td>{fmt_bp(r['s2_bp'])} ({r['s2_pct']:.1f}%)</td>"
            f"<td>{fmt_meansd(r['s3_vals'])}</td>"
            f"<td>{fmt_meansd(r['combined_vals'])}</td>"
            f"<td style='color:{color}'>{pct:.1f}%{warn}</td>"
            f"<td>{bar_html(pct)}</td>"
            f"</tr>"
        )

    total_s1_pct = 100.0 * total_s1 / total_len if total_len else 0
    total_s2_pct = 100.0 * total_s2 / total_len if total_len else 0
    ret_rows.append(
        f"<tr class='total-row'>"
        f"<td>All</td>"
        f"<td>{fmt_bp(total_len)}</td>"
        f"<td>{fmt_bp(total_s1)} ({total_s1_pct:.1f}%)</td>"
        f"<td>{fmt_bp(total_s2)} ({total_s2_pct:.1f}%)</td>"
        f"<td>{fmt_meansd(all_s3)}</td>"
        f"<td>{fmt_meansd(all_comb)}</td>"
        f"<td>{mean_total_pct:.1f}%</td>"
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
Step-3 and combined-mask bp shown as mean ± SD across {len(replicates)} replicates.</p>
<table>
<thead><tr>
  <th>Chrom</th><th>Length</th>
  <th>Step 1 masked (low-rec)</th>
  <th>Step 2 masked (low-access)</th>
  <th>Step 3 masked (mutload)</th>
  <th>Total removed (combined)</th>
  <th>Retained</th><th></th>
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

{val_html}

</div>
</body>
</html>
"""

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(html)
    print(f"Wrote: {args.out}")


if __name__ == "__main__":
    main()
