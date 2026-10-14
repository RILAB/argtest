from pathlib import Path
from types import SimpleNamespace

import numpy as np


import pipeline_summary as ps


def test_collect_outliers_counts_zero_replicates(tmp_path):
    outlier_dir = tmp_path / "step3_mutload" / "chr1"
    outlier_dir.mkdir(parents=True)
    (outlier_dir / "1.outliers.bed").write_text("chr1\t0\t10\tA\n")

    summary = ps.collect_outliers(["chr1"], ["1", "2"], tmp_path)

    assert len(summary) == 1
    assert summary[0]["ind"] == "A"
    assert summary[0]["mean_windows"] == 0.5
    assert summary[0]["mean_bp"] == 5.0


def test_weighted_retained_pct_uses_chrom_lengths():
    retention = [
        {"seq_len": 100, "retained_by_rep": {"1": 100}},
        {"seq_len": 900, "retained_by_rep": {"1": 0}},
    ]

    assert ps.weighted_retained_pct(retention, ["1"]) == 10.0


def test_mask_percentages_use_chromosome_length():
    percentages = ps.percentages_of_length([100, 200], 1000)

    assert percentages == [10.0, 20.0]
    assert ps.fmt_pct_meansd(percentages) == "15.0% ± 7.1%"


def test_retained_bp_from_final_ts_intersects_input_accessibility_and_tree_coverage(
    tmp_path, monkeypatch
):
    class FakeTree:
        def __init__(self, left, right, num_edges):
            self.interval = SimpleNamespace(left=left, right=right)
            self.num_edges = num_edges

    ts = SimpleNamespace(
        metadata={"mu_position": [0, 10, 20, 30, 40], "mu_rate": [1, 0, 1, 0]},
        trees=lambda: iter([
            FakeTree(0, 15, 1),
            FakeTree(15, 25, 0),
            FakeTree(25, 40, 1),
        ]),
    )
    mu = SimpleNamespace(
        position=np.array([0, 10, 20, 30, 40]),
        rate=np.array([1, 0, 1, 0]),
    )
    monkeypatch.setattr(ps, "load_ts", lambda path: ts)
    monkeypatch.setattr(ps, "ratemap_from_metadata", lambda metadata: mu)

    # Accessible [0,10) contributes 10 bp. Accessible [20,30) overlaps an
    # empty final tree on [20,25), so only [25,30) contributes another 5 bp.
    assert ps.retained_bp_from_final_ts(tmp_path / "1.tsz") == 15.0


def test_retained_bp_prefers_kept_intervals_over_mutation_map(tmp_path, monkeypatch):
    class FakeTree:
        def __init__(self, left, right):
            self.interval = SimpleNamespace(left=left, right=right)
            self.num_edges = 1

    ts = SimpleNamespace(
        metadata={
            "kept_intervals": [[5, 15]],
            "mu_position": [0, 20],
            "mu_rate": [1],
        },
        trees=lambda: iter([FakeTree(0, 20)]),
    )
    monkeypatch.setattr(ps, "load_ts", lambda path: ts)

    assert ps.retained_bp_from_final_ts(tmp_path / "1.tsz") == 10.0


def test_retained_bp_uses_tree_coverage_without_accessibility_metadata(
    tmp_path, monkeypatch
):
    class FakeTree:
        def __init__(self, left, right, num_edges):
            self.interval = SimpleNamespace(left=left, right=right)
            self.num_edges = num_edges

    ts = SimpleNamespace(
        metadata={},
        trees=lambda: iter([FakeTree(0, 8, 1), FakeTree(8, 20, 0)]),
    )
    monkeypatch.setattr(ps, "load_ts", lambda path: ts)

    assert ps.retained_bp_from_final_ts(tmp_path / "1.tsz") == 8.0


def test_all_row_totals_sum_chromosomes_within_each_replicate():
    retention = [
        {
            "combined_vals": [10, 20],
            "retained_by_rep": {"1": 80, "2": 70},
            "retained_vals": [80, 70],
        },
        {
            "combined_vals": [30, 40],
            "retained_by_rep": {"1": 60, "2": 50},
            "retained_vals": [60, 50],
        },
    ]

    assert ps.totals_by_replicate(retention, ["1", "2"], "combined_vals") == [40, 60]
    assert ps.totals_by_replicate(retention, ["1", "2"], "retained_vals") == [140, 120]


# --------------------------------------------------------------------------- #
# Genome-wide section coverage labelling
# --------------------------------------------------------------------------- #

def _write_pooled(step6_dir, variant, chroms):
    d = step6_dir / "genomewide" / variant
    d.mkdir(parents=True, exist_ok=True)
    header = "chrom\tdataset\twindow_start\twindow_end\trec_rate_cm_per_mb\t"
    header += "obs_pi\texp_pi\tobs_tajimas_d\texp_tajimas_d"
    rows = [f"{c}\tprimary\t0.0\t100.0\t1.0\t0.1\t0.2\t-1.0\t-2.0" for c in chroms]
    (d / "genomewide-windows.tsv").write_text("\n".join([header] + rows) + "\n")


def test_pooled_chroms_are_read_from_the_data(tmp_path):
    _write_pooled(tmp_path, "cleaned", ["chr1", "chr2", "chr1"])
    assert ps.pooled_chroms(tmp_path / "genomewide") == {"chr1", "chr2"}


def test_all_chromosomes_pooled_is_called_genome_wide(tmp_path):
    _write_pooled(tmp_path, "cleaned", ["chr1", "chr2"])
    _write_pooled(tmp_path, "original", ["chr1", "chr2"])

    html = ps.genomewide_section(["chr1", "chr2"], tmp_path)

    assert "Genome-wide expected vs observed" in html
    assert "partial genome" not in html.lower()


def test_partial_pool_is_never_captioned_as_genome_wide(tmp_path):
    """The shipped default runs step 6 on one chromosome; the report must say so."""
    _write_pooled(tmp_path, "cleaned", ["chr1"])
    _write_pooled(tmp_path, "original", ["chr1"])

    html = ps.genomewide_section(["chr1", "chr2", "chr3"], tmp_path)

    assert "All windows from all" not in html
    assert "partial genome" in html.lower()
    assert "1 of 3 chromosomes" in html
    # The missing chromosomes are named, and the warning is visually flagged.
    assert "chr2" in html and "chr3" in html
    assert "validation_first_chrom_only" in html
    assert "#c33" in html


def test_absent_genomewide_dir_yields_no_section(tmp_path):
    assert ps.genomewide_section(["chr1"], tmp_path) == ""


def test_dropped_mask_intervals_flow_from_step4_log_to_summary(tmp_path, monkeypatch):
    import combine_remove_masks as crm

    bed = tmp_path / "mask.bed"
    bed.write_text("chr1\t5\t5\nchr1\t0\t10\nchr1\t30\t20\n")
    log = tmp_path / "logs" / "step4_masks" / "chr1" / "1.log"
    monkeypatch.setattr(
        crm,
        "parse_args",
        lambda: SimpleNamespace(
            chrom="chr1", out=tmp_path / "out.bed", inputs=[bed], log=log
        ),
    )
    crm.main()

    dropped = ps.collect_dropped_mask_intervals(["chr1"], ["1", "2"], tmp_path)

    assert [(r["chrom"], r["rep"], r["n"]) for r in dropped] == [("chr1", "1", 2)]
    html = ps.dropped_intervals_section(dropped, tmp_path)
    assert "2 mask intervals dropped" in html
    assert "chr1/1 (2)" in html


def test_no_dropped_mask_intervals_yields_no_section(tmp_path):
    log_dir = tmp_path / "logs" / "step4_masks" / "chr1"
    log_dir.mkdir(parents=True)
    (log_dir / "1.log").write_text(
        "# combine_remove_masks summary\nmerged_intervals=3\ndropped_intervals=0\n"
    )
    # Logs written before dropped_intervals existed lack the field entirely.
    (log_dir / "2.log").write_text("# combine_remove_masks summary\nmerged_intervals=3\n")

    dropped = ps.collect_dropped_mask_intervals(["chr1", "chr2"], ["1", "2"], tmp_path)

    assert dropped == []
    assert ps.dropped_intervals_section(dropped, tmp_path) == ""


# --------------------------------------------------------------------------- #
# Run settings and input provenance
# --------------------------------------------------------------------------- #

def _run_config(**step3):
    return {
        "inputs": {"root_dir": "/data", "tree_subdir": "trees", "tree_pattern": "*.tsz",
                   "burnin": 75, "replicates_found": 150, "replicates_used": ["75", "149"],
                   "missing_pairs": [], "hapmap": "/data/map.tsv", "fai": "/data/g.fai"},
        "step1": {"rec_fraction": 0.005},
        "step2": {"low_access_window_size": 50000.0, "low_access_cutoff_bp": 2500.0},
        "step3": {"mutload_cutoff": 0.5, "mutload_fraction": 0.25,
                  "mutation_rate": 3.3e-8, **step3},
        "validation": {"run_validation": True, "run_validation_requested": True},
    }


def test_settings_section_shows_every_setting_including_unknown_keys():
    cfg = _run_config(brand_new_knob=7)
    cfg["future_step"] = {"some_setting": "x"}

    html, flags = ps.settings_section(cfg, ["chr1"], ["75"], {})

    assert "mutload_fraction:</span> 0.25" in html
    assert "rec_fraction:</span> 0.005" in html
    assert "replicates_used:</span> 75–149" in html
    # Sizes are exact, not rounded to whole kb.
    assert "low_access_cutoff_bp:</span> 2,500 bp" in html
    assert "low_access_window_size:</span> 50,000 bp" in html
    # Keys this script does not list explicitly still render.
    assert "brand_new_knob:</span> 7" in html
    assert "some_setting:</span> x" in html
    # Paths belong to the inputs table, not the settings grid.
    assert "root_dir:" not in html
    assert flags == []


def test_missing_run_config_is_flagged():
    html, flags = ps.settings_section(None, ["chr1"], ["1"], {})
    assert "run_config.json not provided" in html
    assert flags


def _write_logs(out_dir, kind, detail):
    (out_dir / "logs" / "step2").mkdir(parents=True)
    (out_dir / "logs" / "step2" / "chr1.log").write_text(
        f"mutation_rate_source={kind}\nmutation_rate_detail={detail}\n")
    (out_dir / "logs" / "step3" / "chr1").mkdir(parents=True)
    (out_dir / "logs" / "step3" / "chr1" / "75.log").write_text(
        f"cutoff=0.5\nmutation_rate_source={kind}\nmutation_rate_detail={detail}\n")
    (out_dir / "logs" / "step1").mkdir(parents=True)
    (out_dir / "logs" / "step1" / "chr1.log").write_text(
        "hapmap=/data/map.tsv\nfai=/data/g.fai\nchrom_match=chr1\thapmap=1\tfai=chr1\n")


def _provenance(out_dir, cfg=None):
    cfg = cfg or _run_config()
    mu = ps.collect_step_mu_sources(["chr1"], ["75"], out_dir, {})
    html, flags = ps.provenance_section(
        cfg, ps.collect_step1_refs(["chr1"], out_dir), mu, {"chr1": "exact"},
        Path("/data/g.fai"))
    settings_html, settings_flags = ps.settings_section(cfg, ["chr1"], ["75"], mu)
    return html, flags, settings_html


def test_flat_scalar_mutation_rate_is_flagged_red(tmp_path):
    _write_logs(tmp_path, "scalar", "FLAT scalar 3.3e-08 (no local rate correction)")

    html, flags, settings_html = _provenance(tmp_path)

    assert html.count('class="flagged"') == 2  # steps 2 and 3
    assert "FLAT SCALAR 1/1" in html
    assert any(f.startswith("Step 3") and "flat scalar" in f for f in flags)
    assert "used 2× — see Inputs" in settings_html
    banner = ps.flags_banner(flags)
    assert "This run used fallbacks" in banner and "Step 2" in banner


def test_map_files_are_listed_and_aliases_noted(tmp_path):
    _write_logs(tmp_path, "sibling", "map file /data/chr1/chr1.mut_rate.p")

    html, flags, _ = _provenance(tmp_path)

    assert 'class="flagged"' not in html
    assert flags == []
    assert "<code>/data/chr1/chr1.mut_rate.p</code>" in html
    assert "<code>/data/map.tsv</code>" in html
    assert "names matched via alias: chr1 → 1" in html
    assert "<td>step 1 · step 6b</td>" in html
    assert "<td>step 1 · step 6b · summary lengths</td>" in html


def test_reference_provenance_omits_step6b_when_validation_did_not_run(tmp_path):
    _write_logs(tmp_path, "sibling", "map file /data/chr1/chr1.mut_rate.p")
    cfg = _run_config()
    cfg["validation"]["run_validation"] = False

    html, _, _ = _provenance(tmp_path, cfg)

    assert "<td>step 1</td>" in html
    assert "<td>step 1 · summary lengths</td>" in html
    assert "step 6b" not in html


def test_logs_without_mu_source_are_unrecorded_not_flagged(tmp_path):
    (tmp_path / "logs" / "step2").mkdir(parents=True)
    (tmp_path / "logs" / "step2" / "chr1.log").write_text("Flagged 1 windows\n")

    html, flags, _ = _provenance(tmp_path)

    assert "not recorded 1/1 (log predates this field)" in html
    assert flags == []


def test_missing_pairs_and_fuzzy_chrom_lengths_are_flagged(tmp_path):
    cfg = _run_config()
    cfg["inputs"]["missing_pairs"] = [["chr2", "80"]]
    mu = ps.collect_step_mu_sources(["chr1"], ["75"], tmp_path, {})

    html, flags = ps.provenance_section(
        cfg, ps.collect_step1_refs(["chr1"], tmp_path), mu,
        {"chr1": "suffix"}, Path("/data/g.fai"))

    assert "1 chromosome/replicate pairs missing" in html
    assert "matched only by trailing number" in html
    assert len(flags) == 2


def test_chrom_len_source_distinguishes_exact_suffix_and_missing():
    lengths = {"chr1": 10, "scaffold_2": 20}
    by_num = ps._fai_num_index(lengths)
    assert ps.chrom_len_source("chr1", lengths, by_num) == "exact"
    assert ps.chrom_len_source("maize.2", lengths, by_num) == "suffix"
    assert ps.chrom_len_source("chrX", lengths, by_num) == "missing"


def test_vcf_reps_show_a_range_only_for_the_full_replicate_set():
    all_reps = [str(r) for r in range(50, 150)]

    def vcf_line(vcf_reps):
        cfg = _run_config()
        cfg["outputs"] = {"emit_vcf": True, "vcf_reps": vcf_reps}
        html, _ = ps.settings_section(cfg, ["chr1"], all_reps, {})
        return html.split("vcf_reps:</span> ", 1)[1].split("</div>", 1)[0]

    assert vcf_line(all_reps) == "all 100 post-burn-in (50–149)"
    # A subset is never written as a range, which would imply it is contiguous.
    assert vcf_line(["50", "60", "70", "80", "90", "100", "149"]) == "7: 50, 60, 70, …, 149"
    assert vcf_line(["50", "149"]) == "50, 149"
