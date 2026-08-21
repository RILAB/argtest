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
