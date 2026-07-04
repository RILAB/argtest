from pathlib import Path

import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

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
        {"seq_len": 100, "retained_vals": [100]},
        {"seq_len": 900, "retained_vals": [0]},
    ]

    assert ps.weighted_retained_pct(retention, ["1"]) == 10.0
