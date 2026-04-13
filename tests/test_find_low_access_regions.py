import pickle
from pathlib import Path

import sys
from types import SimpleNamespace

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import find_low_access_regions as far


def test_infer_mu_path_uses_parent_chromosome_name(tmp_path):
    chrom_dir = tmp_path / "chr1"
    chrom_dir.mkdir()
    ts_path = chrom_dir / "1.tsz"
    ts_path.touch()
    mu_path = tmp_path / "chr1.mut_rate.p"
    with open(mu_path, "wb") as fh:
        pickle.dump(SimpleNamespace(position=[0, 10], rate=[1.0]), fh)

    inferred = far.infer_mu_path(ts_path)
    assert inferred == mu_path


def test_find_low_access_regions_default_log_goes_under_logs(tmp_path, monkeypatch):
    ts_dir = tmp_path / "trees"
    ts_dir.mkdir()
    ts_path = ts_dir / "base.chr1.1.trees"
    ts_path.touch()
    mu_path = ts_dir / "base.chr1.mut_rate.p"
    with open(mu_path, "wb") as fh:
        pickle.dump(SimpleNamespace(position=[0, 3, 10], rate=[1.0, 0.0]), fh)
    out_bed = tmp_path / "low_access.bed"

    monkeypatch.setattr(
        far,
        "parse_args",
        lambda: type(
            "A",
            (),
            {
                "ts": ts_path,
                "window_size": 5.0,
                "cutoff_bp": 4.0,
                "out": out_bed,
                "log": None,
            },
        )(),
    )
    monkeypatch.setattr(
        far,
        "load_ts",
        lambda path: type("TS", (), {"sequence_length": 10.0})(),
    )
    monkeypatch.setattr(
        far,
        "accessible_intervals_from_mu",
        lambda mu: __import__("numpy").array([[0.0, 3.0]]),
    )
    far.main()

    assert out_bed.exists()
    assert (tmp_path / "logs" / "low_access.log").exists()
