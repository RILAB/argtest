import pickle
from pathlib import Path

import msprime


import find_low_access_regions as far
from argtest_common import infer_mu_path


def test_infer_mu_path_uses_parent_chromosome_name(tmp_path):
    chrom_dir = tmp_path / "chr1"
    chrom_dir.mkdir()
    ts_path = chrom_dir / "1.tsz"
    ts_path.touch()
    mu_path = tmp_path / "chr1.mut_rate.p"
    with open(mu_path, "wb") as fh:
        pickle.dump(msprime.RateMap(position=[0, 10], rate=[1.0]), fh)

    inferred = infer_mu_path(ts_path)
    assert inferred == mu_path


def test_infer_mu_path_finds_chrom_map_above_tree_subdir(tmp_path):
    # singer-snakemake layout: <chrom>/trees/<chrom>.<rep>.tsz with the
    # chromosome-wide map at <chrom>/<chrom>.mut_rate.p. Neither the ts stem
    # ("<chrom>.100") nor the parent name ("trees") matches that file.
    chrom = "chr.1.combined.snp.te.sorted"
    trees_dir = tmp_path / chrom / "trees"
    trees_dir.mkdir(parents=True)
    ts_path = trees_dir / f"{chrom}.100.tsz"
    ts_path.touch()
    mu_path = tmp_path / chrom / f"{chrom}.mut_rate.p"
    with open(mu_path, "wb") as fh:
        pickle.dump(msprime.RateMap(position=[0, 10], rate=[1.0]), fh)

    assert infer_mu_path(ts_path) == mu_path


def test_infer_mu_path_prefers_replicate_map_over_chrom_map(tmp_path):
    chrom = "chr1"
    trees_dir = tmp_path / chrom / "trees"
    trees_dir.mkdir(parents=True)
    ts_path = trees_dir / f"{chrom}.7.tsz"
    ts_path.touch()
    rate = msprime.RateMap(position=[0, 10], rate=[1.0])
    for path in (trees_dir / f"{chrom}.7.mut_rate.p", tmp_path / chrom / f"{chrom}.mut_rate.p"):
        with open(path, "wb") as fh:
            pickle.dump(rate, fh)

    assert infer_mu_path(ts_path) == trees_dir / f"{chrom}.7.mut_rate.p"


def test_find_low_access_regions_default_log_goes_under_logs(tmp_path, monkeypatch):
    ts_dir = tmp_path / "trees"
    ts_dir.mkdir()
    ts_path = ts_dir / "base.chr1.1.trees"
    ts_path.touch()
    mu_path = ts_dir / "base.chr1.1.mut_rate.p"
    with open(mu_path, "wb") as fh:
        pickle.dump(msprime.RateMap(position=[0, 3, 10], rate=[1.0, 0.0]), fh)
    out_bed = tmp_path / "low_access.bed"

    monkeypatch.setattr(
        far,
        "parse_args",
        lambda: type(
            "A",
            (),
            {
                "ts": ts_path,
                "chrom": "chr1",
                "window_size": 5.0,
                "cutoff_bp": 4.0,
                "out": out_bed,
                "log": None,
                "mutation_rate": None,
            },
        )(),
    )
    monkeypatch.setattr(
        far,
        "load_ts",
        lambda path: type("TS", (), {"sequence_length": 10.0, "metadata": {}})(),
    )
    monkeypatch.setattr(
        far,
        "accessible_intervals_from_mu",
        lambda mu: __import__("numpy").array([[0.0, 3.0]]),
    )
    far.main()

    assert out_bed.exists()
    log = (tmp_path / "logs" / "low_access.log").read_text()
    assert "mutation_rate_source=sibling\n" in log
    assert f"mutation_rate_detail=map file {mu_path}\n" in log


def test_find_low_access_regions_uses_scalar_mutation_rate_fallback(tmp_path, monkeypatch):
    ts_path = tmp_path / "chr1.tsz"
    ts_path.touch()
    out_bed = tmp_path / "low_access.bed"

    monkeypatch.setattr(
        far,
        "parse_args",
        lambda: type(
            "A",
            (),
            {
                "ts": ts_path,
                "chrom": "chr1",
                "window_size": 5.0,
                "cutoff_bp": 4.0,
                "out": out_bed,
                "log": None,
                "mutation_rate": 1.0,
            },
        )(),
    )
    monkeypatch.setattr(
        far,
        "load_ts",
        lambda path: type("TS", (), {"sequence_length": 10.0, "metadata": {}})(),
    )

    far.main()

    assert out_bed.exists()
    assert out_bed.read_text() == ""
    # The summary reads this line to flag the flat-rate fallback in red.
    log = (tmp_path / "logs" / "low_access.log").read_text()
    assert "mutation_rate_source=scalar\n" in log


def test_find_low_access_regions_nonpositive_rate_flags_whole_sequence(tmp_path, monkeypatch):
    # A non-positive scalar rate means "no accessible sequence", so every window
    # has zero accessible bp and is emitted as low-access.
    ts_path = tmp_path / "chr1.tsz"
    ts_path.touch()
    out_bed = tmp_path / "low_access.bed"

    monkeypatch.setattr(
        far,
        "parse_args",
        lambda: type(
            "A",
            (),
            {
                "ts": ts_path,
                "chrom": "chr1",
                "window_size": 5.0,
                "cutoff_bp": 4.0,
                "out": out_bed,
                "log": None,
                "mutation_rate": 0.0,
            },
        )(),
    )
    monkeypatch.setattr(
        far,
        "load_ts",
        lambda path: type("TS", (), {"sequence_length": 10.0, "metadata": {}})(),
    )

    far.main()

    assert out_bed.read_text().splitlines() == [
        "chr1\t0\t5\t0.000",
        "chr1\t5\t10\t0.000",
    ]


def test_find_low_access_regions_uses_metadata_ratemap(tmp_path, monkeypatch):
    ts_path = tmp_path / "chr1.tsz"
    ts_path.touch()
    out_bed = tmp_path / "low_access.bed"
    metadata = {"mu_position": [0.0, 3.0, 10.0], "mu_rate": [1.0, 0.0]}

    monkeypatch.setattr(
        far,
        "parse_args",
        lambda: type(
            "A",
            (),
            {
                "ts": ts_path,
                "chrom": "chr1",
                "window_size": 5.0,
                "cutoff_bp": 4.0,
                "out": out_bed,
                "log": None,
                "mutation_rate": None,
            },
        )(),
    )
    monkeypatch.setattr(
        far,
        "load_ts",
        lambda path: type("TS", (), {"sequence_length": 10.0, "metadata": metadata})(),
    )

    far.main()

    assert out_bed.read_text().splitlines() == [
        "chr1\t0\t5\t3.000",
        "chr1\t5\t10\t0.000",
    ]
