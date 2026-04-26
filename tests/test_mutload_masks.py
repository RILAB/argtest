from pathlib import Path

import sys
import tskit

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import mutload_masks as mm


def make_simple_ts(length=10):
    tables = tskit.TableCollection(sequence_length=length)
    pop = tables.populations.add_row()
    n0 = tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0, population=pop)
    n1 = tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0, population=pop)
    anc = tables.nodes.add_row(time=1, population=pop)
    tables.edges.add_row(left=0, right=length, parent=anc, child=n0)
    tables.edges.add_row(left=0, right=length, parent=anc, child=n1)
    site = tables.sites.add_row(position=1, ancestral_state="0")
    tables.mutations.add_row(site=site, node=n0, derived_state="1")
    tables.sort()
    return tables.tree_sequence()


def _make_args(**overrides):
    base = {
        "ts": None,
        "chrom": "chr1",
        "window_size": 10.0,
        "snp_window": None,
        "cutoff": 0.5,
        "fraction": None,
        "mutation_rate": 1.0,
        "random_seed": 1,
        "suffix_to_strip": "_anchorwave",
        "outlier_bed": None,
        "masked_bed": None,
        "log": None,
    }
    base.update(overrides)
    return type("A", (), base)()


def test_mutload_masks_writes_outlier_and_empty_masked(tmp_path, monkeypatch):
    ts = make_simple_ts()
    ts_path = tmp_path / "1.trees"
    ts.dump(ts_path)
    outlier = tmp_path / "outliers.bed"
    masked = tmp_path / "masked.bed"

    # Force a deterministic non-zero expected so the single observed mutation
    # falls below 0.5 * expected and gets flagged.
    import numpy as np

    def _fake(ts, windows, names, mu=None, seed=None):
        unique = list(dict.fromkeys(names))
        return np.full((len(windows) - 1, len(unique)), 4.0), unique

    monkeypatch.setattr(mm, "simulate_expected_load", _fake)
    monkeypatch.setattr(
        mm,
        "parse_args",
        lambda: _make_args(
            ts=ts_path,
            outlier_bed=outlier,
            masked_bed=masked,
        ),
    )
    mm.main()

    assert outlier.exists()
    assert outlier.read_text().strip() != ""
    assert masked.exists()
    assert masked.read_text().strip() == ""
    assert (tmp_path / "logs" / "chr1.1.mutload.log").exists()


def test_mutload_masks_outlier_bed_columns(tmp_path, monkeypatch):
    # Verify the BED schema: chrom, start, end, names, observed, expected.
    import numpy as np
    ts = make_simple_ts()
    ts_path = tmp_path / "1.trees"
    ts.dump(ts_path)
    outlier = tmp_path / "outliers.bed"
    masked = tmp_path / "masked.bed"

    # Both individuals get expected=4 in the only window, so observed=1 and
    # observed=0 are both flagged at cutoff=0.5.
    def _fake(ts, windows, names, mu=None, seed=None):
        unique = list(dict.fromkeys(names))
        return np.full((len(windows) - 1, len(unique)), 4.0), unique

    monkeypatch.setattr(mm, "simulate_expected_load", _fake)
    monkeypatch.setattr(
        mm, "parse_args",
        lambda: _make_args(ts=ts_path, outlier_bed=outlier, masked_bed=masked, cutoff=0.5),
    )
    mm.main()

    parts = outlier.read_text().strip().split("\t")
    assert parts[0] == "chr1"
    assert parts[1] == "0"
    assert parts[2] == "10"
    names = parts[3].split(",")
    observed = parts[4].split(",")
    expected = parts[5].split(",")
    assert len(names) == len(observed) == len(expected)
    # Expected col holds the per-individual sim expectation, not the median.
    assert all(float(x) == 4.0 for x in expected)


def test_mutload_masks_seed_determinism(tmp_path, monkeypatch):
    ts = make_simple_ts()
    ts_path = tmp_path / "1.trees"
    ts.dump(ts_path)
    out1 = tmp_path / "o1.bed"
    masked1 = tmp_path / "m1.bed"
    out2 = tmp_path / "o2.bed"
    masked2 = tmp_path / "m2.bed"

    monkeypatch.setattr(
        mm, "parse_args",
        lambda: _make_args(
            ts=ts_path, outlier_bed=out1, masked_bed=masked1,
            mutation_rate=1.0, random_seed=42,
        ),
    )
    mm.main()
    monkeypatch.setattr(
        mm, "parse_args",
        lambda: _make_args(
            ts=ts_path, outlier_bed=out2, masked_bed=masked2,
            mutation_rate=1.0, random_seed=42,
        ),
    )
    mm.main()
    assert out1.read_text() == out2.read_text()
    assert masked1.read_text() == masked2.read_text()


def test_mutload_masks_requires_rate_when_no_ratemap(tmp_path, monkeypatch):
    import pytest
    ts = make_simple_ts()
    ts_path = tmp_path / "1.trees"
    ts.dump(ts_path)
    monkeypatch.setattr(
        mm, "parse_args",
        lambda: _make_args(
            ts=ts_path,
            outlier_bed=tmp_path / "o.bed",
            masked_bed=tmp_path / "m.bed",
            mutation_rate=None,
        ),
    )
    with pytest.raises(FileNotFoundError):
        mm.main()
