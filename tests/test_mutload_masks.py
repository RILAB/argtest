from pathlib import Path

import tskit


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
        "name_substring_to_remove": "_anchorwave",
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

    # Per-window expected = 4 for both individuals; observed (1, 0) both fall
    # below (1-0.5)*4 = 2 in the single window, so both flag and appear as
    # comma-joined entries on one BED row.
    import numpy as np

    def _fake(ts, windows, names, mutation_rate=None, seed=None):
        unique = list(dict.fromkeys(names))
        return np.full((len(windows) - 1, len(unique)), 4.0)

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
    lines = [ln for ln in outlier.read_text().splitlines() if ln.strip()]
    # One BED row per outlier window; here both individuals share the row.
    assert len(lines) == 1
    parts = lines[0].split("\t")
    # Simple TS has no individual table, so sample names come back as "node<u>".
    assert parts[3].split(",") == ["node0", "node1"]
    assert masked.exists()
    assert masked.read_text().strip() == ""
    assert (tmp_path / "logs" / "chr1.1.mutload.log").exists()


def test_mutload_masks_outlier_bed_columns(tmp_path, monkeypatch):
    # BED schema: chrom, window-start, window-end, names, observed, expected.
    import numpy as np
    ts = make_simple_ts()
    ts_path = tmp_path / "1.trees"
    ts.dump(ts_path)
    outlier = tmp_path / "outliers.bed"
    masked = tmp_path / "masked.bed"

    # Per-window expected = 4 for both individuals; observed (1, 0) both fall
    # below 0.5*4 so both flag in the single window.
    def _fake(ts, windows, names, mutation_rate=None, seed=None):
        unique = list(dict.fromkeys(names))
        return np.full((len(windows) - 1, len(unique)), 4.0)

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
    # Expected col carries the per-individual sim expectation in this window.
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


def _run_one_window_half_outlier(tmp_path, monkeypatch, fraction):
    """Run main() on a single window where exactly 1 of 2 individuals is an outlier.

    Observed load is (1, 0); expected is faked as (1, 4) with cutoff 0.5, so
    node0 sits inside [0.5, 1.5] and node1 falls below 2.0. The window's outlier
    fraction is therefore exactly 0.5, which is what makes it a boundary case
    for --fraction.
    """
    import numpy as np

    ts = make_simple_ts()
    ts_path = tmp_path / "1.trees"
    ts.dump(ts_path)
    outlier = tmp_path / "outliers.bed"
    masked = tmp_path / "masked.bed"

    def _fake(ts, windows, names, mutation_rate=None, seed=None):
        unique = list(dict.fromkeys(names))
        expected = np.empty((len(windows) - 1, len(unique)))
        expected[:, 0] = 1.0
        expected[:, 1] = 4.0
        return expected

    monkeypatch.setattr(mm, "simulate_expected_load", _fake)
    monkeypatch.setattr(
        mm,
        "parse_args",
        lambda: _make_args(
            ts=ts_path,
            outlier_bed=outlier,
            masked_bed=masked,
            cutoff=0.5,
            fraction=fraction,
        ),
    )
    mm.main()
    return outlier.read_text(), masked.read_text()


def test_fraction_masks_window_when_outlier_fraction_is_above_threshold(
    tmp_path, monkeypatch
):
    outlier_text, masked_text = _run_one_window_half_outlier(
        tmp_path, monkeypatch, fraction=0.4
    )

    # 0.5 > 0.4, so the whole window is masked: chrom, start, end, fraction,
    # outlier count, individual count.
    assert masked_text.strip().split("\t") == ["chr1", "0", "10", "0.500", "1", "2"]
    # A masked window is dropped from the per-individual outlier BED.
    assert outlier_text.strip() == ""


def test_fraction_keeps_window_when_outlier_fraction_equals_threshold(
    tmp_path, monkeypatch
):
    outlier_text, masked_text = _run_one_window_half_outlier(
        tmp_path, monkeypatch, fraction=0.5
    )

    # The comparison is strict (fraction > threshold), so equality does NOT mask.
    assert masked_text.strip() == ""
    assert outlier_text.strip().split("\t")[3] == "node1"


def test_fraction_keeps_window_when_outlier_fraction_is_below_threshold(
    tmp_path, monkeypatch
):
    outlier_text, masked_text = _run_one_window_half_outlier(
        tmp_path, monkeypatch, fraction=0.6
    )

    assert masked_text.strip() == ""
    assert outlier_text.strip().split("\t")[3] == "node1"


def test_bed_bounds_round_outward_without_collapsing_short_window():
    assert mm.bed_bounds(1.2, 1.8) == (1, 2)
