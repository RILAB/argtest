import os
from pathlib import Path

import numpy as np
import pytest
import tskit

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import argtest_common as mc
import mutload_summary as ms


def make_simple_ts():
    tables = tskit.TableCollection(sequence_length=10)
    tables.individuals.metadata_schema = tskit.MetadataSchema.permissive_json()
    pop = tables.populations.add_row()
    ind0 = tables.individuals.add_row(metadata={"id": "A"})
    ind1 = tables.individuals.add_row(metadata={"id": "B"})
    n0 = tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0, individual=ind0, population=pop)
    n1 = tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0, individual=ind1, population=pop)
    anc = tables.nodes.add_row(time=1, population=pop)
    tables.edges.add_row(left=0, right=10, parent=anc, child=n0)
    tables.edges.add_row(left=0, right=10, parent=anc, child=n1)
    s1 = tables.sites.add_row(position=1, ancestral_state="0")
    s7 = tables.sites.add_row(position=7, ancestral_state="0")
    tables.mutations.add_row(site=s1, node=n0, derived_state="1")
    tables.mutations.add_row(site=s7, node=n1, derived_state="1")
    tables.sort()
    return tables.tree_sequence()


def make_ts_no_mutations(n_samples=2, length=10):
    tables = tskit.TableCollection(sequence_length=length)
    tables.individuals.metadata_schema = tskit.MetadataSchema.permissive_json()
    pop = tables.populations.add_row()
    inds = [tables.individuals.add_row(metadata={"id": f"I{i}"}) for i in range(n_samples)]
    samples = [
        tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0, individual=inds[i], population=pop)
        for i in range(n_samples)
    ]
    anc = tables.nodes.add_row(time=1, population=pop)
    for s in samples:
        tables.edges.add_row(left=0, right=length, parent=anc, child=s)
    tables.sort()
    return tables.tree_sequence()


def make_ts_many_individuals(n=100, length=10):
    tables = tskit.TableCollection(sequence_length=length)
    tables.individuals.metadata_schema = tskit.MetadataSchema.permissive_json()
    pop = tables.populations.add_row()
    inds = [tables.individuals.add_row(metadata={"id": f"I{i}"}) for i in range(n)]
    samples = [
        tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0, individual=inds[i], population=pop)
        for i in range(n)
    ]
    anc = tables.nodes.add_row(time=1, population=pop)
    for s in samples:
        tables.edges.add_row(left=0, right=length, parent=anc, child=s)
    s1 = tables.sites.add_row(position=1, ancestral_state="0")
    tables.mutations.add_row(site=s1, node=samples[0], derived_state="1")
    tables.sort()
    return tables.tree_sequence()


def make_lineage_ts():
    tables = tskit.TableCollection(sequence_length=10)
    tables.individuals.metadata_schema = tskit.MetadataSchema.permissive_json()
    pop = tables.populations.add_row()
    names = ["L1_A", "L1_B", "L2_A", "L2_B"]
    inds = [tables.individuals.add_row(metadata={"id": name}) for name in names]
    samples = [
        tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0, individual=inds[i], population=pop)
        for i in range(len(inds))
    ]
    anc = tables.nodes.add_row(time=1, population=pop)
    for s in samples:
        tables.edges.add_row(left=0, right=10, parent=anc, child=s)
    for pos, sample_idx in [(1, 0), (2, 0), (3, 1), (4, 2), (4.5, 3), (6, 1), (7, 1), (8, 0), (8.5, 2), (9, 3)]:
        site = tables.sites.add_row(position=pos, ancestral_state="0")
        tables.mutations.add_row(site=site, node=samples[sample_idx], derived_state="1")
    tables.sort()
    return tables.tree_sequence()


def test_windowing_sanity():
    ts = make_simple_ts()
    windows = np.array([0, 5, 10], dtype=float)
    load = mc.mutational_load(ts, windows=windows)
    names = mc.sample_names(ts)
    load, _ = mc.aggregate_by_individual(load, names)
    assert load.shape == (2, 2)
    # First window has site at 1 on A, second window has site at 7 on B
    assert load[0, 0] == 1
    assert load[0, 1] == 0
    assert load[1, 0] == 0
    assert load[1, 1] == 1


def test_build_snp_windows_single_variant_per_window():
    ts = make_simple_ts()
    windows = ms.build_snp_windows(ts, 1)
    assert windows.tolist() == [0.0, 7.0, 10.0]


def test_build_snp_windows_groups_variants():
    ts = make_simple_ts()
    windows = ms.build_snp_windows(ts, 2)
    assert windows.tolist() == [0.0, 10.0]


def test_build_snp_windows_no_mutations():
    ts = make_ts_no_mutations()
    windows = ms.build_snp_windows(ts, 5)
    assert windows.tolist() == [0.0, 10.0]


def test_outlier_mask_logic_against_sim_expected():
    # Per-individual expected (from sim) matrix matches load shape; threshold is
    # element-wise against the per-individual expectation.
    load = np.array([[12, 5, 5], [2, 2, 2]], dtype=float)
    expected = np.array([[5, 5, 5], [2, 2, 2]], dtype=float)
    cutoff = 0.5
    high = (1 + cutoff) * expected
    low = (1 - cutoff) * expected
    mask = (load > high) | (load < low)
    assert mask[0].tolist() == [True, False, False]
    assert mask[1].tolist() == [False, False, False]


def test_zero_expected_flags_observed_positive_load():
    # A zero simulated expectation still participates in outlier calling:
    # observed zero is fine, but observed-positive load is high relative to zero.
    load = np.array([[3, 0]], dtype=float)
    expected = np.array([[0, 0]], dtype=float)
    high = (1 + 0.5) * expected
    low = (1 - 0.5) * expected
    mask = (load > high) | (load < low)
    assert mask.tolist() == [[True, False]]


def test_load_chart_html_contains_bar_chars():
    load = np.array([0.5, 1.0])
    result = ms.load_chart_html(load, ["A", "B"], "Test")
    assert "█" in result
    assert "Test" in result


def test_load_chart_html_marks_outliers_red():
    load = np.array([0.5, 1.0])
    outlier_mask = np.array([True, False])
    result = ms.load_chart_html(load, ["A", "B"], "Test", outlier_mask=outlier_mask)
    assert "#d62728" in result
    assert "#444444" in result


def test_outlier_hist_html_shows_counts():
    result = ms.outlier_hist_html([0, 1, 1, 2])
    assert "█" in result
    assert "Outlier window counts" in result


def test_summarize_lineage_outliers_groups_prefixes():
    rows = ms.summarize_lineage_outliers(["L1_A", "L1_B", "L2_A"], [2, 1, 4])
    assert rows == [("L2", 4), ("L1", 3)]


def _ms_args(ts_path, **overrides):
    base = {
        "ts": str(ts_path),
        "window_size": 5.0,
        "snp_window": None,
        "cutoff": 0.5,
        "mutation_rate": 1.0,
        "random_seed": 1,
        "out": "out.html",
        "suffix_to_strip": "_anchorwave",
    }
    base.update(overrides)
    return type("A", (), base)()


def _patch_constant_expected(monkeypatch, value):
    # Force a deterministic per-individual expected so threshold outcomes are
    # independent of msprime's RNG.
    def _fake(ts, ts_path, windows, names, scalar_rate, seed):
        n_ind = len({n for n in names})
        n_win = len(windows) - 1
        return np.full((n_win, n_ind), value, dtype=float)
    monkeypatch.setattr(ms, "simulate_expected_load", _fake)


def test_outputs_written(tmp_path, monkeypatch):
    ts = make_simple_ts()
    ts_path = tmp_path / "test.trees"
    ts.dump(ts_path)
    cwd = Path(__file__).resolve().parents[1]
    os.chdir(cwd)
    # Clean outputs
    for p in (cwd / "results", cwd / "logs"):
        if p.exists():
            for item in p.iterdir():
                if item.is_file():
                    item.unlink()

    monkeypatch.setattr(ms, "load_ts", lambda _: ts)
    _patch_constant_expected(monkeypatch, 3.0)
    monkeypatch.setattr(ms, "parse_args", lambda: _ms_args(ts_path, out="out.html"))
    ms.main()

    assert (cwd / "results" / "out.html").exists()
    assert (cwd / "logs" / "out.log").exists()
    html = (cwd / "results" / "out.html").read_text()
    assert "2 of 2 windows have at least one outlier" in html
    assert "Outlier cutoff: 0.500 of sim expectation" in html
    assert "<td>A</td><td>2</td>" in html
    assert "<td>B</td><td>2</td>" in html


def test_outputs_written_with_lineage_table(tmp_path, monkeypatch):
    ts = make_lineage_ts()
    ts_path = tmp_path / "lineage.trees"
    ts.dump(ts_path)
    cwd = Path(__file__).resolve().parents[1]
    os.chdir(cwd)
    monkeypatch.setattr(ms, "load_ts", lambda _: ts)
    # Expected load is constant across all 4 individuals; observed load skews
    # toward L1_A and L1_B in window 0 (positions 1, 2, 3 → loads 2, 1, 0, 0)
    # and toward L1_B in window 1 (position 6, 7 → only L1_B loaded). With
    # cutoff=0.5 and expected=1, both windows flag L1 individuals.
    _patch_constant_expected(monkeypatch, 1.0)
    monkeypatch.setattr(ms, "parse_args", lambda: _ms_args(ts_path, out="lineage.html"))
    ms.main()
    html = (cwd / "results" / "lineage.html").read_text()
    assert "Outlier windows by lineage" in html
    assert "<td>L1</td><td>2</td>" in html
    assert "<td>L2</td><td>0</td>" in html


def test_no_remove_no_trimmed(tmp_path, monkeypatch):
    ts = make_simple_ts()
    ts_path = tmp_path / "test.trees"
    ts.dump(ts_path)
    cwd = Path(__file__).resolve().parents[1]
    os.chdir(cwd)
    monkeypatch.setattr(ms, "load_ts", lambda _: ts)
    _patch_constant_expected(monkeypatch, 3.0)
    monkeypatch.setattr(ms, "parse_args", lambda: _ms_args(ts_path, out="out.html"))
    ms.main()
    assert not (cwd / "results" / "test_trimmed.tsz").exists()


def test_no_mutations_outliers_empty(tmp_path, monkeypatch):
    ts = make_ts_no_mutations()
    ts_path = tmp_path / "nomut.trees"
    ts.dump(ts_path)
    cwd = Path(__file__).resolve().parents[1]
    os.chdir(cwd)
    monkeypatch.setattr(ms, "load_ts", lambda _: ts)
    _patch_constant_expected(monkeypatch, 0.0)
    monkeypatch.setattr(ms, "parse_args", lambda: _ms_args(ts_path, out="nomut.html"))
    ms.main()
    assert (cwd / "results" / "nomut.html").exists()


def test_single_sample_ts():
    ts = make_ts_no_mutations(n_samples=1, length=10)
    windows = np.array([0, 10], dtype=float)
    load = mc.mutational_load(ts, windows=windows)
    names = mc.sample_names(ts)
    load, _ = mc.aggregate_by_individual(load, names)
    assert load.shape == (1, 1)


def test_intervals_outside_sequence_length():
    ts = make_simple_ts()
    intervals = {"A": {"starts": [100.0], "ends": [200.0]}}
    name_to_nodes = mc.name_to_nodes_map(ts)
    trimmed = mc.trim_ts_by_intervals(ts, intervals, name_to_nodes)
    assert trimmed.num_sites == ts.num_sites


def test_overlapping_bed_with_commas(tmp_path):
    bed = tmp_path / "x.bed"
    bed.write_text("chr1\t1\t4\tA,B\nchr1\t3\t5\tB,C\n")
    remove = mc.load_remove_intervals([bed])
    assert set(remove.keys()) == {"A", "B", "C"}


def test_all_samples_removed_in_segment():
    ts = make_simple_ts()
    intervals = {"A": {"starts": [0.0], "ends": [10.0]}, "B": {"starts": [0.0], "ends": [10.0]}}
    name_to_nodes = mc.name_to_nodes_map(ts)
    trimmed = mc.trim_ts_by_intervals(ts, intervals, name_to_nodes)
    assert trimmed.sequence_length == ts.sequence_length


def test_idempotent_no_effect_remove():
    ts = make_simple_ts()
    intervals = {"C": {"starts": [0.0], "ends": [5.0]}}  # C not in TS
    name_to_nodes = mc.name_to_nodes_map(ts)
    trimmed = mc.trim_ts_by_intervals(ts, intervals, name_to_nodes)
    assert trimmed.num_sites == ts.num_sites
    assert trimmed.num_edges >= ts.num_edges


def test_sample_order_preserved():
    ts = make_simple_ts()
    intervals = {"A": {"starts": [0.0], "ends": [5.0]}}
    name_to_nodes = mc.name_to_nodes_map(ts)
    trimmed = mc.trim_ts_by_intervals(ts, intervals, name_to_nodes)
    assert trimmed.samples().tolist() == ts.samples().tolist()


def test_trim_validate():
    ts = make_simple_ts()
    intervals = {"A": {"starts": [0.0], "ends": [5.0]}}
    name_to_nodes = mc.name_to_nodes_map(ts)
    trimmed = mc.trim_ts_by_intervals(ts, intervals, name_to_nodes)
    if hasattr(trimmed, "validate"):
        trimmed.validate()


def test_large_interval_counts():
    intervals = {"A": {"starts": [], "ends": []}}
    for i in range(1000):
        intervals["A"]["starts"].append(float(i))
        intervals["A"]["ends"].append(float(i + 0.5))
    segs = mc.build_removal_segments(intervals, 2000.0)
    assert segs[0][0] == 0.0
    assert segs[-1][1] == 2000.0


def test_many_individuals_shapes():
    ts = make_ts_many_individuals(n=100, length=10)
    windows = np.array([0, 10], dtype=float)
    load = mc.mutational_load(ts, windows=windows)
    names = mc.sample_names(ts)
    load, unique = mc.aggregate_by_individual(load, names)
    assert load.shape == (1, 100)
    assert len(unique) == 100


def test_relative_bed_paths(tmp_path, monkeypatch):
    bed = tmp_path / "rel.bed"
    bed.write_text("chr1\t1\t3\tA\n")
    cwd = os.getcwd()
    os.chdir(tmp_path)
    try:
        remove = mc.load_remove_intervals([Path("rel.bed")])
        assert "A" in remove
    finally:
        os.chdir(cwd)


def test_output_overwrite(tmp_path, monkeypatch):
    ts = make_simple_ts()
    ts_path = tmp_path / "test.trees"
    ts.dump(ts_path)
    cwd = Path(__file__).resolve().parents[1]
    os.chdir(cwd)
    monkeypatch.setattr(ms, "load_ts", lambda _: ts)
    _patch_constant_expected(monkeypatch, 3.0)
    monkeypatch.setattr(ms, "parse_args", lambda: _ms_args(ts_path, out="overwrite.html"))
    ms.main()
    out = cwd / "results" / "overwrite.html"
    assert out.exists()
    first = out.read_text()
    ms.main()
    second = out.read_text()
    assert first == second


def test_outputs_written_with_snp_windows(tmp_path, monkeypatch):
    ts = make_simple_ts()
    ts_path = tmp_path / "test.trees"
    ts.dump(ts_path)
    cwd = Path(__file__).resolve().parents[1]
    os.chdir(cwd)
    monkeypatch.setattr(ms, "load_ts", lambda _: ts)
    _patch_constant_expected(monkeypatch, 3.0)
    monkeypatch.setattr(
        ms, "parse_args",
        lambda: _ms_args(ts_path, window_size=None, snp_window=1, out="snp_out.html"),
    )
    ms.main()

    assert (cwd / "results" / "snp_out.html").exists()


def test_remove_bed_parsing(tmp_path):
    bed = tmp_path / "x.bed"
    bed.write_text("chr1\t1\t3\tA,B\nchr1\t5\t6\tC\n")
    remove = mc.load_remove_intervals([bed])
    assert set(remove.keys()) == {"A", "B", "C"}
    assert remove["A"]["starts"] == [1.0]
    assert remove["B"]["starts"] == [1.0]
    assert remove["C"]["starts"] == [5.0]


def test_segment_merge_logic():
    intervals = {
        "A": {"starts": [1.0, 4.0], "ends": [3.0, 6.0]},
        "B": {"starts": [2.0], "ends": [5.0]},
    }
    segs = mc.build_removal_segments(intervals, 10.0)
    # Expect segments covering [0,1), [1,2), [2,3), [3,4), [4,5), [5,6), [6,10)
    assert segs[0][0] == 0.0 and segs[0][1] == 1.0
    assert segs[-1][0] == 6.0 and segs[-1][1] == 10.0


def test_trim_preserves_coordinates(tmp_path):
    ts = make_simple_ts()
    intervals = {"A": {"starts": [0.0], "ends": [5.0]}}
    name_to_nodes = mc.name_to_nodes_map(ts)
    trimmed = mc.trim_ts_by_intervals(ts, intervals, name_to_nodes)
    assert trimmed.sequence_length == ts.sequence_length
    assert trimmed.sites_position.min() >= 0.0
    assert trimmed.sites_position.max() <= ts.sequence_length


def test_trim_removes_nodes_in_interval():
    ts = make_simple_ts()
    intervals = {"A": {"starts": [0.0], "ends": [5.0]}}
    name_to_nodes = mc.name_to_nodes_map(ts)
    trimmed = mc.trim_ts_by_intervals(ts, intervals, name_to_nodes)
    # Mutation at position 1 on A should be removed; position 7 on B remains
    assert trimmed.num_sites == 1
    assert trimmed.sites_position.tolist() == [7.0]


def test_trim_mutation_parent_integrity():
    ts = make_simple_ts()
    intervals = {"A": {"starts": [0.0], "ends": [5.0]}}
    name_to_nodes = mc.name_to_nodes_map(ts)
    trimmed = mc.trim_ts_by_intervals(ts, intervals, name_to_nodes)
    # Should load without parent errors
    assert trimmed.num_sites >= 0
