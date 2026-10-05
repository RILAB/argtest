import os
from pathlib import Path

import numpy as np
import pytest
import tskit

import argtest_common as mc
import mutload_summary as ms
import trim_samples as tsamp


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


def make_stacked_mutation_ts():
    """3 haploid individuals A, B, C (nodes 0, 1, 2).

    Tree over [0, 10): root 4 -> (n0, node 3), node 3 -> (n1, n2).
    Site at 1: private mutation on n0 (A).
    Site at 7: 0->1 on node 3 (shared by B and C), then a back mutation 1->0
    on n2, whose parent is the node-3 mutation.
    """
    tables = tskit.TableCollection(sequence_length=10)
    tables.individuals.metadata_schema = tskit.MetadataSchema.permissive_json()
    pop = tables.populations.add_row()
    inds = [tables.individuals.add_row(metadata={"id": name}) for name in "ABC"]
    n0, n1, n2 = (
        tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0, individual=i, population=pop)
        for i in inds
    )
    mid = tables.nodes.add_row(time=1, population=pop)
    root = tables.nodes.add_row(time=2, population=pop)
    tables.edges.add_row(left=0, right=10, parent=mid, child=n1)
    tables.edges.add_row(left=0, right=10, parent=mid, child=n2)
    tables.edges.add_row(left=0, right=10, parent=root, child=n0)
    tables.edges.add_row(left=0, right=10, parent=root, child=mid)
    s1 = tables.sites.add_row(position=1, ancestral_state="0")
    tables.mutations.add_row(site=s1, node=n0, derived_state="1")
    s7 = tables.sites.add_row(position=7, ancestral_state="0")
    m_mid = tables.mutations.add_row(site=s7, node=mid, derived_state="1", time=1.5)
    tables.mutations.add_row(site=s7, node=n2, derived_state="0", parent=m_mid, time=0.5)
    tables.sort()
    tables.build_index()
    tables.compute_mutation_parents()
    return tables.tree_sequence()


@pytest.fixture
def summary_root(tmp_path, monkeypatch):
    """Redirect mutload_summary's repo-root results/ and logs/ into tmp_path.

    ``ms.main`` derives its output dirs from ``Path(__file__).parent.parent``,
    so pointing ``__file__`` at ``tmp_path/scripts/...`` keeps the tests from
    writing into (or cleaning) the real repo's results/ and logs/.
    """
    monkeypatch.setattr(ms, "load_ts", lambda path: tskit.load(path))
    monkeypatch.setattr(ms, "__file__", str(tmp_path / "scripts" / "mutload_summary.py"))
    return tmp_path


def _run_summary(monkeypatch, root, ts, out, **overrides):
    ts_path = root / f"{Path(out).stem}.trees"
    ts.dump(ts_path)
    monkeypatch.setattr(ms, "parse_args", lambda: _ms_args(ts_path, out=out, **overrides))
    ms.main()
    return (root / "results" / out).read_text()


def test_windowing_sanity():
    ts = make_simple_ts()
    windows = np.array([0, 5, 10], dtype=float)
    load = mc.mutational_load(ts, windows=windows)
    names = mc.sample_names(ts)
    load, unique = mc.aggregate_by_individual(load, names)
    assert unique == ["A", "B"]
    assert load.shape == (2, 2)
    # First window has site at 1 on A, second window has site at 7 on B
    assert load[0, 0] == 1
    assert load[0, 1] == 0
    assert load[1, 0] == 0
    assert load[1, 1] == 1


def test_trim_samples_single_pass_removes_target_node_mutations():
    ts = make_simple_ts()
    intervals = {"A": {"starts": [0.0], "ends": [5.0]}}

    trimmed, summary = tsamp.trim_samples_single_pass(ts, intervals)

    assert summary["names_removed"] == {"A"}
    assert trimmed.sites_position.tolist() == [7.0]


def test_trim_samples_single_pass_preserves_structure():
    ts = make_stacked_mutation_ts()
    intervals = {"A": {"starts": [0.0], "ends": [5.0]}}

    trimmed, summary = tsamp.trim_samples_single_pass(ts, intervals)

    assert summary == {
        "names_removed": {"A"},
        "intervals_applied": 1,
        "sample_nodes_removed": 1,
    }
    # Sample ids, order and names survive simplify.
    assert trimmed.samples().tolist() == [0, 1, 2]
    assert mc.sample_names(trimmed) == ["A", "B", "C"]

    # A is isolated inside the removed interval and attached outside it.
    inside = trimmed.at(2.0)
    assert inside.parent(0) == tskit.NULL and inside.num_children(0) == 0
    outside = trimmed.at(7.0)
    assert outside.parent(0) != tskit.NULL
    # B and C keep their shared parent throughout.
    for pos in (2.0, 7.0):
        tree = trimmed.at(pos)
        assert tree.parent(1) == tree.parent(2) != tskit.NULL

    # A's private site inside the interval is gone; the stacked site remains.
    assert trimmed.sites_position.tolist() == [7.0]
    site = trimmed.site(0)
    assert [m.derived_state for m in site.mutations] == ["1", "0"]
    shared, back = site.mutations
    assert back.parent == shared.id
    assert shared.parent == tskit.NULL
    assert next(trimmed.variants()).genotypes.tolist() == [0, 1, 0]

    # Mutation parents are consistent with the trimmed topology.
    recomputed = trimmed.dump_tables()
    recomputed.compute_mutation_parents()
    assert recomputed.mutations.parent.tolist() == trimmed.tables.mutations.parent.tolist()


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


def test_outside_band_is_elementwise_against_expected():
    # Each cell is compared to its own expectation, high and low.
    load = np.array([[12, 5, 1], [2, 2, 2]], dtype=float)
    expected = np.array([[5, 5, 5], [2, 2, 4]], dtype=float)
    mask = mc.outside_band(load, expected, 0.5)
    assert mask.tolist() == [[True, False, True], [False, False, False]]


def test_outside_band_edges_are_not_flagged():
    # Band is [2.5, 7.5]; values exactly on either edge stay inside.
    expected = np.array([5.0, 5.0, 5.0, 5.0])
    load = np.array([2.5, 7.5, 2.4, 7.6])
    assert mc.outside_band(load, expected, 0.5).tolist() == [False, False, True, True]


def test_zero_expected_flags_observed_positive_load():
    # A zero simulated expectation still participates in outlier calling:
    # observed zero is fine, but observed-positive load is high relative to zero.
    load = np.array([[3, 0]], dtype=float)
    expected = np.array([[0, 0]], dtype=float)
    assert mc.outside_band(load, expected, 0.5).tolist() == [[True, False]]


def test_load_chart_html_contains_bar_chars():
    load = np.array([0.5, 1.0])
    result = ms.load_chart_html(load, ["A", "B"], "Test")
    assert "█" in result
    assert "Test" in result


def test_load_chart_html_marks_outliers_red():
    load = np.array([0.5, 1.0])
    outlier_mask = np.array([True, False])
    result = ms.load_chart_html(load, ["A", "B"], "Test", outlier_mask=outlier_mask)
    lines = result.splitlines()
    row_a = next(ln for ln in lines if ">A<" in ln)
    row_b = next(ln for ln in lines if ">B<" in ln)
    # The flagged individual's whole row is red; the other row is not.
    assert "#d62728" in row_a and "#444444" not in row_a
    assert "#444444" in row_b and "#d62728" not in row_b


def test_summarize_lineage_flags_groups_prefixes():
    rows = ms.summarize_lineage_flags(
        ["L1_A", "L1_B", "L2_A"],
        [True, False, True],
    )
    # Sorted by flagged count desc, then by lineage name asc.
    assert rows == [("L1", 2, 1), ("L2", 1, 1)]


def _ms_args(ts_path, **overrides):
    base = {
        "ts": str(ts_path),
        "window_size": 5.0,
        "snp_window": None,
        "cutoff": 0.5,
        "mutation_rate": 1.0,
        "random_seed": 1,
        "out": "out.html",
        "name_substring_to_remove": "_anchorwave",
    }
    base.update(overrides)
    return type("A", (), base)()


def _patch_constant_expected(monkeypatch, value):
    # Force a deterministic per-individual expected so threshold outcomes are
    # independent of msprime's RNG.
    def _fake(ts, windows, names, mutation_rate, seed):
        n_ind = len({n for n in names})
        n_win = len(windows) - 1
        return np.full((n_win, n_ind), value, dtype=float)
    monkeypatch.setattr(ms, "simulate_expected_load", _fake)


def test_outputs_written(summary_root, monkeypatch):
    _patch_constant_expected(monkeypatch, 3.0)
    html = _run_summary(monkeypatch, summary_root, make_simple_ts(), "out.html")

    assert (summary_root / "logs" / "out.log").exists()
    # Two windows × two individuals = 4 (window, individual) pairs.
    # Per-window expected = 3, band = [1.5, 4.5]. Observed per cell is 1 or 0,
    # all below 1.5, so all 4 pairs are flagged for trimming.
    assert "4 of 4 (window, individual) pairs flagged for trimming" in html
    # With everything pruned, residual obs and exp are zero for each individual
    # → residual flag never fires.
    assert "All 2 individuals within the cutoff band after pruning" in html
    assert "Outlier cutoff: 0.500 of sim expectation" in html
    # Lineage table shows flagged + total per lineage. No residual flags.
    assert "<td>A</td><td>0</td><td>1</td>" in html
    assert "<td>B</td><td>0</td><td>1</td>" in html


def test_outputs_written_with_lineage_table(summary_root, monkeypatch):
    # Window 0 (pos 0-5): L1_A=2, L1_B=1, L2_A=1, L2_B=1.
    # Window 1 (pos 5-10): L1_A=1, L1_B=2, L2_A=1, L2_B=1.
    # Per-window expected = 0.75, band = [0.375, 1.125]. L1_A's 2 in W0 and
    # L1_B's 2 in W1 are above the band → those pairs flag. After pruning,
    # residuals lie inside the cutoff band by construction.
    _patch_constant_expected(monkeypatch, 0.75)
    html = _run_summary(monkeypatch, summary_root, make_lineage_ts(), "lineage.html")
    assert "Flagged individuals by lineage" in html
    assert "2 of 8 (window, individual) pairs flagged for trimming" in html
    assert "All 4 individuals within the cutoff band after pruning" in html
    # Residual flag never fires here → all lineage flagged counts are 0.
    assert "<td>L1</td><td>0</td><td>2</td>" in html
    assert "<td>L2</td><td>0</td><td>2</td>" in html


def test_summary_writes_only_report_and_log(summary_root, monkeypatch):
    # mutload_summary is report-only: it must not emit a trimmed tree sequence.
    _patch_constant_expected(monkeypatch, 3.0)
    _run_summary(monkeypatch, summary_root, make_simple_ts(), "out.html")
    assert sorted(p.name for p in (summary_root / "results").iterdir()) == ["out.html"]
    assert sorted(p.name for p in (summary_root / "logs").iterdir()) == ["out.log"]


def test_no_mutations_outliers_empty(summary_root, monkeypatch):
    # Zero observed and zero expected everywhere: nothing is outside the band.
    _patch_constant_expected(monkeypatch, 0.0)
    html = _run_summary(monkeypatch, summary_root, make_ts_no_mutations(), "nomut.html")
    assert "0 of 4 (window, individual) pairs flagged for trimming" in html
    assert "All 2 individuals within the cutoff band after pruning" in html


def test_single_sample_ts():
    ts = make_ts_no_mutations(n_samples=1, length=10)
    windows = np.array([0, 10], dtype=float)
    load = mc.mutational_load(ts, windows=windows)
    names = mc.sample_names(ts)
    load, unique = mc.aggregate_by_individual(load, names)
    assert unique == ["I0"]
    assert load.tolist() == [[0]]


def test_overlapping_bed_with_commas(tmp_path):
    bed = tmp_path / "x.bed"
    bed.write_text("chr1\t1\t4\tA,B\nchr1\t3\t5\tB,C\n")
    remove = mc.load_remove_intervals([bed])
    assert remove == {
        "A": {"starts": [1.0], "ends": [4.0]},
        # Overlapping spans are kept separate here; merging happens in trim.
        "B": {"starts": [1.0, 3.0], "ends": [4.0, 5.0]},
        "C": {"starts": [3.0], "ends": [5.0]},
    }


def test_many_individuals_shapes():
    ts = make_ts_many_individuals(n=100, length=10)
    windows = np.array([0, 10], dtype=float)
    load = mc.mutational_load(ts, windows=windows)
    names = mc.sample_names(ts)
    load, unique = mc.aggregate_by_individual(load, names)
    assert unique == [f"I{i}" for i in range(100)]
    assert load.shape == (1, 100)
    # Only I0 carries a mutation.
    assert load[0, 0] == 1 and load[0, 1:].sum() == 0


def test_relative_bed_paths(tmp_path, monkeypatch):
    bed = tmp_path / "rel.bed"
    bed.write_text("chr1\t1\t3\tA\n")
    cwd = os.getcwd()
    os.chdir(tmp_path)
    try:
        remove = mc.load_remove_intervals([Path("rel.bed")])
        assert remove == {"A": {"starts": [1.0], "ends": [3.0]}}
    finally:
        os.chdir(cwd)


def test_output_overwrite(summary_root, monkeypatch):
    _patch_constant_expected(monkeypatch, 3.0)
    first = _run_summary(monkeypatch, summary_root, make_simple_ts(), "overwrite.html")
    out = summary_root / "results" / "overwrite.html"
    ms.main()
    second = out.read_text()
    assert first == second


def test_outputs_written_with_snp_windows(summary_root, monkeypatch):
    _patch_constant_expected(monkeypatch, 3.0)
    html = _run_summary(
        monkeypatch, summary_root, make_simple_ts(), "snp_out.html",
        window_size=None, snp_window=1,
    )
    # One SNP per window -> windows [0, 7) and [7, 10): 2 windows x 2 individuals,
    # each observed 0 or 1 against a band of [1.5, 4.5], so all 4 flag.
    assert "4 of 4 (window, individual) pairs flagged for trimming" in html


def test_remove_bed_parsing(tmp_path):
    bed = tmp_path / "x.bed"
    bed.write_text("chr1\t1\t3\tA,B\nchr1\t5\t6\tC\n")
    remove = mc.load_remove_intervals([bed])
    assert set(remove.keys()) == {"A", "B", "C"}
    assert remove["A"]["starts"] == [1.0]
    assert remove["B"]["starts"] == [1.0]
    assert remove["C"]["starts"] == [5.0]
