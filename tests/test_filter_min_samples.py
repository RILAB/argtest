from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest
import tskit

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import filter_min_samples as fms

# Reuse the integration-test dataset/runner for the Snakemake-level checks.
from test_pipeline_integration import (  # noqa: E402
    REPO_ROOT,
    SNK_AVAILABLE,
    SNK_SKIP_REASON,
    build_dataset,
    run_snakemake,
)


def _ts_with_isolated_middle(length=30):
    """Build a 3-tree ts: full samples on [0,10) and [20,30), but on [10,20)
    two of the four samples are isolated (no parent edge).

    A recombination breakpoint is created by giving sample nodes different
    parent edges over the three intervals.
    """
    tables = tskit.TableCollection(sequence_length=length)
    pop = tables.populations.add_row()
    s = [tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0, population=pop) for _ in range(4)]
    anc = tables.nodes.add_row(time=1, population=pop)
    anc2 = tables.nodes.add_row(time=1, population=pop)

    # Interval [0,10): all 4 samples attached to anc.
    for u in s:
        tables.edges.add_row(left=0, right=10, parent=anc, child=u)
    # Interval [10,20): only s0, s1 attached (s2, s3 isolated -> 2 retained).
    tables.edges.add_row(left=10, right=20, parent=anc, child=s[0])
    tables.edges.add_row(left=10, right=20, parent=anc, child=s[1])
    # Interval [20,30): all 4 samples attached to anc2.
    for u in s:
        tables.edges.add_row(left=20, right=30, parent=anc2, child=u)

    tables.sort()
    tables.build_index()
    return tables.tree_sequence()


def test_retained_counts_full_partial_and_no_edges():
    ts = _ts_with_isolated_middle()
    counts = fms.retained_sample_intervals(ts)
    assert counts == [(0.0, 10.0, 4), (10.0, 20.0, 2), (20.0, 30.0, 4)]

    # A tree with no edges must count 0 retained samples.
    tables = tskit.TableCollection(sequence_length=5)
    pop = tables.populations.add_row()
    tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0, population=pop)
    tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0, population=pop)
    tables.build_index()
    empty = tables.tree_sequence()
    assert fms.retained_sample_intervals(empty) == [(0.0, 5.0, 0)]


def _ts_with_tie_and_gap(length=40):
    """3 samples exercising coinciding edge endpoints and a zero-coverage tail.

    Coverage: [0,5)=2, then three adjacent unit-coverage spans that meet at
    coinciding endpoints (s0 ends / s1 starts at 10; s1 ends / s0 starts at 20)
    which must merge into one [5,30) segment, then a [30,40) span with no
    covering sample edge (coverage 0).
    """
    tables = tskit.TableCollection(sequence_length=length)
    pop = tables.populations.add_row()
    s = [tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0, population=pop) for _ in range(3)]
    anc = tables.nodes.add_row(time=1, population=pop)
    tables.edges.add_row(0, 10, anc, s[0])
    tables.edges.add_row(20, 30, anc, s[0])
    tables.edges.add_row(10, 20, anc, s[1])  # coincides with s0 endpoints at 10 and 20
    tables.edges.add_row(0, 5, anc, s[2])
    tables.sort()
    tables.build_index()
    return tables.tree_sequence()


def test_coverage_handles_ties_and_zero_gap():
    ts = _ts_with_tie_and_gap()
    # Coinciding endpoints must not spuriously spike/dip coverage, and adjacent
    # equal-coverage spans merge; the edge-less tail counts 0.
    assert fms.retained_sample_intervals(ts) == [
        (0.0, 5.0, 2),
        (5.0, 30.0, 1),
        (30.0, 40.0, 0),
    ]
    counts = fms.retained_sample_intervals(ts)
    assert fms.low_sample_intervals(counts, min_samples=2) == [[5.0, 40.0]]


def test_low_sample_intervals_merges_adjacent():
    counts = [
        (0.0, 10.0, 4),
        (10.0, 15.0, 1),
        (15.0, 20.0, 2),  # both below 3 -> merge into one [10,20)
        (20.0, 30.0, 4),
    ]
    assert fms.low_sample_intervals(counts, min_samples=3) == [[10.0, 20.0]]
    # Nothing below threshold.
    assert fms.low_sample_intervals(counts, min_samples=1) == []


def test_filter_drops_expected_span_and_preserves_coordinates():
    ts = _ts_with_isolated_middle()
    filtered, dropped, counts, dropped_bp = fms.filter_min_samples(ts, min_samples=3)

    assert dropped == [[10.0, 20.0]]
    assert dropped_bp == 10.0
    # Coordinates preserved (no trim/compaction).
    assert filtered.sequence_length == ts.sequence_length == 30.0
    # The dropped span is now an empty (edge-less) tree.
    mid_tree = filtered.at(15.0)
    assert mid_tree.num_edges == 0
    # Flanks retained.
    assert filtered.at(5.0).num_edges > 0
    assert filtered.at(25.0).num_edges > 0
    # Excision is structurally correct: the flanks keep all four samples while
    # the dropped middle span retains none.
    assert fms.retained_sample_intervals(filtered) == [
        (0.0, 10.0, 4),
        (10.0, 20.0, 0),
        (20.0, 30.0, 4),
    ]


def test_no_drops_returns_input_unchanged():
    ts = _ts_with_isolated_middle()
    filtered, dropped, counts, dropped_bp = fms.filter_min_samples(ts, min_samples=2)
    assert dropped == []
    assert dropped_bp == 0.0
    assert filtered is ts


def test_kept_intervals_intersected_with_drops():
    ts = _ts_with_isolated_middle()
    # Seed a prior kept_intervals metadata covering the whole sequence.
    tables = ts.dump_tables()
    tables.metadata_schema = tskit.MetadataSchema({"codec": "json"})
    tables.metadata = {"kept_intervals": [[0.0, 30.0]]}
    ts = tables.tree_sequence()

    filtered, dropped, _, _ = fms.filter_min_samples(ts, min_samples=3)
    assert dropped == [[10.0, 20.0]]
    kept = filtered.metadata["kept_intervals"]
    assert kept == [[0.0, 10.0], [20.0, 30.0]]


def test_kept_intervals_intersection_helper():
    # keep [0,100); drop [10,20) and [50,60) -> three surviving fragments.
    out = fms._intersect_keep_with_drops(
        [[0.0, 100.0]], [[10.0, 20.0], [50.0, 60.0]], 100.0
    )
    assert out == [[0.0, 10.0], [20.0, 50.0], [60.0, 100.0]]

    # Prior keep already excludes part of a drop: keep [0,15); drop [10,20).
    out = fms._intersect_keep_with_drops([[0.0, 15.0]], [[10.0, 20.0]], 30.0)
    assert out == [[0.0, 10.0]]


def test_main_writes_ts_and_bed(tmp_path, monkeypatch):
    ts = _ts_with_isolated_middle()
    ts_path = tmp_path / "chrX.1.trees"
    ts.dump(ts_path)
    out_ts = tmp_path / "out" / "chrX.1.trees"
    out_bed = tmp_path / "out" / "chrX.1.low_sample.bed"

    monkeypatch.setattr(
        fms,
        "parse_args",
        lambda: type(
            "A",
            (),
            {
                "ts": ts_path,
                "min_samples": 3,
                "out": out_ts,
                "out_mask": out_bed,
                "chrom": "chrX",
                "log": None,
            },
        )(),
    )
    fms.main()

    assert out_ts.exists()
    assert out_bed.exists()
    lines = out_bed.read_text().strip().splitlines()
    assert lines == ["chrX\t10\t20\t2\t3"]
    out = tskit.load(out_ts)
    assert out.sequence_length == 30.0
    assert out.at(15.0).num_edges == 0


# --- Snakemake-level routing checks ------------------------------------------


@pytest.mark.skipif(not SNK_AVAILABLE, reason=SNK_SKIP_REASON)
def test_min_samples_null_skips_step(tmp_path):
    dataset = build_dataset(tmp_path)
    # No min_samples key -> step 5b must not appear and merge reads step 5.
    result = run_snakemake(REPO_ROOT, dataset["config"], "-n", "-p")
    assert "step5b_filter_min_samples" not in result.stdout
    assert "step5b_min_samples" not in result.stdout


@pytest.mark.skipif(not SNK_AVAILABLE, reason=SNK_SKIP_REASON)
def test_min_samples_enabled_inserts_step(tmp_path):
    dataset = build_dataset(tmp_path)
    with open(dataset["config"], "a") as fh:
        fh.write("min_samples: 2\n")
    result = run_snakemake(REPO_ROOT, dataset["config"], "-n", "-p")
    assert "step5b_filter_min_samples" in result.stdout
    # Merge/step6 must now consume the step-5b output directory.
    assert "step5b_min_samples" in result.stdout
