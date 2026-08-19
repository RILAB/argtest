from __future__ import annotations

from pathlib import Path

import numpy as np
import tskit


import coalescence_ne_plots_from_ts as coal


class _RecordingAxes:
    def __init__(self):
        self.calls = []

    def step(self, x, y, **kwargs):
        self.calls.append((np.asarray(x), np.asarray(y), kwargs))


def test_find_tree_files_uses_natural_replicate_order(tmp_path):
    for replicate in [10, 2, 1, 0]:
        (tmp_path / f"run.combined.{replicate}.tsz").touch()

    files = coal.find_tree_files(tmp_path, "*.tsz")

    assert [path.name for path in files] == [
        "run.combined.0.tsz",
        "run.combined.1.tsz",
        "run.combined.2.tsz",
        "run.combined.10.tsz",
    ]


def _partially_isolated_ts():
    """Eight equal-span trees with only one of six sample pairs connected."""
    tables = tskit.TableCollection(sequence_length=8)
    pop = tables.populations.add_row()
    samples = [
        tables.nodes.add_row(
            flags=tskit.NODE_IS_SAMPLE,
            time=0,
            population=pop,
        )
        for _ in range(4)
    ]
    for j in range(8):
        ancestor = tables.nodes.add_row(time=j + 1, population=pop)
        for sample in samples[:2]:
            tables.edges.add_row(j, j + 1, ancestor, sample)
    tables.sort()
    tables.build_index()
    return tables.tree_sequence()


def test_native_pair_quantiles_condition_out_isolated_pair_spans():
    ts = _partially_isolated_ts()
    quantiles = np.linspace(0, 1, 5)

    result = ts.pair_coalescence_quantiles(quantiles)
    assert np.isfinite(result).all()
    assert np.all(np.diff(result) > 0)
    np.testing.assert_array_equal(result[[0, -1]], [1, 8])


def test_compute_pair_coal_uses_native_partial_missing_normalisation():
    ts = _partially_isolated_ts()
    time_windows = np.array([0.0, 2.0, 4.0, 6.0, 8.0, np.inf])

    pdf, rates = coal.compute_pair_coal(ts, time_windows, tail_cutoff=1e-12)

    expected_pdf = ts.pair_coalescence_counts(
        time_windows=time_windows,
        pair_normalise=True,
    )
    expected_rates = ts.pair_coalescence_rates(time_windows=time_windows)
    np.testing.assert_allclose(pdf, expected_pdf)
    np.testing.assert_allclose(rates, expected_rates, equal_nan=True)


def test_replicate_plot_excludes_burnin_values():
    ax = _RecordingAxes()
    x = np.array([1, 2])
    values = np.array([[10, 20], [30, 40], [50, 60]])

    coal.plot_postburn_replicates(
        ax,
        x,
        values,
        np.array([1, 2]),
        color="gray",
    )

    assert len(ax.calls) == 2
    np.testing.assert_array_equal(ax.calls[0][1], values[1])
    np.testing.assert_array_equal(ax.calls[1][1], values[2])


def test_write_estimates_saves_only_postburn_replicates_and_mean(tmp_path):
    out = tmp_path / "estimates.tsv"
    values = np.array([[1.0, 2.0], [3.0, 4.0]])

    coal.write_coalescence_estimates(
        out,
        [Path("rep0.tsz"), Path("rep1.tsz")],
        np.array([1]),
        np.array([0.0, 1.0, 2.0, np.inf]),
        np.array([False, True, True]),
        values,
        values + 10,
        values + 20,
        values + 30,
        np.array([3.0, 4.0]),
        np.array([13.0, 14.0]),
        np.array([23.0, 24.0]),
        np.array([33.0, 34.0]),
    )

    rows = [line.split("\t") for line in out.read_text().splitlines()]
    assert len(rows) == 5
    assert [row[0] for row in rows[1:]] == [
        "replicate",
        "replicate",
        "posterior_mean",
        "posterior_mean",
    ]
    assert {row[1] for row in rows[1:3]} == {"1"}
    assert {row[2] for row in rows[1:3]} == {"rep1.tsz"}
    assert rows[1][4:6] == ["1", "2"]
