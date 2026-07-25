from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import tskit

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import coalescence_ne_plots_from_ts as coal


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


def test_connected_pair_quantiles_condition_out_isolated_pair_spans():
    ts = _partially_isolated_ts()
    quantiles = np.linspace(0, 1, 5)

    standard = ts.pair_coalescence_quantiles(quantiles)
    assert np.isnan(standard[1:-1]).all()

    conditional = coal.connected_pair_coalescence_quantiles(ts, quantiles)
    assert np.isfinite(conditional).all()
    assert np.all(np.diff(conditional) > 0)
    np.testing.assert_array_equal(conditional[[0, -1]], [1, 8])
