from pathlib import Path

import sys

import pytest
import tskit

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import merge_treefiles_by_replicate as merger
from argtest_common import (
    chrom_position_from_genome,
    genome_position,
    tree_at_chrom_position,
)


def make_simple_ts(length=10, site_pos=1):
    tables = tskit.TableCollection(sequence_length=length)
    pop = tables.populations.add_row()
    n0 = tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0, population=pop)
    n1 = tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0, population=pop)
    anc = tables.nodes.add_row(time=1, population=pop)
    tables.edges.add_row(left=0, right=length, parent=anc, child=n0)
    tables.edges.add_row(left=0, right=length, parent=anc, child=n1)
    tables.sort()
    return tables.tree_sequence()


def make_merged(tmp_path):
    ts1 = make_simple_ts(length=5)
    ts2 = make_simple_ts(length=7)
    p1 = tmp_path / "base.chr1.1.trees"
    p2 = tmp_path / "base.chr2.1.trees"
    ts1.dump(p1)
    ts2.dump(p2)
    merged, _ = merger.merge_group([("chr1", p1), ("chr2", p2)])
    return merged


def test_genome_position_maps_chrom_offset(tmp_path):
    merged = make_merged(tmp_path)
    assert genome_position(merged, "chr1", 3) == 3.0
    assert genome_position(merged, "chr2", 0) == 5.0
    assert genome_position(merged, "chr2", 2) == 7.0


def test_genome_position_roundtrip(tmp_path):
    merged = make_merged(tmp_path)
    gpos = genome_position(merged, "chr2", 4)
    assert chrom_position_from_genome(merged, gpos) == ("chr2", 4.0)


def test_genome_position_out_of_range(tmp_path):
    merged = make_merged(tmp_path)
    with pytest.raises(ValueError):
        genome_position(merged, "chr1", 5)  # length is 5, so [0, 5)
    with pytest.raises(ValueError):
        genome_position(merged, "chr2", -1)
    with pytest.raises(ValueError):
        genome_position(merged, "chr1", float("nan"))


def test_genome_position_unknown_chrom(tmp_path):
    merged = make_merged(tmp_path)
    with pytest.raises(KeyError):
        genome_position(merged, "chr9", 0)


def test_genome_position_missing_metadata():
    ts = make_simple_ts(length=5)  # no chrom_offsets metadata
    with pytest.raises(KeyError):
        genome_position(ts, "chr1", 0)


def test_tree_at_chrom_position(tmp_path):
    merged = make_merged(tmp_path)
    tree = tree_at_chrom_position(merged, "chr2", 3)
    assert tree.interval.left <= 8.0 < tree.interval.right
