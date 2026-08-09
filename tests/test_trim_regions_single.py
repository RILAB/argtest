import sys
from pathlib import Path

import tskit


import trim_regions_single as trs


def test_mask_loading_and_complement_clip_and_merge(tmp_path):
    mask_path = tmp_path / "mask.bed"
    mask_path.write_text("chr1\t-2\t3\nchr1\t2\t5\nchr1\t9\t12\n")

    masked = trs.load_mask_intervals(mask_path, 10)

    assert masked == [[0.0, 5.0], [9.0, 10.0]]
    assert trs.complement_intervals(masked, 10) == [[5.0, 9.0]]


def test_trim_preserves_embedded_input_ratemap(tmp_path, monkeypatch):
    tables = tskit.TableCollection(sequence_length=10)
    parent = tables.nodes.add_row(time=1)
    child = tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0)
    tables.edges.add_row(0, 10, parent=parent, child=child)
    tables.metadata_schema = tskit.MetadataSchema({"codec": "json"})
    tables.metadata = {
        "mu_position": [0.0, 4.0, 10.0],
        "mu_rate": [0.0, 1e-8],
    }
    tables.sort()

    input_path = tmp_path / "input.trees"
    output_path = tmp_path / "output.trees"
    mask_path = tmp_path / "empty.bed"
    tables.tree_sequence().dump(input_path)
    mask_path.write_text("")

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "trim_regions_single.py",
            "--ts",
            str(input_path),
            "--remove",
            str(mask_path),
            "--out",
            str(output_path),
        ],
    )
    trs.main()

    metadata = tskit.load(output_path).metadata
    assert metadata["mu_position"] == [0.0, 4.0, 10.0]
    assert metadata["mu_rate"] == [0.0, 1e-8]


def test_trim_preserves_unrelated_metadata_and_embeds_scalar_rate(tmp_path, monkeypatch):
    tables = tskit.TableCollection(sequence_length=10)
    parent = tables.nodes.add_row(time=1)
    child = tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0)
    tables.edges.add_row(0, 10, parent=parent, child=child)
    tables.metadata_schema = tskit.MetadataSchema({"codec": "json"})
    tables.metadata = {"provenance_label": "keep-me"}
    tables.sort()

    input_path = tmp_path / "input.trees"
    output_path = tmp_path / "nested" / "output.trees"
    mask_path = tmp_path / "mask.bed"
    tables.tree_sequence().dump(input_path)
    mask_path.write_text("chr1\t2\t4\n")

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "trim_regions_single.py",
            "--ts", str(input_path),
            "--remove", str(mask_path),
            "--out", str(output_path),
            "--mutation-rate", "2e-8",
        ],
    )
    trs.main()

    metadata = tskit.load(output_path).metadata
    assert metadata["provenance_label"] == "keep-me"
    assert metadata["kept_intervals"] == [[0.0, 2.0], [4.0, 10.0]]
    assert metadata["mu_position"] == [0.0, 10.0]
    assert metadata["mu_rate"] == [2e-8]
    # The scalar fallback must be recorded, not silent: it replaces a spatially
    # varying map with a flat rate and so removes step 3's local-rate correction.
    assert metadata["mu_source"] == {"kind": "scalar", "rate": 2e-8}


def test_sibling_map_wins_over_scalar_and_is_recorded(tmp_path, monkeypatch):
    """A discoverable map takes precedence and the scalar is never consulted."""
    import pickle

    import msprime

    from argtest_common import resolve_mu_rate_with_source

    tables = tskit.TableCollection(sequence_length=10)
    parent = tables.nodes.add_row(time=1)
    child = tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0)
    tables.edges.add_row(0, 10, parent=parent, child=child)
    tables.sort()

    chrom_dir = tmp_path / "chr1"
    chrom_dir.mkdir()
    input_path = chrom_dir / "7.trees"
    tables.tree_sequence().dump(input_path)
    map_path = chrom_dir / "chr1.mut_rate.p"
    with open(map_path, "wb") as fh:
        pickle.dump(msprime.RateMap(position=[0, 5, 10], rate=[1e-8, 5e-8]), fh)

    rate, source = resolve_mu_rate_with_source(
        tskit.load(input_path), input_path, scalar_fallback=2e-8
    )
    assert source == {"kind": "sibling", "path": str(map_path)}
    assert list(rate.rate) == [1e-8, 5e-8]

    mask_path = tmp_path / "empty.bed"
    mask_path.write_text("")
    output_path = tmp_path / "out.trees"
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "trim_regions_single.py",
            "--ts", str(input_path),
            "--remove", str(mask_path),
            "--out", str(output_path),
            "--mutation-rate", "2e-8",
        ],
    )
    trs.main()

    metadata = tskit.load(output_path).metadata
    assert metadata["mu_source"]["kind"] == "sibling"
    assert metadata["mu_rate"] == [1e-8, 5e-8]  # map preserved, scalar unused


def test_missing_map_and_no_scalar_fails_loudly(tmp_path, monkeypatch):
    import pytest

    tables = tskit.TableCollection(sequence_length=10)
    parent = tables.nodes.add_row(time=1)
    child = tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0)
    tables.edges.add_row(0, 10, parent=parent, child=child)
    tables.sort()

    input_path = tmp_path / "input.trees"
    tables.tree_sequence().dump(input_path)
    mask_path = tmp_path / "empty.bed"
    mask_path.write_text("")

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "trim_regions_single.py",
            "--ts", str(input_path),
            "--remove", str(mask_path),
            "--out", str(tmp_path / "out.trees"),
        ],
    )
    with pytest.raises(FileNotFoundError, match="No mutation rate available"):
        trs.main()
