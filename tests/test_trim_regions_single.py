from pathlib import Path

import sys
import tskit

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import trim_regions_single as trs


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
