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


def test_mutload_masks_writes_outlier_and_empty_masked(tmp_path, monkeypatch):
    ts = make_simple_ts()
    ts_path = tmp_path / "1.trees"
    ts.dump(ts_path)
    outlier = tmp_path / "outliers.bed"
    masked = tmp_path / "masked.bed"

    monkeypatch.setattr(
        mm,
        "parse_args",
        lambda: type(
            "A",
            (),
            {
                "ts": ts_path,
                "chrom": "chr1",
                "window_size": 10.0,
                "snp_window": None,
                "cutoff": 0.25,
                "fraction": None,
                "suffix_to_strip": "_anchorwave",
                "outlier_bed": outlier,
                "masked_bed": masked,
                "log": None,
            },
        )(),
    )
    mm.main()

    assert outlier.exists()
    assert outlier.read_text().strip() != ""
    assert masked.exists()
    assert masked.read_text().strip() == ""
    assert (tmp_path / "logs" / "chr1.1.mutload.log").exists()
