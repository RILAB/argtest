import pickle
from pathlib import Path
from types import SimpleNamespace

import sys
import tskit

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import find_low_access_regions as finder
import trim_regions as tr


def make_simple_ts(length=10, site_pos=1):
    tables = tskit.TableCollection(sequence_length=length)
    pop = tables.populations.add_row()
    n0 = tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0, population=pop)
    n1 = tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0, population=pop)
    anc = tables.nodes.add_row(time=1, population=pop)
    tables.edges.add_row(left=0, right=length, parent=anc, child=n0)
    tables.edges.add_row(left=0, right=length, parent=anc, child=n1)
    site = tables.sites.add_row(position=site_pos, ancestral_state="0")
    tables.mutations.add_row(site=site, node=n0, derived_state="1")
    tables.sort()
    return tables.tree_sequence()


def test_find_low_access_regions_writes_bed(tmp_path, monkeypatch):
    ts = make_simple_ts(length=10, site_pos=6)
    ts_dir = tmp_path / "trees"
    ts_dir.mkdir()
    ts_path = ts_dir / "base.chr1.1.trees"
    ts.dump(ts_path)
    mu_path = ts_dir / "base.chr1.mut_rate.p"
    with open(mu_path, "wb") as fh:
        pickle.dump(SimpleNamespace(position=[0, 3, 10], rate=[1.0, 0.0]), fh)

    out_bed = tmp_path / "low_access.bed"
    monkeypatch.setattr(
        finder,
        "parse_args",
        lambda: type(
            "A",
            (),
            {
                "ts_dir": ts_dir,
                "window_size": 5.0,
                "cutoff_bp": 4.0,
                "pattern": "*",
                "out": out_bed,
            },
        )(),
    )
    finder.main()

    assert out_bed.exists()
    assert out_bed.read_text().strip().splitlines() == [
        "base.chr1.1\t0\t5\t3.000",
        "base.chr1.1\t5\t10\t0.000",
    ]


def test_trim_regions_applies_bed_mask(tmp_path, monkeypatch):
    ts = make_simple_ts(length=10, site_pos=5)
    ts_dir = tmp_path / "trees"
    ts_dir.mkdir()
    ts_path = ts_dir / "base.chr1.1.trees"
    ts.dump(ts_path)
    bed = tmp_path / "mask.bed"
    bed.write_text("chr1\t2\t4\n")

    monkeypatch.setattr(
        tr,
        "parse_args",
        lambda: type(
            "A",
            (),
            {
                "ts_dir": ts_dir,
                "remove": bed,
                "out_dir": tmp_path / "trimmed",
                "pattern": "*",
                "log": None,
            },
        )(),
    )
    tr.main()

    out = tmp_path / "trimmed" / "base.chr1.1.trimmed.mask.trees"
    assert out.exists()
    trimmed = tskit.load(out)
    assert trimmed.sequence_length == 8
