import pickle
from pathlib import Path

import msprime
import tskit


import find_low_access_regions as finder


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
    mu_path = ts_dir / "base.chr1.1.mut_rate.p"
    with open(mu_path, "wb") as fh:
        pickle.dump(msprime.RateMap(position=[0, 3, 10], rate=[1.0, 0.0]), fh)

    out_bed = tmp_path / "low_access.bed"
    monkeypatch.setattr(
        finder,
        "parse_args",
        lambda: type(
            "A",
            (),
            {
                "ts": ts_path,
                "chrom": "chr1",
                "window_size": 5.0,
                "cutoff_bp": 4.0,
                "out": out_bed,
                "log": None,
                "mutation_rate": None,
            },
        )(),
    )
    finder.main()

    assert out_bed.exists()
    assert out_bed.read_text().strip().splitlines() == [
        "chr1\t0\t5\t3.000",
        "chr1\t5\t10\t0.000",
    ]


def test_find_low_access_regions_creates_output_parent(tmp_path, monkeypatch):
    ts = make_simple_ts(length=10)
    ts_path = tmp_path / "replicate.trees"
    ts.dump(ts_path)
    out_bed = tmp_path / "new" / "nested" / "low_access.bed"

    monkeypatch.setattr(
        finder,
        "parse_args",
        lambda: type(
            "A",
            (),
            {
                "ts": ts_path,
                "chrom": "chr7",
                "window_size": 5.0,
                "cutoff_bp": 4.0,
                "out": out_bed,
                "log": None,
                "mutation_rate": 1.0,
            },
        )(),
    )
    finder.main()

    assert out_bed.exists()
