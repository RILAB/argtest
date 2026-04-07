from pathlib import Path

import sys
import tskit

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import merge_treefiles_by_replicate as merger


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


def test_parse_input_name():
    parsed = merger.parse_input_name(Path("maize.chr10.1.tsz"))
    assert parsed == ("maize", "chr10", "1")
    assert merger.parse_input_name(Path("badname.tsz")) is None


def test_group_tree_files():
    paths = [
        Path("base.chr1.1.trees"),
        Path("base.chr2.1.trees"),
        Path("base.chr1.2.trees"),
    ]
    grouped, skipped = merger.group_tree_files(paths)
    assert skipped == []
    assert sorted(grouped.keys()) == [("base", "1"), ("base", "2")]


def test_merge_group_concatenates_sequence_length(tmp_path, monkeypatch):
    ts1 = make_simple_ts(length=5, site_pos=1)
    ts2 = make_simple_ts(length=7, site_pos=2)
    p1 = tmp_path / "base.chr2.1.trees"
    p2 = tmp_path / "base.chr10.1.trees"
    ts1.dump(p1)
    ts2.dump(p2)

    merged, ordered = merger.merge_group([("chr10", p2), ("chr2", p1)])

    assert [chrom for chrom, _ in ordered] == ["chr2", "chr10"]
    assert merged.sequence_length == 12
    assert merged.num_sites == 2


def test_main_writes_combined_file(tmp_path, monkeypatch):
    ts1 = make_simple_ts(length=5, site_pos=1)
    ts2 = make_simple_ts(length=7, site_pos=2)
    in_dir = tmp_path / "trees"
    in_dir.mkdir()
    ts1.dump(in_dir / "base.chr1.1.trees")
    ts2.dump(in_dir / "base.chr2.1.trees")

    monkeypatch.setattr(
        merger,
        "parse_args",
        lambda: type(
            "A",
            (),
            {
                "ts_dir": in_dir,
                "out_dir": tmp_path / "combined",
                "pattern": "*",
                "out_suffix": None,
            },
        )(),
    )
    merger.main()

    out = tmp_path / "combined" / "base.combined.1.trees"
    assert out.exists()
    merged = tskit.load(out)
    assert merged.sequence_length == 12
