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
    grouped, skipped, _ = merger.group_tree_files(paths)
    assert skipped == []
    assert sorted(grouped.keys()) == [("base", "1"), ("base", "2")]


def make_simple_ts_with_kept(length, site_pos, kept_intervals):
    ts = make_simple_ts(length=length, site_pos=site_pos)
    tables = ts.dump_tables()
    tables.metadata_schema = tskit.MetadataSchema({"codec": "json"})
    tables.metadata = {"kept_intervals": [[float(l), float(r)] for l, r in kept_intervals]}
    return tables.tree_sequence()


def test_merge_group_offsets_kept_intervals(tmp_path):
    ts1 = make_simple_ts_with_kept(length=5, site_pos=1, kept_intervals=[[0, 3]])
    ts2 = make_simple_ts_with_kept(length=7, site_pos=2, kept_intervals=[[1, 4], [5, 7]])
    p1 = tmp_path / "base.chr1.1.trees"
    p2 = tmp_path / "base.chr2.1.trees"
    ts1.dump(p1)
    ts2.dump(p2)

    merged, _ = merger.merge_group([("chr1", p1), ("chr2", p2)])

    assert merged.sequence_length == 12
    assert merged.metadata["kept_intervals"] == [[0.0, 3.0], [6.0, 9.0], [10.0, 12.0]]


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
                "layout": "flat",
                "base_name": None,
                "out_suffix": None,
                "replicate": None,
            },
        )(),
    )
    merger.main()

    out = tmp_path / "combined" / "base.combined.1.trees"
    assert out.exists()
    merged = tskit.load(out)
    assert merged.sequence_length == 12


def test_group_tree_files_nested_layout(tmp_path):
    p1 = tmp_path / "chr1" / "1.trees"
    p2 = tmp_path / "chr2" / "1.trees"
    p3 = tmp_path / "chr1" / "2.trees"
    p1.parent.mkdir()
    p2.parent.mkdir()
    p1.touch()
    p2.touch()
    p3.touch()

    grouped, skipped, detected = merger.group_tree_files(
        [p1, p2, p3],
        ts_dir=tmp_path,
        layout="nested",
        base_name="base",
    )

    assert skipped == []
    assert detected == "nested"
    assert sorted(grouped.keys()) == [("base", "1"), ("base", "2")]
    assert sorted(chrom for chrom, _ in grouped[("base", "1")]) == ["chr1", "chr2"]


def test_main_writes_combined_file_nested_layout(tmp_path, monkeypatch):
    ts1 = make_simple_ts(length=5, site_pos=1)
    ts2 = make_simple_ts(length=7, site_pos=2)
    in_dir = tmp_path / "trees"
    (in_dir / "chr1").mkdir(parents=True)
    (in_dir / "chr2").mkdir(parents=True)
    ts1.dump(in_dir / "chr1" / "1.trees")
    ts2.dump(in_dir / "chr2" / "1.trees")

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
                "layout": "nested",
                "base_name": "base",
                "out_suffix": None,
                "replicate": None,
            },
        )(),
    )
    merger.main()

    out = tmp_path / "combined" / "base.combined.1.trees"
    assert out.exists()
    merged = tskit.load(out)
    assert merged.sequence_length == 12


def test_main_replicate_filter_writes_only_matching(tmp_path, monkeypatch):
    ts1 = make_simple_ts(length=5, site_pos=1)
    ts2 = make_simple_ts(length=7, site_pos=2)
    in_dir = tmp_path / "trees"
    in_dir.mkdir()
    ts1.dump(in_dir / "base.chr1.1.trees")
    ts2.dump(in_dir / "base.chr2.1.trees")
    ts1.dump(in_dir / "base.chr1.2.trees")
    ts2.dump(in_dir / "base.chr2.2.trees")

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
                "layout": "flat",
                "base_name": None,
                "out_suffix": None,
                "replicate": "1",
            },
        )(),
    )
    merger.main()

    out_dir = tmp_path / "combined"
    assert (out_dir / "base.combined.1.trees").exists()
    assert not (out_dir / "base.combined.2.trees").exists()
