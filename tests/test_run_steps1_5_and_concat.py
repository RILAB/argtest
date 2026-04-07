from pathlib import Path

import sys
import tskit

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import run_steps1_5_and_concat as runner


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


def test_discover_chromosome_dirs(tmp_path):
    (tmp_path / "chr1").mkdir()
    (tmp_path / "chr2").mkdir()
    (tmp_path / "misc").mkdir()
    (tmp_path / "chr1" / "1.trees").touch()
    (tmp_path / "chr2" / "1.trees").touch()
    (tmp_path / "misc" / "note.txt").write_text("x")

    found = runner.discover_chromosome_dirs(tmp_path, "*")
    assert sorted(found.keys()) == ["chr1", "chr2"]
    assert found["chr1"][0].name == "1.trees"


def test_concatenate_by_replicate_writes_output(tmp_path):
    ts1 = make_simple_ts(length=5, site_pos=1)
    ts2 = make_simple_ts(length=7, site_pos=2)
    p1 = tmp_path / "chr1_1.trees"
    p2 = tmp_path / "chr2_1.trees"
    ts1.dump(p1)
    ts2.dump(p2)

    out_dir = tmp_path / "combined"
    written = runner.concatenate_by_replicate(
        replicate_map={"1": [("chr1", p1), ("chr2", p2)]},
        out_dir=out_dir,
        base_name="base",
        out_suffix=None,
        n_chromosomes=2,
        allow_missing=False,
    )

    assert len(written) == 1
    assert written[0].name == "base.combined.1.trees"
    merged = tskit.load(written[0])
    assert merged.sequence_length == 12


def test_concatenate_by_replicate_raises_on_missing_chromosome(tmp_path):
    ts1 = make_simple_ts(length=5, site_pos=1)
    p1 = tmp_path / "chr1_1.trees"
    ts1.dump(p1)

    try:
        runner.concatenate_by_replicate(
            replicate_map={"1": [("chr1", p1)]},
            out_dir=tmp_path / "combined",
            base_name="base",
            out_suffix=None,
            n_chromosomes=2,
            allow_missing=False,
        )
    except RuntimeError as exc:
        assert "missing" in str(exc)
    else:
        raise AssertionError("Expected RuntimeError for missing chromosome replicate")
