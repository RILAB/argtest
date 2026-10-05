from pathlib import Path

import pytest


import combine_remove_masks as crm


def _run_main(monkeypatch, out, inputs):
    monkeypatch.setattr(
        crm,
        "parse_args",
        lambda: type(
            "A",
            (),
            {
                "chrom": "chr1",
                "out": out,
                "inputs": inputs,
                "log": None,
            },
        )(),
    )
    crm.main()


def test_combine_remove_masks_merges_overlaps(tmp_path, monkeypatch):
    a = tmp_path / "a.bed"
    b = tmp_path / "b.bed"
    a.write_text("chr1\t0\t10\nchr1\t20\t30\n")
    b.write_text("chr1\t5\t15\nchr1\t30\t40\n")
    out = tmp_path / "out.bed"

    _run_main(monkeypatch, out, [a, b])

    lines = [line.strip() for line in out.read_text().splitlines() if line.strip()]
    assert lines == ["chr1\t0\t15", "chr1\t20\t40"]
    assert (tmp_path / "logs" / "chr1.combine_masks.log").exists()


def test_read_intervals_rejects_missing_input(tmp_path):
    missing = tmp_path / "missing.bed"
    with pytest.raises(FileNotFoundError, match=str(missing)):
        crm.read_intervals(missing)


def test_read_intervals_rejects_malformed_row_with_context(tmp_path):
    bed = tmp_path / "bad.bed"
    bed.write_text("# header\nchr1 0\n")
    with pytest.raises(ValueError, match=r"bad\.bed at line 2"):
        crm.read_intervals(bed)


def test_read_intervals_accepts_empty_file(tmp_path):
    bed = tmp_path / "empty.bed"
    bed.write_text("")
    assert crm.read_intervals(bed) == []


def test_rounds_outward_before_merging(tmp_path):
    bed = tmp_path / "float.bed"
    bed.write_text("chr1\t0.1\t1.1\nchr1\t1.9\t2.1\n")
    # The rounded intervals [0, 2] and [1, 3] overlap and must merge.
    assert crm.merge_intervals(crm.read_intervals(bed)) == [[0, 3]]


def test_read_intervals_skips_comments_and_blank_lines(tmp_path):
    bed = tmp_path / "commented.bed"
    bed.write_text("# header\n\nchr1\t0\t10\n   \n# mid comment\nchr1\t20\t30\n")
    assert crm.read_intervals(bed) == [[0, 10], [20, 30]]


def test_read_intervals_drops_zero_and_negative_with_warning(tmp_path, capsys):
    bed = tmp_path / "degenerate.bed"
    bed.write_text("chr1\t5\t5\nchr1\t0\t10\nchr1\t30\t20\n")
    dropped = []
    assert crm.read_intervals(bed, dropped=dropped) == [[0, 10]]
    assert dropped == [(bed, 1, "5", "5"), (bed, 3, "30", "20")]
    err = capsys.readouterr().err
    assert "zero-length interval" in err and "line 1" in err
    assert "negative-length interval" in err and "line 3" in err


def test_dropped_intervals_recorded_in_log(tmp_path, monkeypatch):
    bed = tmp_path / "degenerate.bed"
    bed.write_text("chr1\t5\t5\nchr1\t0\t10\nchr1\t30\t20\n")
    out = tmp_path / "out.bed"
    _run_main(monkeypatch, out, [bed])

    assert out.read_text() == "chr1\t0\t10\n"
    log = (tmp_path / "logs" / "chr1.combine_masks.log").read_text()
    assert "dropped_intervals=2\n" in log
    assert f"dropped\t{bed}:1\tstart=5\tend=5\n" in log
    assert f"dropped\t{bed}:3\tstart=30\tend=20\n" in log


def test_empty_merged_output_writes_empty_bed(tmp_path, monkeypatch):
    empty = tmp_path / "empty.bed"
    empty.write_text("")
    comments_only = tmp_path / "comments.bed"
    comments_only.write_text("# nothing here\n\n")
    degenerate = tmp_path / "degenerate.bed"
    degenerate.write_text("chr1\t7\t7\n")
    out = tmp_path / "out.bed"
    _run_main(monkeypatch, out, [empty, comments_only, degenerate])

    assert out.exists()
    assert out.read_text() == ""
    log = (tmp_path / "logs" / "chr1.combine_masks.log").read_text()
    assert "merged_intervals=0\n" in log
    assert "dropped_intervals=1\n" in log
