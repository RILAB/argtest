from pathlib import Path

import sys
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import combine_remove_masks as crm


def test_combine_remove_masks_merges_overlaps(tmp_path, monkeypatch):
    a = tmp_path / "a.bed"
    b = tmp_path / "b.bed"
    a.write_text("chr1\t0\t10\nchr1\t20\t30\n")
    b.write_text("chr1\t5\t15\nchr1\t30\t40\n")
    out = tmp_path / "out.bed"

    monkeypatch.setattr(
        crm,
        "parse_args",
        lambda: type(
            "A",
            (),
            {
                "chrom": "chr1",
                "out": out,
                "inputs": [a, b],
                "log": None,
            },
        )(),
    )
    crm.main()

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
