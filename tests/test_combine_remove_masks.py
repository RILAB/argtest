from pathlib import Path

import sys

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
