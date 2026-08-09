from pathlib import Path

import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import hapmap_low_rec_mask as hlrm


def test_rec_fraction_zero_writes_empty_low_rec_mask(tmp_path, monkeypatch):
    hapmap = tmp_path / "hapmap.tsv"
    hapmap.write_text(
        "Chromosome\tPosition(bp)\tRate(cM/Mb)\tMap(cM)\n"
        "chr1\t0\t1.0\t0.0\n"
        "chr1\t100\t0.1\t0.1\n"
        "chr1\t200\t2.0\t0.2\n"
    )
    fai = tmp_path / "genome.fai"
    fai.write_text("chr1\t300\t0\t0\t0\n")
    out_dir = tmp_path / "out"

    monkeypatch.setattr(
        hlrm,
        "parse_args",
        lambda: type(
            "A",
            (),
            {
                "hapmap": hapmap,
                "fai": fai,
                "rec_fraction": 0.0,
                "out_dir": out_dir,
                "chrom": None,
                "log": None,
            },
        )(),
    )

    hlrm.main()

    bed = out_dir / "chr1.low_rec.mask.bed"
    assert bed.exists()
    assert bed.read_text() == ""
    assert "Total intervals written: 0" in (
        out_dir / "logs" / "hapmap_low_rec_mask.log"
    ).read_text()


def test_chromosome_filter_scans_hapmap_once_and_preserves_requested_alias(
    tmp_path, monkeypatch
):
    hapmap = tmp_path / "hapmap.tsv"
    hapmap.write_text(
        "Chromosome\tPosition(bp)\tRate(cM/Mb)\tMap(cM)\n"
        "chr1\t0\t1.0\t0.0\n"
        "chr1\t100\t0.1\t0.1\n"
        "chr2\t0\t2.0\t0.0\n"
    )
    fai = tmp_path / "genome.fai"
    fai.write_text("chr1\t200\t0\t0\t0\nchr2\t200\t0\t0\t0\n")
    out_dir = tmp_path / "out"
    real_load = hlrm.load_hapmap
    calls = []

    def recording_load(*args, **kwargs):
        calls.append((args, kwargs))
        return real_load(*args, **kwargs)

    monkeypatch.setattr(hlrm, "load_hapmap", recording_load)
    monkeypatch.setattr(
        hlrm,
        "parse_args",
        lambda: type(
            "A",
            (),
            {
                "hapmap": hapmap,
                "fai": fai,
                "rec_fraction": 1.0,
                "out_dir": out_dir,
                "chrom": "combined.1",
                "log": None,
            },
        )(),
    )

    hlrm.main()

    assert len(calls) == 1
    bed = out_dir / "combined.1.low_rec.mask.bed"
    assert bed.exists()
    assert {line.split("\t")[0] for line in bed.read_text().splitlines()} == {
        "combined.1"
    }
