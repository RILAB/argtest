from pathlib import Path


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


# --------------------------------------------------------------------------- #
# Chromosome-name resolution
#
# Pipeline chroms ("combined.1"), bare-numeric maps ("1") and "chr"-prefixed
# maps ("chr1") all coexist in real datasets. Resolution must work in BOTH
# directions: a mismatch silently attaches the wrong recombination map or the
# wrong sequence length to a chromosome.
# --------------------------------------------------------------------------- #

import pytest

from scripts import hapmap_low_rec_mask as hlrm


@pytest.mark.parametrize("query, available, expected", [
    # Exact match always wins.
    ("chr1", {"chr1", "chr2"}, "chr1"),
    ("1", {"1", "2"}, "1"),
    # Bare query against a chr-prefixed file (was already supported).
    ("1", {"chr1"}, "chr1"),
    ("1", {"chr_1"}, "chr_1"),
    ("1", {"chr-1"}, "chr-1"),
    # chr-prefixed query against a bare file (the direction that used to fail).
    ("chr1", {"1"}, "1"),
    ("chr_1", {"1"}, "1"),
    ("chr-1", {"1"}, "1"),
    # Case is not significant either way.
    ("Chr1", {"1"}, "1"),
    ("1", {"CHR1"}, "CHR1"),
    # Pipeline base-name prefixes are stripped from the query, and the result
    # still resolves across the chr convention in both directions.
    ("combined.1", {"1"}, "1"),
    ("combined.1", {"chr1"}, "chr1"),
    ("amaranth.16", {"chr16"}, "chr16"),
    ("combined.chr1", {"1"}, "1"),
    ("combined.chr1", {"chr1"}, "chr1"),
])
def test_chromosome_resolution_is_symmetric(query, available, expected):
    assert hlrm._resolve_chrom(query, available) == expected


@pytest.mark.parametrize("query, available", [
    ("chr3", {"1", "2"}),
    ("combined.3", {"chr1", "chr2"}),
    ("scaffold_9", {"chr1"}),
])
def test_unresolvable_chromosomes_return_none(query, available):
    assert hlrm._resolve_chrom(query, available) is None


def test_unrelated_contigs_are_not_collapsed_onto_a_chromosome():
    """Only the query gets its base-name prefix stripped, never the available names.

    Stripping 'chrUn_random.1' down to '1' would make it collide with 'chr1',
    so a query of 'combined.1' would become ambiguous on a perfectly ordinary
    reference FAI.
    """
    available = {"chr1", "chrUn_random.1"}
    assert hlrm._resolve_chrom("combined.1", available) == "chr1"
    assert hlrm._resolve_chrom("1", available) == "chr1"


def test_two_spellings_of_one_chromosome_raise_instead_of_guessing():
    """Which one you got would otherwise depend on set iteration order."""
    with pytest.raises(ValueError, match="matches more than one entry"):
        hlrm._resolve_chrom("combined.1", {"1", "chr1"})


def test_ambiguity_error_names_the_offending_file():
    with pytest.raises(ValueError, match="in /data/ref.fa.fai"):
        hlrm._resolve_chrom("combined.1", {"1", "chr1"}, source="/data/ref.fa.fai")


def test_exact_match_wins_over_an_otherwise_ambiguous_set():
    """A file using one convention consistently is never reinterpreted."""
    assert hlrm._resolve_chrom("chr1", {"1", "chr1"}) == "chr1"
    assert hlrm._resolve_chrom("1", {"1", "chr1"}) == "1"
