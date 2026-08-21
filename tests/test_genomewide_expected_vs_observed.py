"""Tests for the genome-wide expected-vs-observed scatter and its input dumps."""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import msprime
import numpy as np
import pytest

from scripts import genomewide_expected_vs_observed as gw
from scripts import validation_plots_from_ts as validation

REPO_ROOT = Path(__file__).resolve().parents[1]


# --------------------------------------------------------------------------- #
# Per-window recombination rate
# --------------------------------------------------------------------------- #

def test_mean_rate_is_length_weighted_not_a_plain_interval_mean():
    """A window spanning a long slow interval and a short fast one is slow.

    Averaging the intervals themselves would report (1 + 11) / 2 = 6; weighting
    by the bp each contributes gives (90*1 + 10*11) / 100 = 2.
    """
    intervals = [(0, 90, 1.0), (90, 100, 11.0)]
    rates = gw.mean_rate_per_window(
        np.array([0.0]), np.array([100.0]), intervals
    )
    assert rates[0] == pytest.approx(2.0)


def test_windows_outside_the_map_are_nan_not_zero():
    """"No map here" must not be plotted as "no recombination here"."""
    intervals = [(0, 100, 5.0)]
    rates = gw.mean_rate_per_window(
        np.array([200.0]), np.array([300.0]), intervals
    )
    assert np.isnan(rates[0])


def test_partial_overlap_uses_only_the_overlapping_bp():
    intervals = [(0, 50, 2.0), (50, 200, 8.0)]
    rates = gw.mean_rate_per_window(
        np.array([40.0]), np.array([60.0]), intervals
    )
    # 10 bp at rate 2 and 10 bp at rate 8.
    assert rates[0] == pytest.approx(5.0)


def test_empty_map_gives_all_nan():
    rates = gw.mean_rate_per_window(np.array([0.0, 10.0]), np.array([10.0, 20.0]), [])
    assert np.isnan(rates).all()


# --------------------------------------------------------------------------- #
# Colour normalisation
# --------------------------------------------------------------------------- #

def test_wide_rate_spread_uses_a_log_scale_with_zeros_clipped():
    rates = np.array([0.0, 0.01, 0.1, 1.0, 10.0])
    colors, norm, label = gw.rate_norm(rates)

    assert "log scale" in label
    assert "clipped" in label
    # The zero is lifted to the smallest positive rate so it still gets a colour.
    assert colors[0] == pytest.approx(0.01)
    assert norm.vmin == pytest.approx(0.01)
    assert norm.vmax == pytest.approx(10.0)


def test_narrow_rate_spread_stays_linear():
    rates = np.array([1.0, 1.5, 2.0])
    _colors, norm, label = gw.rate_norm(rates)

    assert "log scale" not in label
    assert norm.vmin == pytest.approx(1.0)
    assert norm.vmax == pytest.approx(2.0)


def test_all_nan_rates_do_not_crash_the_norm():
    _colors, norm, _label = gw.rate_norm(np.array([np.nan, np.nan]))
    assert norm.vmax > norm.vmin


# --------------------------------------------------------------------------- #
# windows.tsv parsing
# --------------------------------------------------------------------------- #

def _write_windows_tsv(path: Path, rows: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    header = "\t".join(
        ["dataset", "window_start", "window_end", "window_mid",
         "obs_pi", "exp_pi", "exp_pi_lo", "exp_pi_hi",
         "obs_tajimas_d", "exp_tajimas_d", "exp_tajimas_d_lo", "exp_tajimas_d_hi",
         "obs_segsites", "exp_segsites", "exp_segsites_lo", "exp_segsites_hi"]
    )
    path.write_text("\n".join([header] + rows) + "\n")


def test_only_the_first_dataset_is_pooled(tmp_path, capsys):
    """A --compare run puts two datasets in one file; mixing them would be wrong."""
    path = tmp_path / "chr1" / "cleaned" / "windows.tsv"
    _write_windows_tsv(path, [
        "primary\t0\t10\t5\t0.1\t0.2\t\t\t-1\t-2\t\t\t\t\t\t",
        "compare\t0\t10\t5\t0.9\t0.8\t\t\t1\t2\t\t\t\t\t\t",
    ])

    label, rows = gw.read_windows_tsv(path)

    assert label == "primary"
    assert len(rows) == 1
    assert "compare" in capsys.readouterr().err


def test_blank_cells_read_back_as_nan(tmp_path):
    path = tmp_path / "chr1" / "cleaned" / "windows.tsv"
    _write_windows_tsv(path, ["primary\t0\t10\t5\t0.1\t\t\t\t-1\t\t\t\t\t\t\t"])

    _label, rows = gw.read_windows_tsv(path)

    assert gw._cell(rows[0], "obs_pi") == pytest.approx(0.1)
    assert np.isnan(gw._cell(rows[0], "exp_pi"))


def test_chromosome_is_taken_from_the_grandparent_directory(tmp_path):
    """Matches the step-6 layout <step6_dir>/<chrom>/<variant>/windows.tsv."""
    path = tmp_path / "combined.7" / "cleaned" / "windows.tsv"
    _write_windows_tsv(path, ["primary\t0\t100\t50\t0.1\t0.2\t\t\t-1\t-2\t\t\t\t\t\t"])

    pooled = gw.collect([path], {"7": [(0, 1.0)]}, {"7": 100})

    assert [e["chrom"] for e in pooled] == ["combined.7"]


def test_chromosome_absent_from_the_map_is_dropped_with_a_warning(tmp_path, capsys):
    path = tmp_path / "combined.9" / "cleaned" / "windows.tsv"
    _write_windows_tsv(path, ["primary\t0\t100\t50\t0.1\t0.2\t\t\t-1\t-2\t\t\t\t\t\t"])

    pooled = gw.collect([path], {"1": [(0, 1.0)]}, {"1": 100})

    assert pooled == []
    assert "combined.9" in capsys.readouterr().err


def test_no_poolable_windows_is_an_error_not_an_empty_plot(tmp_path, monkeypatch):
    """Silently writing a blank scatter would read as "the ARGs look fine"."""
    path = tmp_path / "combined.9" / "cleaned" / "windows.tsv"
    _write_windows_tsv(path, ["primary\t0\t100\t50\t0.1\t0.2\t\t\t-1\t-2\t\t\t\t\t\t"])
    hapmap = tmp_path / "map.tsv"
    hapmap.write_text("Chromosome\tPosition\tRate\tMap\n1\t0\t1.0\t0.0\n")
    fai = tmp_path / "ref.fai"
    fai.write_text("1\t100\n")
    monkeypatch.setattr(sys, "argv", [
        "genomewide_expected_vs_observed.py",
        "--windows", str(path), "--hapmap", str(hapmap),
        "--fai", str(fai), "--out-dir", str(tmp_path / "out"),
    ])

    with pytest.raises(RuntimeError, match="No windows could be pooled"):
        gw.main()


# --------------------------------------------------------------------------- #
# End-to-end: validation dumps feed the genome-wide plot
# --------------------------------------------------------------------------- #

def _make_ts(path: Path, seed: int, length: int = 20_000):
    ts = msprime.sim_ancestry(
        samples=6, sequence_length=length, recombination_rate=1e-6,
        population_size=1000, random_seed=seed,
    )
    ts = msprime.sim_mutations(ts, rate=1e-6, random_seed=seed + 500)
    path.parent.mkdir(parents=True, exist_ok=True)
    ts.dump(str(path))


def test_validation_dumps_then_genomewide_plot(tmp_path):
    """The two features are one chain: step 6 writes the TSV, step 6b reads it."""
    step6 = tmp_path / "step6_validation"
    chroms = ["combined.1", "combined.2"]
    for i, chrom in enumerate(chroms):
        ts_dir = tmp_path / "trees" / chrom
        for rep in range(2):
            _make_ts(ts_dir / f"{rep}.trees", seed=1 + i * 10 + rep)
        out_dir = step6 / chrom / "cleaned"
        subprocess.run(
            [sys.executable, str(REPO_ROOT / "scripts" / "validation_plots_from_ts.py"),
             "--ts", *[str(p) for p in sorted(ts_dir.glob("*.trees"))],
             "--burnin-frac", "0", "--window-size", "5000",
             "--mutation-rate", "1e-6", "--sim-branch",
             "--out-dir", str(out_dir)],
            check=True, cwd=REPO_ROOT,
        )

    # Request 1: the data behind each plot axis is on disk and parseable.
    for chrom in chroms:
        col = step6 / chrom / "cleaned"
        for name in ("windows.tsv", "samples.tsv", "sfs.tsv"):
            assert (col / name).exists(), f"{chrom}/{name} was not written"
        header = (col / "windows.tsv").read_text().splitlines()[0].split("\t")
        assert "obs_pi" in header and "exp_pi" in header
        assert "obs_tajimas_d" in header and "exp_tajimas_d" in header
        assert "window_start" in header and "window_end" in header

    hapmap = tmp_path / "map.tsv"
    hapmap.write_text(
        "Chromosome\tPosition\tRate\tMap\n"
        # A slow first half and a fast second half on each chromosome, so the
        # colour axis has real spread to normalise over.
        "1\t0\t0.05\t0.0\n1\t10000\t5.0\t0.5\n"
        "2\t0\t0.05\t0.0\n2\t10000\t5.0\t0.5\n"
    )
    fai = tmp_path / "ref.fai"
    fai.write_text("1\t20000\n2\t20000\n")

    out_dir = step6 / "genomewide" / "cleaned"
    subprocess.run(
        [sys.executable,
         str(REPO_ROOT / "scripts" / "genomewide_expected_vs_observed.py"),
         "--windows", *[str(step6 / c / "cleaned" / "windows.tsv") for c in chroms],
         "--hapmap", str(hapmap), "--fai", str(fai), "--out-dir", str(out_dir)],
        check=True, cwd=REPO_ROOT,
    )

    # Request 2: both plots exist, and the pooled table backing them covers
    # every chromosome with a real rate attached.
    assert (out_dir / "diversity-expected-vs-observed-by-rec.png").exists()
    assert (out_dir / "tajima-d-expected-vs-observed-by-rec.png").exists()

    pooled = (out_dir / "genomewide-windows.tsv").read_text().splitlines()
    header = pooled[0].split("\t")
    rows = [dict(zip(header, line.split("\t"))) for line in pooled[1:]]
    assert {r["chrom"] for r in rows} == set(chroms)
    rates = np.array([float(r["rec_rate_cm_per_mb"]) for r in rows])
    assert np.isfinite(rates).all()
    # The map's slow and fast halves must both be represented.
    assert rates.min() < 1.0 < rates.max()
    assert all(r["obs_pi"] for r in rows)


# --------------------------------------------------------------------------- #
# Genome coverage labelling
#
# validation_first_chrom_only defaults to true, so step 6b routinely pools ONE
# chromosome. Every one of these guards the claim that the result is genome-wide.
# --------------------------------------------------------------------------- #

def test_full_coverage_is_labelled_genome_wide():
    note = gw.coverage_note({"chr1", "chr2"}, ["chr1", "chr2"])
    assert "genome-wide" in note
    assert "PARTIAL" not in note


def test_partial_coverage_is_labelled_partial_and_counts_both_sides():
    note = gw.coverage_note({"chr1"}, ["chr1", "chr2", "chr3"])
    assert note.startswith("PARTIAL GENOME")
    assert "1 of 3" in note
    assert "chr1" in note


def test_many_pooled_chromosomes_are_truncated_not_dumped():
    covered = {f"chr{i}" for i in range(1, 8)}
    note = gw.coverage_note(covered, [f"chr{i}" for i in range(1, 20)])
    assert "7 of 19" in note
    assert note.endswith("...)")


def test_coverage_is_not_claimed_when_the_full_list_is_unknown():
    """Without --all-chroms the script cannot know what it is missing."""
    note = gw.coverage_note({"chr1"}, None)
    assert "genome-wide" not in note
    assert "PARTIAL" not in note
    assert "1 chromosome pooled" == note


def test_partial_run_warns_on_stderr_and_titles_the_plot(tmp_path, monkeypatch, capsys):
    """A one-chromosome pool must not produce an unlabelled 'genome-wide' plot."""
    path = tmp_path / "chr1" / "cleaned" / "windows.tsv"
    _write_windows_tsv(path, [
        "primary\t0\t100\t50\t0.1\t0.2\t\t\t-1\t-2\t\t\t\t\t\t",
        "primary\t100\t200\t150\t0.3\t0.25\t\t\t-0.5\t-1\t\t\t\t\t\t",
    ])
    # Map and FAI are keyed by the chromosome directory name itself; unlike
    # "combined.1", a bare "chr1" has no numeric suffix to strip.
    hapmap = tmp_path / "map.tsv"
    hapmap.write_text("Chromosome\tPosition\tRate\tMap\nchr1\t0\t1.0\t0.0\n")
    fai = tmp_path / "ref.fai"
    fai.write_text("chr1\t200\n")
    out_dir = tmp_path / "out"
    monkeypatch.setattr(sys, "argv", [
        "genomewide_expected_vs_observed.py",
        "--windows", str(path), "--hapmap", str(hapmap), "--fai", str(fai),
        "--all-chroms", "chr1", "chr2", "chr3",
        "--out-dir", str(out_dir),
    ])

    gw.main()

    err = capsys.readouterr().err
    assert "NOT genome-wide" in err
    assert "chr2" in err and "chr3" in err
    assert (out_dir / "diversity-expected-vs-observed-by-rec.png").exists()
