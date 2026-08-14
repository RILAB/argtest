from pathlib import Path

import pytest

from scripts import validation_plots_from_ts as validation


def test_explicit_tree_files_preserve_caller_order():
    args = validation.parse_args(["--ts", "rep10.tsz", "rep2.tsz"])

    assert args.ts == [Path("rep10.tsz"), Path("rep2.tsz")]
    assert args.ts_dir is None


def test_tree_directory_and_explicit_files_are_mutually_exclusive():
    with pytest.raises(SystemExit):
        validation.parse_args(["--ts", "rep1.tsz", "--ts-dir", "trees"])


def test_explicit_comparison_files_preserve_caller_order():
    args = validation.parse_args(
        ["--ts", "primary.tsz", "--compare-ts", "clean10.tsz", "clean2.tsz"]
    )

    assert args.compare_ts == [Path("clean10.tsz"), Path("clean2.tsz")]
    assert args.compare is None


def test_optional_ratemap_only_swallows_absence(monkeypatch):
    def missing(*args, **kwargs):
        raise FileNotFoundError("not present")

    monkeypatch.setattr(validation, "resolve_mu_rate", missing)
    assert validation.optional_mu_ratemap(object(), Path("rep1.tsz")) is None


def test_optional_ratemap_propagates_corrupt_map_error(monkeypatch):
    def corrupt(*args, **kwargs):
        raise RuntimeError("bad pickle")

    monkeypatch.setattr(validation, "resolve_mu_rate", corrupt)
    with pytest.raises(RuntimeError, match="bad pickle"):
        validation.optional_mu_ratemap(object(), Path("rep1.tsz"))


def test_optional_ratemap_falls_back_when_msprime_is_absent(monkeypatch):
    """msprime is optional in this module, but resolve_mu_rate requires it.

    Without the guard, `_require_msprime` raises RuntimeError before the helper
    can report absence, so a non-sim-branch run over ARGs with no kept_intervals
    would crash instead of falling back to whole-sequence accessibility.
    """
    def unreachable(*args, **kwargs):  # pragma: no cover - must not be called
        raise AssertionError("resolve_mu_rate must not run without msprime")

    monkeypatch.setattr(validation, "msprime", None)
    monkeypatch.setattr(validation, "resolve_mu_rate", unreachable)
    assert validation.optional_mu_ratemap(object(), Path("rep1.tsz")) is None
