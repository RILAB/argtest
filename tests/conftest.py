"""Shared pytest configuration.

Puts the repository root and ``scripts/`` on ``sys.path`` so test modules can
import pipeline scripts either as top-level modules (``import argtest_common``)
or through the package path (``from scripts import validation_plots_from_ts``).

pytest loads this before collecting any test module, so it works for a single
file, a ``-k`` selection, and parallel or randomized ordering alike. Test
modules must NOT re-do this themselves: a per-file ``sys.path.insert`` only runs
once that particular module is imported, which silently makes correctness depend
on collection order.
"""
from __future__ import annotations

import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[1]

for _path in (_REPO_ROOT / "scripts", _REPO_ROOT):
    _entry = str(_path)
    if _entry not in sys.path:
        sys.path.insert(0, _entry)
