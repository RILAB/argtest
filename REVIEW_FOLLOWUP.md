# Post-overhaul review follow-ups

Findings from reviewing the working-tree overhaul against `IMPLEMENTATION_PLAN.md`
on 2026-08-09. Test baseline at time of review: **89 passed, 1 failed, 5 skipped**.

**Status: cross-reviewed by codex (2026-08-09).** Issues 1, 2, 3, 4, 6 confirmed;
issue 5 was refuted and is withdrawn; the fix for issue 3 was corrected. Both
corrections are marked inline below. Issue 2 is confirmed **deliberate** — it
needs documentation and a pre-release path check, not a revert.

**Resolution (2026-08-09):**

| # | State | Notes |
|---|-------|-------|
| 1 | **FIXED** | `tests/conftest.py` added; per-file `sys.path` bootstraps removed from 15 modules. Every test file now imports in isolation. |
| 2 | **NO CODE CHANGE — deliberate** | Documented as breaking in `CHANGELOG.md`. Still needs the real-corpus `infer_mu_path` check below before tagging. |
| 3 | **FIXED** | `optional_mu_ratemap` returns `None` when `msprime is None`; corrupt-map errors still propagate. Regression test added. |
| 4 | **FIXED** | `require_pinned_tskit` now runs a behavioural probe instead of matching `".dev"`. |
| 5 | **WITHDRAWN** | Finding was wrong; current code is correct. |
| 6 | **OPEN — needs HPC** | Cannot be measured on this machine. |

Suite after the fixes: **90 passed, 1 failed, 5 skipped** — the single failure is
the environmental tskit-fork check described above, unchanged from baseline.

v1.10 was tagged on 2026-08-09 with two items still open. They are tracked in
`TODO.md` at the repo root; both require the pinned environment and real data:

1. Verify `infer_mu_path` resolves on the real amaranth/admix layout (issue 2).
   Note the risk is NOT a crash when `mutation_rate` is set — the scalar is a
   last-resort fallback, so a discovery miss silently substitutes a flat rate and
   drops step 3's local-rate correction. Run `scripts/check_mu_paths.py`
   **without** `--mutation-rate` so misses show as MISSING rather than `scalar`.
   Step 4 now also stamps `mu_source` into each output ARG and the pipeline
   summary reports it, so a completed run can be audited after the fact.
2. Benchmark merge peak RSS; revert to incremental concatenation if batching is
   worse (issue 6). Codex flagged that batching *can* raise peak RSS, since it
   holds every input tree sequence at once — treat this as a live possibility.

The single failure (`tests/test_coalescence_ne_plots.py::test_native_pair_quantiles_condition_out_isolated_pair_spans`)
is **environmental, not a regression** — the test file is unmodified and its body
is pure tskit. The machine has stock tskit 1.0.0 rather than the pinned
`nspope/tskit@73d8cd922482475020ae01180cae95bf5abbf067` from `environment.yml`.
Rebuild the env before treating it as a code defect.

---

## 1. `tests/test_validation_plots.py` fails when run in isolation

**Severity: medium** (breaks `-k`, `-x`, `pytest-xdist`, random ordering, and any
single-file run; currently passes only because another test module runs first and
mutates `sys.path`).

`tests/test_validation_plots.py:5` uses:

```python
from scripts import validation_plots_from_ts as validation
```

Every other test module instead does:

```python
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
```

There is no `scripts/__init__.py` and no `conftest.py`, so `scripts` resolves as a
namespace package whose modules cannot satisfy their own top-level
`from argtest_common import ...`.

Reproduce:

```
$ python -m pytest tests/test_validation_plots.py -q
E   ModuleNotFoundError: No module named 'argtest_common'
```

**Fix (pick one):**
- Add `tests/conftest.py` containing the `sys.path.insert` above — fixes this
  permanently for all present and future test modules, and lets every existing
  module drop its private copy of the incantation.
- Or make `tests/test_validation_plots.py` match the existing convention.

Preference: the `conftest.py` route, since it removes the duplicated boilerplate
rather than adding a ninth copy of it.

---

## 2. Step 4 gained a hard failure with no fallback

**Severity: high** (most likely change to break an existing working config).

`scripts/trim_regions_single.py` now calls `resolve_mu_rate(...)` unguarded, and
`resolve_mu_rate` raises `FileNotFoundError` when nothing is available.
Previously the rate map was optional:

```python
mu = ratemap_from_metadata(ts.metadata or {})
if mu is None:
    try: ...
    except Exception: mu = None
if mu is not None:
    metadata.update(ratemap_to_metadata(mu))
```

Two independent tightenings landed on this same path in one release:
1. `infer_mu_path` narrowed to exact candidates only (no suffix-derived bases, no
   glob).
2. Step 4 now requires the result.

Net: a dataset with no `mutation_rate:` in config, no embedded map, and no
*exact*-path sibling `*.mut_rate.p` fails at step 4 where v1.9 ran to completion.

**Action required — confirm intent.** If deliberate (the plan's Phase 4 implies
it is), no code change is needed; it is documented as breaking in the changelog.
If not, restore an opt-out. Either way, sanity-check the real amaranth/admix
layout resolves under the new exact-path rules **before** tagging:

```
python -c "
import sys; sys.path.insert(0,'scripts')
from argtest_common import infer_mu_path
from pathlib import Path
infer_mu_path(Path('<a real input treefile>'))
"
```

---

## 3. New uncaught exception path in validation plots

**Severity: low** (requires msprime to be absent, which `environment.yml` makes
unlikely — but the old code explicitly tolerated it).

`scripts/validation_plots_from_ts.py`:

```python
def optional_mu_ratemap(ts, ts_path: Path):
    try:
        return resolve_mu_rate(ts, ts_path)
    except FileNotFoundError:
        return None
```

`resolve_mu_rate` opens with `_require_msprime()`, which raises **`RuntimeError`**,
not `FileNotFoundError`. The module deliberately treats msprime as optional
(`except ImportError: msprime = None`), and the helper this replaced
(`_try_sibling_mu_ratemap`) swallowed every exception.

So without msprime, a run that is *not* `--sim-branch` and whose ARGs carry no
`kept_intervals` now crashes instead of falling back to whole-sequence
accessibility.

**Fix:** return `None` early when `msprime is None`, before calling
`resolve_mu_rate`.

Do **not** broaden the handler to `except (FileNotFoundError, RuntimeError)`.
That was this document's original suggestion and it is wrong:
`tests/test_validation_plots.py:38-43` asserts that a `RuntimeError` from
`resolve_mu_rate` (corrupt pickle, "Unrecognized mutation-rate object")
propagates out of `optional_mu_ratemap`. Catching it would break that test and
restore exactly the error-swallowing fallback this overhaul removed. Credit to
codex for catching this.

---

## 4. `require_pinned_tskit()` uses a weak version marker

**Severity: low.**

`scripts/coalescence_ne_plots_from_ts.py`:

```python
if missing or ".dev" not in tskit.__version__:
    raise RuntimeError(...)
```

`".dev" in __version__` is a proxy, not a capability check: any stock development
build of tskit passes it, and the guard breaks if the fork ever publishes a
non-`.dev` version. The `hasattr` half of the check is sound.

**Fix options:** pin on the known fork version string, or run the actual
partial-missing-data probe (the four-node/eight-tree fixture already in
`tests/test_coalescence_ne_plots.py::_partially_isolated_ts`) once at startup.
The probe is more honest and costs microseconds; the plan explicitly warned
against "a substantial synthetic calculation at every import," and this fixture
is not substantial.

---

## 5. ~~`merge_group` deletes a loop variable~~ — WITHDRAWN

**No action. This finding was wrong; the current code is correct as written.**

`scripts/merge_treefiles_by_replicate.py`:

```python
first_ts, *remaining = tree_sequences
merged = first_ts.concatenate(*remaining) if remaining else first_ts
del tree_sequences, remaining, first_ts, ts
```

This document originally proposed reducing that to `del tree_sequences`, on the
theory that the list held the last strong reference to each input. It does not.
Starred unpacking binds **independent** references to every element: `first_ts`
to element 0, and a fresh `remaining` list to elements 1..N-1. The loop variable
`ts` binds the last element as well. Verified:

```
after starred unpack, refcounts: [4, 4, 4, 4]
after `del tree_sequences` alone -> reachable via first_ts + remaining: [0, 1, 2, 3]
```

So `del tree_sequences` alone would free **nothing** and would silently defeat
the memory optimization the statement exists to provide. All four names must be
dropped. `del ts` is likewise necessary, not incidental — it is the only other
reference to the final input beyond `remaining`.

The one real (cosmetic) wrinkle is that `del ts` would raise `NameError` if the
loop never ran. That cannot happen: `grouped` is a `defaultdict(list)` populated
only via `.append`, so no group is ever empty. If that coupling is worth
removing, restructure metadata extraction so no loop-local reference survives —
do not shorten the `del`.

Caught by codex during review of this document.

---

## 6. Merge memory improvement is asserted but unmeasured

**Severity: process, not code.**

`merge_group` now carries the comment:

> This avoids materialising a progressively larger intermediate tree sequence at
> every chromosome boundary.

The plan set a target of ~34% lower peak RSS and ~15% lower wall time, explicitly
flagged as "a target to reproduce, not an assumed guarantee." That measurement has
not been taken on the pinned environment or on real data, and cannot be taken on
this machine (stock tskit, no real inputs).

Note the theoretical peak is similar either way — incremental holds
`merged_{N-1} + merged_N + ts_N`; batched holds all inputs plus the output. The
win, if real, comes from tskit's internal allocation behavior, which is exactly
why it needs measuring rather than reasoning.

**Action:** run `measure_merge_mem.py` (or the reproducible profiling script) on a
representative replicate under the pinned env before the changelog claims a number.
The changelog currently states the measurement is outstanding — keep it that way
until it is done.

---

## Also worth doing before tagging

- `scripts/validation_plots_from_ts.py`: `collect_stats`'s `mu_ratemap` parameter
  is now passed `None` at both call sites. Dead parameter; drop it or use it.
- Rebuild the conda env and re-run the full suite so the coalescence test is
  actually exercised, rather than failing for an unrelated reason.
- Run `scripts/audit_arg_contract.py` across the real amaranth/admix corpus. That
  audit is the gate on every deferred enforcement decision (uniform ploidy, leaf
  samples, VCF one-column-per-node fallback removal), and it has not been run on
  anything but fixtures.
