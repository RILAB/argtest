# argtest code review — 2026-08-09

Reviewed at commit `6691acb` (v1.9). Scope: `Snakefile`, `scripts/*.py`, `tests/*.py`,
`README.md`, `CHANGELOG.md`, `measure_merge_mem.py`. ~9,700 lines of Python.

**No code was changed.** Every claim marked *(measured)* was verified by running code;
everything else is from reading.

---

## Executive summary

The pipeline is in good shape: the hot paths (`trim_samples`, `filter_min_samples`)
have already been vectorized thoughtfully, the Snakefile config validation is sound,
and test coverage is real. The problems are concentrated in three places:

1. **~700 lines (≈7% of the Python) are dead** — superseded implementations kept alive
   only by tests, plus one guard function that provably does nothing but print a
   spurious warning on every run.
2. **The merge step's memory profile can be cut ~34% with a one-line change** *(measured)*.
   This is the step CLAUDE.md documents as the OOM risk.
3. **Documentation has drifted** — the README's workflow section still points at the
   superseded `trim_regions.py` for step 4, which has *different semantics* than the
   script the pipeline actually runs.

Priority ranking is at the end.

---

## A. Bugs and correctness

### A1. `validate_trimmed_ts` is a no-op that prints a false warning on every run — HIGH

[argtest_common.py:275-303](scripts/argtest_common.py#L275-L303)

Called by [trim_samples.py:361](scripts/trim_samples.py#L361) and
[filter_min_samples.py:263](scripts/filter_min_samples.py#L263) — i.e. every step-5 and
step-5b job. Traced against the installed tskit *(measured)*:

- `hasattr(ts, "check_index")` → `True`, so the outer branch is entered.
- `ts.check_index()` raises `TypeError` (real signature is `check_index(index, length)`).
- In the `except TypeError` handler, `hasattr(ts, "index")` → `False` and
  `hasattr(ts.tables, "index")` → `False`, so `index` stays `None`.
- It prints `WARNING: trimmed ts skipped check_index (no index available)` and **returns**.
- `ts.validate()` is never reached — and `TreeSequence` has no `validate` method anyway.

So 29 lines of triple-nested defensive branching validate nothing, and every step-5 log
in the repo contains a warning that reads like a real problem. Recommend: delete the
function and its two call sites. If validation is genuinely wanted, `tables.tree_sequence()`
(already called upstream in both paths) raises on a malformed table collection, which is
the actual guarantee you have today.

### A2. `find_low_access_regions.py` writes the replicate ID into the BED chrom column

[find_low_access_regions.py:105](scripts/find_low_access_regions.py#L105) — `chrom = ts_path.stem`.

In the nested layout the Snakefile feeds it (`<root>/<chrom>/<rep>.tsz` via
[Snakefile:427-434](Snakefile#L427-L434)), the stem is the **replicate ID**, not the
chromosome. Verified by running it on `example_data/sim_2chr_5rep/trees/chr1/1.trees`
*(measured)*:

```
1	0	10000	10000.000     ← should be "chr1"
```

Harmless inside the pipeline because `combine_remove_masks.py` rewrites column 1 from
`--chrom {wildcards.chrom}`, but the step-2 artifact under `step2_low_access/<chrom>/`
is mislabeled for anyone reading it directly or feeding it to bedtools. Fix: add a
`--chrom` argument (matching every other script in the pipeline) and pass
`{wildcards.chrom}`.

### A3. `trim_regions_single.py` discards pre-existing tree-sequence metadata

[trim_regions_single.py:46-58](scripts/trim_regions_single.py#L46-L58)

```python
metadata: dict = {"kept_intervals": [...]}
...
tables.metadata = metadata          # wholesale replacement
```

Any top-level metadata the input ARG carried (e.g. SINGER provenance fields) is dropped
at step 4. Compare [merge_treefiles_by_replicate.py:245-247](scripts/merge_treefiles_by_replicate.py#L245-L247)
and [filter_min_samples.py:246-249](scripts/filter_min_samples.py#L246-L249), which both
do `{**existing, **extra}`. Fix: merge instead of replace, for consistency.

### A4. `combine_remove_masks.py` truncates float mask boundaries inward

[combine_remove_masks.py:65](scripts/combine_remove_masks.py#L65) —
`f"{args.chrom}\t{int(left)}\t{int(right)}"`.

For a *removal* mask, truncating `right` loses up to 1 bp of masked sequence per interval.
Should be `floor(left)` / `ceil(right)` to stay conservative. In practice step-1/2/3
boundaries are integral so this rarely fires, but SNP-windowed mutload masks
([mutload_masks.py:164-168](scripts/mutload_masks.py#L164-L168)) emit float window edges
derived from site positions, and those *do* truncate — and can even produce a zero-width
BED row when two consecutive window edges land inside the same integer.

### A5. Inconsistent BED-parsing strictness

- [argtest_common.py:320](scripts/argtest_common.py#L320) — raises `ValueError` on <3 columns.
- [trim_regions.py:71](scripts/trim_regions.py#L71) — raises `ValueError`.
- [combine_remove_masks.py:46-47](scripts/combine_remove_masks.py#L46-L47) — **silently skips**.

The README ([line 322](README.md#L322)) documents the raising behavior as universal. A
malformed step-3 BED reaching step 4 would be silently ignored rather than failing the
job. Pick one (recommend: raise everywhere, it's one shared helper) and delete the other two.

### A6. Pinned-tskit dependency fails open, not closed

`environment.yml` pins `tskit @ git+nspope/tskit@73d8cd9` because
[coalescence_ne_plots_from_ts.py](scripts/coalescence_ne_plots_from_ts.py) depends on that
fork's partial-missing-data normalization (documented at [README.md:304](README.md#L304)).

Running the test suite against **upstream** tskit 1.0.0 *(measured)*:

```
FAILED tests/test_coalescence_ne_plots.py::
  test_native_pair_quantiles_condition_out_isolated_pair_spans
  assert np.isfinite(result).all()
  → array([1., nan, nan, nan, 8.])
```
*(78 passed, 5 skipped — snakemake absent. The `argtest` conda env on this machine is
not provisioned, so this ran in base anaconda; treat the failure as "wrong tskit
installed", not necessarily a code defect.)*

The point stands regardless: with the wrong tskit the script produces NaN rates and Ne
values **without erroring**, and `validate_trimmed_ts` (A1) silently degrades to a
warning. A one-line import-time assertion on the required tskit behavior/version would
turn a silent-wrong-answer into a startup failure. Cheap insurance for a pinned fork.

### A7. Minor

- [filter_min_samples.py:171](scripts/filter_min_samples.py#L171) — `_intersect_keep_with_drops(keep, dropped, seq_len)`: `seq_len` is never used. Dead parameter.
- [trim_regions_single.py:50-54](scripts/trim_regions_single.py#L50-L54) — bare `except Exception: mu = None` around the pickle load swallows corrupt-file errors. Also, if the pickle holds a *scalar* rate, `ratemap_to_metadata(mu)` on line 56 throws outside the try.
- [measure_merge_mem.py:38](measure_merge_mem.py#L38) — `ru_maxrss / 1e6` assumes kB (Linux). On macOS `ru_maxrss` is bytes, so the reported peak is 1000× off. Fine on the HPC target, worth a comment.
- [Snakefile:202](Snakefile#L202) — `TRIM_REMOVE_ARG` builds shell args by string interpolation; a BED path containing `"` breaks the command. Same for `--chroms {params.chroms}` at [Snakefile:813](Snakefile#L813).

---

## B. Over-engineering — what can be deleted

This is the main ask. Everything below is **dead in production**, verified by grep across
`scripts/`, `Snakefile`, and `tests/`.

### B1. Dead-code inventory

| What | Where | Lines | Only referenced by |
|---|---|---:|---|
| `trim_samples_chunked.py` (whole module) | [scripts/trim_samples_chunked.py](scripts/trim_samples_chunked.py) | 146 | nothing at all |
| `trim_ts_by_intervals` | [argtest_common.py:196-258](scripts/argtest_common.py#L196-L258) | 63 | `test_mutload.py` |
| `validate_trimmed_ts` | [argtest_common.py:275-303](scripts/argtest_common.py#L275-L303) | 29 | see A1 — a no-op |
| `collapse_masked_intervals` | [argtest_common.py:631-659](scripts/argtest_common.py#L631-L659) | 29 | `trim_regions.py` only |
| `ratemap_from_keep_intervals` | [argtest_common.py:593-618](scripts/argtest_common.py#L593-L618) | 26 | `trim_regions.py` only |
| `build_shared_mask` | [argtest_common.py:674-697](scripts/argtest_common.py#L674-L697) | 24 | `collapse_masked_and_low_access_windows` only |
| `mutational_load` removal branch | [argtest_common.py:119-159](scripts/argtest_common.py#L119-L159) | ~25 | nothing — no caller ever passes `remove_intervals` |
| `build_removal_segments` | [argtest_common.py:162-181](scripts/argtest_common.py#L162-L181) | 20 | the dead branch above + tests |
| `assert_sample_ids_preserved` | [argtest_common.py:261-272](scripts/argtest_common.py#L261-L272) | 12 | nothing |
| `build_segments_with_drop_nodes` | [argtest_common.py:184-193](scripts/argtest_common.py#L184-L193) | 10 | `trim_ts_by_intervals` only |
| `collapse_masked_and_low_access_windows` | [argtest_common.py:662-671](scripts/argtest_common.py#L662-L671) | 10 | nothing |
| `and_ratemaps_binary` | [argtest_common.py:621-628](scripts/argtest_common.py#L621-L628) | 8 | `build_shared_mask` only |
| `trim_regions.py` CLI (`main`, `find_tree_files`, `output_name`, `parse_args`) | [scripts/trim_regions.py](scripts/trim_regions.py) | ~110 | nothing — superseded by `trim_regions_single.py` |

**≈ 510 lines of `scripts/`, plus ~150 lines of `test_mutload.py`** that exist solely to
exercise the dead `trim_ts_by_intervals` / `build_removal_segments` path
(`test_trim_preserves_coordinates`, `test_trim_removes_nodes_in_interval`,
`test_trim_mutation_parent_integrity`, `test_all_samples_removed_in_segment`,
`test_idempotent_no_effect_remove`, `test_segment_merge_logic`, …).

Notes on the two entries that need care:

- **`trim_regions.py`** — keep `complement_intervals` and `load_mask_intervals`
  (27 lines; `trim_regions_single.py` imports them at
  [line 18](scripts/trim_regions_single.py#L18)). Move them into `argtest_common.py`
  and delete the rest of the module. Deleting the module also removes the
  `collapse_masked_intervals` / `ratemap_from_keep_intervals` chain above — those two
  implement *coordinate-compacting* trimming, which the pipeline deliberately does not
  do any more (it preserves coordinates via `keep_intervals`).
- **`trim_samples_chunked.py`** — its docstring cites
  `reports/profiling/verify_trim_chunked.py` for the equivalence proof; that file does
  not exist (the closest is `verify_trim_samples.py`). If the chunked path is a
  hedge against a future slow `simplify()`, say so in a one-line comment and keep it;
  otherwise it's 146 lines of untested, unreferenced duplicate edge-surgery logic that
  will silently drift from `trim_samples.py`.

### B2. Guards that can never fire

- [mutload_masks.py:149-150](scripts/mutload_masks.py#L149-L150) —
  `if expected_names != unique_names: raise`. Both lists come from
  `aggregate_by_individual(..., names)` with the *same* `names` object, so they are
  always equal by construction.
- [pipeline_summary.py:413-415](scripts/pipeline_summary.py#L413-L415) and
  [:437-444](scripts/pipeline_summary.py#L437-L444) — `if row_by_rep is None` and the
  `key == "retained_vals" and "retained_by_rep" in row` / positional-index fallback.
  `collect_retention` always populates `retained_by_rep`
  ([line 347](scripts/pipeline_summary.py#L347)), so the fallback branches are
  unreachable. `totals_by_replicate` collapses to about six lines without them.
- [find_low_access_regions.py:96](scripts/find_low_access_regions.py#L96) and
  [:120](scripts/find_low_access_regions.py#L120) — `getattr(args, "mutation_rate", None)`
  / `getattr(args, "log", None)`. `argparse` guarantees both attributes exist. Just use
  `args.mutation_rate` / `args.log`.

### B3. `infer_mu_path` — the fallback is worse than no fallback

[argtest_common.py:519-566](scripts/argtest_common.py#L519-L566)

`infer_mu_base` generates candidate stems by stripping `.N`, `_N`, **and** `-N` suffixes
from both the file stem and the parent directory name, then `infer_mu_path` does an exact
lookup, and if that fails, a **glob** across two directories filtered by `startswith`.

CLAUDE.md already documents the failure mode this creates: for `amaranth/amaranth.16.tsz`
with a missing `amaranth.16.mut_rate.p`, the broad base `"amaranth"` matches all 15 other
chromosomes' files and the error reads *"Ambiguous mutation maps"* instead of *"file not
found"*. That is the fallback actively making diagnosis harder — a documented,
already-cost-you-time bug caused by defensive code.

Recommendation: keep the two exact lookups (`<ts_stem>.mut_rate.p`, `<parent_name>.mut_rate.p`
in `ts.parent` and `ts.parent.parent`), delete `infer_mu_base`'s suffix-stripping and the
glob/ambiguity path entirely, and raise `FileNotFoundError` listing the exact paths tried.
~35 lines → ~12, and the error message becomes actionable.

### B4. `resolve_mu_rate` carries a test-fixture shim in production

[argtest_common.py:455-458](scripts/argtest_common.py#L455-L458)

```python
# Some fixtures pickle a SimpleNamespace with .position/.rate; rebuild a RateMap.
if hasattr(obj, "position") and hasattr(obj, "rate"):
    return msprime.RateMap(position=list(obj.position), rate=list(obj.rate))
```

The comment says the plural quiet part out loud. Fix the fixtures to pickle real
`RateMap`s and delete the branch — otherwise a genuinely wrong pickled object gets
duck-typed into a silently wrong ratemap.

### B5. Snakefile step 6 stages files through `/tmp` with 35 lines of bash

[Snakefile:722-756](Snakefile#L722-L756)

Two `mktemp -d`, a trap, two symlink loops, and an inline Python heredoc that imports
`argtest_common` just to resolve a `mut_rate.p` path — all because
`validation_plots_from_ts.py` takes `--ts-dir` + `--pattern` instead of a list of files.

Giving that script a `--ts FILE [FILE ...]` argument (which is how every other script in
the pipeline is invoked) collapses the whole rule body to two `python` calls. It also
removes a hard dependency on `/tmp` having room, which is a real hazard on diskless
compute nodes.

### B6. Duplicated helpers

- [compare_trees_html.py:8-21](scripts/compare_trees_html.py#L8-L21) and
  [trees_gallery_html.py:8-21](scripts/trees_gallery_html.py#L8-L21) each re-implement
  `load_ts` verbatim instead of importing it from `argtest_common`.
- [score_realistic_example.py:64-74](scripts/score_realistic_example.py#L64-L74)
  re-implements `merge_intervals`.
- [mutload_masks.py:101-115](scripts/mutload_masks.py#L101-L115) and
  [mutload_summary.py:140-158](scripts/mutload_summary.py#L140-L158) have near-identical
  `build_bp_windows` / `build_snp_windows` (the only difference is where the `<= 0`
  validation lives).
- `genome_windows`-style "arange then clamp the last edge" appears **five** times:
  [argtest_common.py:684-686](scripts/argtest_common.py#L684-L686),
  [mutload_masks.py:103-106](scripts/mutload_masks.py#L103-L106),
  [find_low_access_regions.py:98-100](scripts/find_low_access_regions.py#L98-L100),
  [validation_plots_from_ts.py:200-203](scripts/validation_plots_from_ts.py#L200-L203),
  [coalescence_ne_plots_from_ts.py:489-491](scripts/coalescence_ne_plots_from_ts.py#L489-L491).
  One shared helper.

### B7. `validation_plots_from_ts.py` `main()` is 480 lines of near-identical plot blocks

[validation_plots_from_ts.py:535-1011](scripts/validation_plots_from_ts.py#L535-L1011)

Diversity / segregating-sites / Tajima's D each get an identical *scatter → skyline →
trace* triple, differing only in the dict keys and axis labels. That's ~330 lines that a
single loop over a spec table

```python
STATS = [("div", "Diversity", ...), ("s", "Segregating sites / base", ...), ("td", "Tajima's D", ...)]
```

reduces to roughly 80. Not urgent, but this is the file most likely to grow a
copy-paste bug (each block already has slightly different NaN handling — compare
`site_td = np.where(np.isnan(div_scale), ...)` at
[line 390](scripts/validation_plots_from_ts.py#L390) against
`sim_td = np.where(rep_acc <= 0, ...)` at
[line 432](scripts/validation_plots_from_ts.py#L432)).

### B8. Things I looked at and would **keep**

Flagging these so they don't get swept up:

- `safe_log_yscale` ([validation_plots_from_ts.py:228](scripts/validation_plots_from_ts.py#L228)) —
  guards a real matplotlib crash-at-savefig on empty spectra. Cheap, keep.
- The Snakefile `resources:` validation block ([Snakefile:65-88](Snakefile#L65-L88)) —
  catches config typos at parse time rather than after a 6-hour job. Keep.
- `export_vcf.resolve_samples` mixed-ploidy fallback ([export_vcf.py:53-80](scripts/export_vcf.py#L53-L80)) —
  small, and the haploid-vs-diploid distinction is genuinely load-bearing for the admix data.
- `hapmap_low_rec_mask._resolve_chrom` chr-prefix matching ([hapmap_low_rec_mask.py:119-133](scripts/hapmap_low_rec_mask.py#L119-L133)) —
  cross-naming-convention hapmaps are a real recurring input problem (CLAUDE.md documents it). Keep.
- `filter_min_samples._sample_edge_coverage` — the sweep is more code than a per-tree
  loop but is asymptotically better and correctly reasoned. Keep.
- `mutload_seed_for` ([Snakefile:146-151](Snakefile#L146-L151)) — deterministic per-(chrom, rep) seeds. Keep.

---

## C. Efficiency

### C1. Batch the merge concatenate — 34% lower peak RSS, 15% faster *(measured)*

[merge_treefiles_by_replicate.py:210-219](scripts/merge_treefiles_by_replicate.py#L210-L219)

```python
merged = first_ts
for chrom, path in chrom_paths[1:]:
    ts = load_ts(path)
    ...
    merged = merged.concatenate(ts)      # rebuilds the growing table N-1 times
```

`tskit.TreeSequence.concatenate` accepts varargs — `first.concatenate(*rest)` — and does
it in one pass. Benchmarked on 16 synthetic chromosomes (10 Mb, 100 diploids,
~76k edges each, 1.2M edges merged), loading from disk exactly as `merge_group` does:

| | wall time | peak RSS |
|---|---:|---:|
| incremental (current) | 4.27 s / 4.33 s | **1.48 GB / 1.46 GB** |
| batch `concatenate(*rest)` | 3.67 s / 3.69 s | **0.98 GB / 0.96 GB** |

Identical output (1,199,806 edges, matching edge tables). Batch wins on memory *even
though* it holds all inputs resident, because the incremental loop keeps both the old and
new merged table alive at every step — and that intermediate is what dominates by the
last few chromosomes.

Given CLAUDE.md documents `merge_replicates` as the OOM step and recommends budgeting
96–128 GB, a ~⅓ reduction from a one-line change is the highest-value item in this review.

### C2. Free the intermediate before the metadata rebuild

[merge_treefiles_by_replicate.py:243-248](scripts/merge_treefiles_by_replicate.py#L243-L248)

```python
tables = merged.dump_tables()      # full copy #2
...
merged = tables.tree_sequence()    # full copy #3
```

`extra` now *always* contains `chrom_offsets` ([line 229](scripts/merge_treefiles_by_replicate.py#L229)),
so this rebuild always runs and momentarily holds three whole-genome copies. Adding
`del merged` immediately after `dump_tables()` (and reading `existing` before it) drops
one copy at the exact peak. `measure_merge_mem.py` already isolates this stage — worth
re-running with the change.

### C3. `mutational_load` copies the whole tree sequence for nothing

[argtest_common.py:131](scripts/argtest_common.py#L131)

```python
sub = ts.keep_intervals([(left, right)], simplify=False)
```

In the only code path that exists (no `remove_intervals`, so `segments == [(0, L, ∅)]`),
this is `keep_intervals` over the *entire* sequence — a full table copy and re-index that
produces a tree sequence identical to the input. Called twice per `mutload_masks` job
(observed + simulated) and twice per replicate in `validation_plots_from_ts`. Skipping it
when `left == 0 and right == ts.sequence_length` is free.

### C4. `mutational_load`'s Python loop is ~2.7× slower than a `variants()` sweep *(measured)*

[argtest_common.py:145-158](scripts/argtest_common.py#L145-L158)

On 200 samples / 5 Mb / 10,481 trees / 11,832 sites: **0.16 s** for the current
tree→site→mutation→`tree.samples()` loop vs **0.06 s** for

```python
wi = np.digitize(ts.sites_position, windows) - 1
for var in ts.variants():
    out[wi[var.site.id]] += (var.genotypes != 0)
```

Note `sample_lists=True` on [line 145](scripts/argtest_common.py#L145) also makes tree
construction itself more expensive.

**Caveat — not a drop-in.** The two disagree by 523 of 399,794 counts (0.13%) on that
test, entirely from 17 sites with recurrent/back mutations: the current code counts a
sample once *per mutation* it sits below, the genotype version counts *derived state
present*. Both are defensible definitions of "mutational load"; the genotype one is
arguably the more standard. If you take the speedup, state the semantic change in the
CHANGELOG — the outlier cutoff in step 3 is calibrated against a simulated expectation
computed by the *same* function, so the change is self-consistent and shifts thresholds
only at the 0.1% level.

### C5. `coalescence_ne_plots_from_ts.py` loads every tree sequence twice

[compute_quantile_time_windows:237-243](scripts/coalescence_ne_plots_from_ts.py#L237-L243)
and [compute_logspaced_time_windows:297-303](scripts/coalescence_ne_plots_from_ts.py#L297-L303)
each `load_ts()` all post-burnin replicates to build the time grid, then
[main:584-597](scripts/coalescence_ne_plots_from_ts.py#L584-L597) loads all of them again
for the statistics. For whole-genome merged ARGs that doubles the I/O and decompression
of the most expensive input in the repo. Restructure to load once and pass the objects
(or accept the grid-building cost only for the first replicate).

### C6. `hapmap_low_rec_mask.py` reads the HapMap twice per invocation

[hapmap_low_rec_mask.py:143-146](scripts/hapmap_low_rec_mask.py#L143-L146) — `hapmap_chroms()`
scans the whole file to resolve the chromosome name, then `load_hapmap()` scans it again.
Since `step1_low_rec_masks` is a per-chromosome rule, a 16-chromosome run makes **32 full
passes** over the map. `load_hapmap` could return the chrom set it already saw, or the
rule could become a single job producing all chromosomes (which is what the script does
natively when `--chrom` is omitted).

### C7. `plt.clf()` leaks figures

`validation_plots_from_ts.py` calls `plt.subplots()` ~15 times and `plt.clf()` after each.
`clf()` clears the current figure but does **not** deregister it from pyplot, so figures
accumulate and matplotlib emits its `More than 20 figures have been opened` RuntimeWarning
once the `--compare` / `--sim-sfs` paths add more. Use `plt.close(fig)`.
Same pattern in `coalescence_ne_plots_from_ts.py`.

### C8. Micro

- [argtest_common.py:77-95](scripts/argtest_common.py#L77-L95) `aggregate_by_individual` —
  Python loop over samples; `np.add.at(agg, idx_array, load)` with a precomputed index
  array is a one-liner and vectorized.
- [argtest_common.py:342-354](scripts/argtest_common.py#L342-L354) `merge_intervals` —
  converts to numpy, sorts, then builds a Python list. Called **once per target node**
  in [trim_samples.py:265-267](scripts/trim_samples.py#L265-L267), so with thousands of
  trimmed nodes you pay thousands of tiny array round-trips. A pure-Python sort path for
  small inputs would be faster.
- [coalescence_ne_plots_from_ts.py:595](scripts/coalescence_ne_plots_from_ts.py#L595) —
  `if i in keep_post` on a numpy array is O(n) per iteration; use a set or `i >= burnin`.

---

## D. Documentation

### D1. README points step 4 at the wrong script — and the two differ semantically

[README.md:79](README.md#L79) (Suggested Workflow, which claims "The Snakemake pipeline
automates exactly these steps"):

> **Remove affected regions.** Combine the step 1–3 BED masks per chromosome and trim
> those regions from each tree sequence with `scripts/trim_regions.py`.

The pipeline runs [`trim_regions_single.py`](scripts/trim_regions_single.py)
([Snakefile:539](Snakefile#L539)). This is not just a naming nit:

| | `trim_regions.py` | `trim_regions_single.py` (actual) |
|---|---|---|
| mechanism | `collapse_masked_intervals` | `ts.keep_intervals(..., simplify=False)` |
| coordinates | **compacted** to unmasked length | **preserved** |

Every downstream coordinate assumption in the repo (`chrom_offsets`, `kept_intervals`,
`locate_tree.py`, the step-5b `delete_intervals` design) depends on coordinates being
preserved. `reports/upstream_stats_comparison.md:28,404` already flags this discrepancy —
it just never made it back into the README. Same issue at
[README.md:285](README.md#L285), [:319](README.md#L319), [:322](README.md#L322).

**If you take the B1 recommendation and delete `trim_regions.py`, all four references
must go.**

### D2. Stale / redundant docs

- Three completed design docs sit at repo root: `min_samples_filter_plan.md` (316 lines),
  `per_sample_trim_plan.md` (80), `sim_based_mutload_plan.md` (95). All three describe
  work that shipped (v1.7–v1.9 per the CHANGELOG). Move to `reports/` or delete — as
  root-level files they read like open work.
- `NOTES.md` (30 lines) is referenced only from
  [README.md:287](README.md#L287) for sample-ID matching rules. That content belongs in
  the README section that links to it.
- [trim_samples_chunked.py:13](scripts/trim_samples_chunked.py#L13) cites
  `reports/profiling/verify_trim_chunked.py`, which does not exist.
- [mutload_summary.py:126-131](scripts/mutload_summary.py#L126-L131) — `--out` docs say
  "Only the filename part is used; the file is always written to `<repo-root>/results/`".
  That's surprising behavior for a `--out` flag (it silently ignores any directory the
  user passes, and writes into the *repo*, not the cwd). Either honor the path or rename
  the flag to `--out-name`.

### D3. `pipeline_summary` and `validation_plots` measure "accessible" differently

- [pipeline_summary.py:316-321](scripts/pipeline_summary.py#L316-L321) derives accessible
  intervals from the **mutation ratemap** (`rate > 0`), and falls back to `None`
  (= whole sequence) when no ratemap is embedded.
- [validation_plots_from_ts.py:354-360](scripts/validation_plots_from_ts.py#L354-L360)
  prefers **`kept_intervals`** metadata, falling back to the ratemap, then to `None`.

Both are labeled "accessible bp" in their outputs. They usually agree (trimmed regions
become edge-less trees either way), but a chromosome whose metadata lost its ratemap gets
a silently inflated retention percentage in the HTML summary while the validation plots
use the correct mask. Worth either unifying on `kept_intervals`-first or documenting the
distinction in the README's retention paragraph ([README.md:83](README.md#L83)).

---

## Priority

**Do first (small, high value):**

1. **C1** — batch `concatenate(*rest)` in `merge_group`. One line, 34% peak RSS *(measured)*. Directly targets the documented OOM step.
2. **A1** — delete `validate_trimmed_ts` and its two call sites. Removes a permanent false warning from every step-5 log.
3. **D1** — fix the four README references to `trim_regions.py`. Actively misleading about coordinate semantics.
4. **C2** — `del merged` before the metadata rebuild. One line, one fewer whole-genome copy at peak.

**Then (the simplification pass you asked about):**

5. **B1** — delete the dead-code table. ~510 lines from `scripts/`, ~150 from `test_mutload.py`. Do `trim_samples_chunked.py` and `trim_regions.py`'s CLI first (largest, cleanest), then the `argtest_common` chain.
6. **B3** — cut `infer_mu_path`'s glob fallback. This one is not just cleanup — the fallback currently converts "file missing" into "ambiguous", which has already cost debugging time.
7. **C3** — skip the whole-sequence `keep_intervals` in `mutational_load`. Free.
8. **B2, B4, B6, A7** — unreachable guards, the fixture shim, duplicated helpers, dead params.

**Worth doing, larger:**

9. **B5** — give `validation_plots_from_ts.py` a `--ts FILE...` argument; delete 35 lines of Snakefile bash and the `/tmp` dependency.
10. **C4** — the `variants()` rewrite of `mutational_load`, *if* you accept the 0.13% recurrent-mutation semantic change (document it).
11. **B7** — table-drive the plot blocks in `validation_plots_from_ts.py`.
12. **A6** — assert the pinned-tskit behavior at import so a wrong install fails loudly.

**Lower priority:** A2–A5, C5–C8, D2, D3.

---

## Appendix: how the measurements were made

- **Merge benchmark (C1):** 16 tree sequences written to disk (`msprime.sim_ancestry`,
  100 diploids, `population_size=1e4`, 10 Mb, `recombination_rate=1e-8`, plus
  `sim_mutations` at 1e-8), then two subprocesses — one replicating `merge_group`'s
  incremental loop (holding `first_ts`, loading each subsequent chromosome inside the
  loop), one loading all and calling `concatenate(*rest)`. Peak from
  `resource.getrusage(RUSAGE_SELF).ru_maxrss`; two runs each, results within 2%.
- **`validate_trimmed_ts` (A1):** called directly on a simulated tree sequence; traced
  each branch with `hasattr` / `inspect.signature` against the installed tskit.
- **`mutational_load` (C3, C4):** 200 samples, 5 Mb, 10,481 trees, 11,832 sites, 50 windows;
  `time.time()` around each implementation, outputs compared with `np.allclose`; the
  discrepancy attributed to recurrent mutations by counting sites with >1 mutation (17)
  and mutations with `edge == NULL` (0).
- **BED chrom column (A2):** ran `find_low_access_regions.py` on
  `example_data/sim_2chr_5rep/trees/chr1/1.trees` with a cutoff forcing all windows to be
  flagged, and read column 1 of the output.
- **Test suite (A6):** `python -m pytest tests -q` in base anaconda (tskit 1.0.0 upstream);
  the project's `argtest` conda env exists but has no `msprime`/`pytest` installed on this
  machine.
- **Dead-code inventory (B1):** `grep -rn` for each symbol across `scripts/`, `Snakefile`,
  `tests/`, and `*.md`, excluding `.git/`.
