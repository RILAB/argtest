# Minimum-Retained-Samples Filter Plan

## Goal

Add a `--min-samples N` filter that prevents trees or intervals with too few
retained samples after pruning from distorting per-window statistics.

The motivating case is post-`trim_samples.py` tree sequences where some local
trees retain only a small subset of samples. Per-window statistics such as
mutational load, diversity, Tajima's D, SFS, and coalescence summaries can be
badly biased or unstable in those low-sample intervals.

## Decisions (locked 2026-06-26)

These supersede the exploratory options below where they conflict (notably the
"start with mask mode" recommendation):

- **Placement: new pipeline step after step 5** (Option B). A new script run
  between `trim_samples.py` (step 5) and `merge_treefiles_by_replicate.py`.
- **Mode: drop via `ts.delete_intervals(...)`** (not mask, and not coordinate-
  compacting `trim()`). `delete_intervals` removes the low-sample intervals'
  content but **preserves sequence coordinates** (the dropped spans become empty
  gaps), so there is no coordinate-compaction / provenance headache. This makes
  the "excise vs mask" question moot: we excise content but keep coordinates.
- **Unit: sample nodes (haploids), not individuals.** `min_samples` = minimum
  number of non-isolated sample *nodes* present in a local tree. On the current
  haploid `admix` data this equals individuals (1 node/individual). An optional
  `--min-samples-unit individuals` mode is deferred — do not build it now (it
  also entangles with the still-undecided per-individual semantics in todo #6).
- **`kept_intervals` coupling (required):** after dropping, intersect the
  existing `kept_intervals` metadata with the complement of the dropped spans,
  or downstream accessibility is overestimated (see `project_kept_intervals`
  memory and the `validation_plots_from_ts.py` accessibility path).
- **Difficulty: easy.** The building blocks already exist — non-isolated
  sample-per-tree counting is already done in `score_realistic_example.py`
  (~line 370), interval merging via `argtest_common.merge_intervals`, and the
  `keep/delete_intervals` + `kept_intervals` pattern in
  `trim_regions_single.py:38-40`. Core transform is ~20-40 lines plus Snakefile
  rule + a `min_samples` config knob + tests.

**Deferred:** agreed 2026-06-26 to save this refactor for another day.

## Core Behavior

For each local tree, compute the number of retained samples represented in that
tree after sample pruning. If the count is below `N`, treat the tree interval as
unusable for downstream statistics.

The main policy decision is whether unusable intervals are:

- physically excised from the tree sequence
- marked in metadata or a mask file and ignored by downstream stat functions

## Placement Options

### Option A: `trim_samples.py`

Add `--min-samples N` directly to `trim_samples.py`.

Pros:

- Local-sample loss is introduced by this step, so the filter is close to the
  cause.
- Downstream scripts receive cleaner tree sequences automatically.
- Easier to reason about pipeline outputs: step 5 becomes "sample pruning plus
  minimum retained-sample cleanup."

Cons:

- `trim_samples.py` already has a focused job: remove individuals globally or
  over BED intervals. Adding interval filtering makes it more complex.
- If intervals are excised, it changes coordinates again after sample pruning,
  which can complicate provenance.
- If intervals are only marked, downstream scripts still need to learn the mask.

Best if the desired behavior is a hard cleanup immediately after pruning.

### Option B: New Pipeline Step

Add a new script, for example `scripts/filter_min_samples.py`, and make it a
separate Snakemake step after `trim_samples.py` and before validation/merge.

Pros:

- Clean separation of concerns.
- Can produce explicit outputs: filtered tree sequence, BED of low-sample
  intervals, log/summary table.
- Easier to test independently.
- Leaves `trim_samples.py` simple and makes the filtering policy visible in the
  workflow.

Cons:

- Adds another pipeline step and output directory.
- More Snakefile/config plumbing.
- Existing users of `trim_samples.py` outside Snakemake need to opt into the new
  script separately.

Best default design if the feature should be auditable and configurable.

### Option C: `merge_treefiles_by_replicate.py`

Add `--min-samples N` to the merge step.

Pros:

- Filtering happens before per-replicate combined outputs are written.
- Avoids storing an intermediate filtered tree sequence per chromosome.

Cons:

- Too late for chromosome-level validation summaries.
- The merge script currently concatenates chromosome-specific tree sequences;
  local filtering is a different responsibility.
- Low-sample intervals could already have distorted earlier step outputs.
- Harder to report per-chromosome low-sample masks before concatenation.

This is not the preferred home unless the only desired consumer is the merged
treefile output.

## Recommended Design

Implement this as a new post-pruning step:

```text
trim_samples.py -> filter_min_samples.py -> merge_treefiles_by_replicate.py
```

The new script should support both output modes:

- `--mode mask`: write a BED mask of low-sample intervals and preserve the tree
  sequence unchanged
- `--mode excise`: remove low-sample intervals and compact coordinates

Default recommendation: start with `--mode mask`.

Masking is safer because it preserves coordinates and provenance. Once all
downstream stat scripts know how to consume the mask, `--mode excise` can be
added or enabled for workflows that prefer physically compacted outputs.

## Proposed CLI

```bash
python scripts/filter_min_samples.py \
  --ts input.tsz \
  --min-samples 20 \
  --mode mask \
  --out-mask low_sample.bed \
  --out-ts filtered.tsz \
  --log logs/filter_min_samples.log
```

Suggested flags:

- `--ts PATH`: input tree sequence
- `--min-samples INT`: minimum retained sample count per local tree
- `--mode {mask,excise}`: mark or remove low-sample intervals
- `--out-mask PATH`: BED of intervals with retained samples below threshold
- `--out-ts PATH`: output tree sequence for excise mode, or optional metadata
  annotated copy for mask mode
- `--log PATH`: summary log

For Snakemake, add config keys:

```yaml
min_samples: null
min_samples_mode: mask
```

If `min_samples` is `null`, skip the step.

## Retained-Sample Counting

The implementation needs a robust definition of "retained samples" per tree.

Candidate approach:

```python
for tree in ts.trees(sample_lists=True):
    retained = sum(1 for u in ts.samples() if tree.is_descendant(u, tree.root))
```

But this may be too slow and may not handle multiple roots cleanly. Better:

- use `tree.samples(root)` across all roots
- count unique sample nodes with at least one path in the local tree
- treat trees with no edges as zero retained samples

The helper should return merged half-open intervals:

```python
[(left, right, retained_count), ...]
```

Adjacent intervals with the same pass/fail status can be merged for BED output.

## Downstream Integration

In mask mode, downstream scripts should accept an optional low-sample mask and
exclude those windows or intervals from stat calculations. Candidate consumers:

- `validation_plots_from_ts.py`
- `mutload_masks.py`
- `mutload_summary.py`
- `coalescence_ne_plots_from_ts.py`
- `pipeline_summary.py`

In excise mode, downstream scripts can mostly operate on the filtered tree
sequence as they do now, but they should report that coordinates have been
compacted.

## Outputs

For each input tree sequence:

- BED file of low-sample intervals:

  ```text
  chrom start end retained_samples min_samples
  ```

- log file with:

  - input path
  - threshold
  - mode
  - total bp filtered
  - fraction of sequence filtered
  - min/median/max retained sample count across trees

- optional output tree sequence in excise mode or annotated mask mode

## Tests

Add tests for:

- counting retained samples on a small tree sequence with full and partial
  intervals
- merging adjacent low-sample intervals
- `--mode mask` writes the expected BED without changing the input tree sequence
- `--mode excise` removes the expected intervals and preserves valid tree
  sequence structure
- Snakemake behavior when `min_samples: null` versus an integer threshold

## Review Backlog (2026-07-04)

Efficiency follow-ups from the read-only review:

- `step1_low_rec_masks` runs once per chromosome but each job reparses the full
  HapMap before filtering to one chromosome. Stream only the requested
  chromosome in `hapmap_low_rec_mask.py`, or make the workflow emit all
  chromosome masks from one parse.
- `merge_replicates` invokes `merge_treefiles_by_replicate.py` once per
  replicate, but the script scans the whole filtered tree directory each time.
  Pass exact Snakemake inputs or use a replicate-specific pattern.
- `merge_group()` loads all chromosome tree sequences at once before
  concatenating. Reduce peak RAM by extracting metadata first and loading /
  concatenating iteratively or in a balanced tree.
- `argtest_common.mutational_load()` slices a sub-tree-sequence for every
  removal segment and expands descendant sample lists per mutation. Rework as a
  single pass over trees/sites with precomputed segment membership, or use
  tskit stats/sample-count primitives where possible.
- Step 3 computes observed mutational load and then simulation-based expected load
  through the same heavy traversal. Consider an analytic expectation, batched
  simulation, or shared traversal state.
- `filter_min_samples.py` has avoidable nested interval scans in
  `_intersect_keep_with_drops()` and dropped-BED annotation. Replace with
  two-pointer sweeps over sorted intervals/count segments.
- Validation and coalescence plotting scripts repeatedly load large tree
  sequences and keep many per-replicate arrays in memory. Stream summaries,
  cache intermediate metrics, or use memmaps where exact quantiles are needed.
- `pipeline_summary.py` reads BED-like files with `read_text().splitlines()`.
  Stream line-by-line for large masks/outlier files.

Test follow-ups from the read-only review:

- Add a regression that `trim_samples.py` drops mutations carried by targeted
  sample nodes inside trimmed intervals.
- Add coverage for `find_low_access_regions.py` with metadata ratemaps and with
  scalar `--mutation-rate` fallback.
- Add `mutload_masks.py --fraction` threshold tests for below, equal, and above
  the cutoff, since the code uses `>` rather than `>=`.
- Strengthen the real Snakemake integration test with content assertions:
  expected outlier names, mask contents, nonempty merged tree sequences, and
  evidence that downstream rules consumed upstream outputs.
- Extend VCF export tests for diploid grouping, mixed-ploidy fallback to
  per-node columns, missing individual metadata fallback, and name-substring removal.
- Add `combine_remove_masks.py` tests for missing inputs, comments/blank lines,
  malformed short lines, zero/negative intervals, and empty merged output.
- Replace tautological/implementation-detail assertions in `tests/test_mutload.py`
  with structural checks: sample IDs/order, expected sites/mutations, valid
  tree-sequence construction, and mutation-parent integrity. The ineffective
  `validate_trimmed_ts` compatibility probe is intentionally removed.

## Open Questions

- ~~Should the first implementation support both `mask` and `excise`, or only
  `mask`?~~ **Resolved 2026-06-26:** neither — drop via `delete_intervals`
  (excise content, preserve coordinates). See Decisions above.
- Should low-sample intervals be combined with existing remove masks, or kept as
  a separate diagnostic mask? *(still open — likely also emit a diagnostic BED
  alongside the dropped-interval `delete_intervals`.)*
- ~~Should `--min-samples` mean haploid sample nodes or diploid individuals?~~
  **Resolved 2026-06-26: sample nodes.** Individual mode deferred.
- Should windows be dropped when any overlapping tree is below threshold, or
  weighted by the fraction of the window that passes? *(still open — note the
  filter drops at tree-interval granularity, so window stats inherit gaps.)*
- Should the filter use absolute `N` only, or also support a fraction such as
  `--min-sample-fraction 0.8`? *(still open.)*

## Difficulty

Moderate if implemented as mask-only: one new script, a small shared helper,
Snakefile/config wiring, and focused tests.

Higher if implemented as excision immediately, because coordinate compaction,
metadata provenance, and downstream interpretation all need careful handling.
