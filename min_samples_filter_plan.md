# Minimum-Retained-Samples Filter Plan

## Goal

Add a `--min-samples N` filter that prevents trees or intervals with too few
retained samples after pruning from distorting per-window statistics.

The motivating case is post-`trim_samples.py` tree sequences where some local
trees retain only a small subset of samples. Per-window statistics such as
mutational load, diversity, Tajima's D, SFS, and coalescence summaries can be
badly biased or unstable in those low-sample intervals.

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

## Open Questions

- Should the first implementation support both `mask` and `excise`, or only
  `mask`?
- Should low-sample intervals be combined with existing remove masks, or kept as
  a separate diagnostic mask?
- Should `--min-samples` mean haploid sample nodes or diploid individuals?
- Should windows be dropped when any overlapping tree is below threshold, or
  weighted by the fraction of the window that passes?
- Should the filter use absolute `N` only, or also support a fraction such as
  `--min-sample-fraction 0.8`?

## Difficulty

Moderate if implemented as mask-only: one new script, a small shared helper,
Snakefile/config wiring, and focused tests.

Higher if implemented as excision immediately, because coordinate compaction,
metadata provenance, and downstream interpretation all need careful handling.
