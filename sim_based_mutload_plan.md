# Simulation-Based Mutational Load Expectation Plan

## Goal

Replace the per-window median mutational-load reference with a per-sample
simulation-based expected load drawn from a single seeded mutation simulation on the input ARG. Apply to
both `scripts/mutload_masks.py` and `scripts/mutload_summary.py`. Keep the
existing `(1 ± cutoff) * expected` outlier semantics.

“Expected” here means one reproducible simulation estimate, not an analytic
expectation or a mean over simulations. The estimate therefore retains Monte
Carlo variance even when the seed makes repeated runs identical.

## Current behavior

Both scripts compute observed load with `argtest_common.mutational_load()`,
aggregate to individuals, then threshold against the per-window median across
individuals (a single scalar per window):

```python
window_medians = np.median(load, axis=1)
high = (1 + cutoff) * window_medians
low  = (1 - cutoff) * window_medians
```

The median drifts when many samples in a window are affected by introgression
or assembly artifacts, so real outliers can hide.

## Implementation

1. **Mutation-rate resolution**, in order:
   - `ratemap_from_metadata(ts.metadata)`
   - `infer_mu_path(ts_path)` → load `*.mut_rate.p`
   - scalar `--mutation-rate` fallback

   Rate may be an `msprime.RateMap` or scalar, matching `--sim-branch` in
   `validation_plots_from_ts.py:362`.

2. **Single-simulation expected load** (use the shared focused mutation-load
   helper so masks and summaries cannot drift):

   ```python
   sim_ts = msprime.sim_mutations(ts, rate=mu, keep=False, random_seed=seed)
   expected = mutational_load(sim_ts, windows=windows)
   expected, _ = aggregate_by_individual(expected, names)
   ```

3. **Threshold replacement** in `mutload_masks.py` (around line 114):

   ```python
   valid = expected > 0
   high  = (1 + args.cutoff) * expected
   low   = (1 - args.cutoff) * expected
   outlier_mask = ((load > high) | (load < low)) & valid
   ```

   Both `expected` and `load` are now `windows × individuals` matrices, so the
   broadcast is element-wise (no `[:, None]` needed).

4. **Mirror edit** in `mutload_summary.py` at lines 206–209.

5. **Output schema (resolved):** the BED's column 5 (`window_medians[w]`) is
   informational only — `pipeline_summary.parse_outlier_bed` reads cols 0–3,
   and `trim_samples.py` via `load_remove_intervals` likewise. Replace col 5
   with a comma-separated list of expected values parallel to the outlier list
   in col 4 (1:1 with the outlier individuals). No compat risk.

6. **CLI options**:
   - `--mutation-rate FLOAT` — scalar fallback when no ratemap is found.
   - `--random-seed INT` — base seed; combined with replicate name (or
     chrom/rep) to derive a deterministic per-replicate seed. Required for
     Snakemake idempotence; without it, re-running step 3 produces different
     outlier sets.

7. **Tests**:
   - Determinism: same `--random-seed` → identical outlier BEDs.
   - Scalar fallback rate path.
   - Metadata ratemap path (if easy to fixture; otherwise rely on the existing
     ratemap-metadata coverage in `test_validation_plots*` style tests).
   - Outlier behavior: a constructed ts where one individual has an obvious
     load excess produces exactly that individual in the outlier BED.

## Non-goals (deferred)

- **Multi-sim averaging / variance estimation.** Mean-of-N sims only trims
  noise on `expected`; the cutoff stays a fixed fraction. Worth doing only
  when adding a variance-aware (z-score) cutoff, which is a separate change.
- **Back-compat `--reference {sim,median}` toggle.** No consumer depends on
  the median path. Replace outright.
- **Analytic or multi-simulation expectation.** The shared helper deliberately
  retains one seeded simulation and the existing scientific semantics.
- **User-supplied demography model.** The simulation rides on the input ARG;
  no separate demography needed.

## Difficulty

Low–moderate. ~30–60 lines net across the two scripts plus tests. The bulk is
rate resolution and the seed plumbing; the threshold change itself is a
handful of lines.
