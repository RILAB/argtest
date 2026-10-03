# Pipeline scoring report

Comparison of the Snakemake pipeline output against `ground_truth.json` for the synthetic dataset emitted by `scripts/make_realistic_example.py`.

**Run summary**

- Chromosomes: 3
- Replicates per chromosome: 8
- Diploid individuals: 16
- Sequence length: 10.0 Mb / chrom
- Contaminated individuals: ind1 (sev 2×), ind12 (sev 5×)

## 1. Accessibility detection (step2_low_access)

Does the pipeline recover the regions zeroed in `mut_rate.p`? BP-level precision / recall of the union of per-replicate `step2_low_access/chr*/chr*.low_access.bed` against the ground-truth `masked_intervals` for each chromosome.

| chrom | truth (Mb) | called (Mb) | overlap (Mb) | precision | recall |
|---|---:|---:|---:|---:|---:|
| chr1 | 3.30 | 3.15 | 3.15 | 1.000 | 0.955 |
| chr2 | 3.30 | 3.15 | 3.15 | 1.000 | 0.955 |
| chr3 | 3.30 | 3.15 | 3.15 | 1.000 | 0.955 |

Precision is essentially perfect; recall is bounded by the 50 kb window grid not tiling the exact mask boundaries — windows partially overlapping the mask miss the `low_access_cutoff_bp` cutoff.

## 2. Mutload outlier individuals (step3_mutload outliers.bed)

How often is each individual flagged as a per-window outlier? Summed across the 8 replicates per chromosome. Contaminated individuals should dominate by orders of magnitude.

| chrom | individual | flag count | role |
|---|---|---:|---|
| chr1 | ind12 | 11167 | **CONTAM sev=5×** |
| chr1 | ind1 | 8750 | **CONTAM sev=2×** |
| chr1 | ind10 | 42 | clean |
| chr1 | ind14 | 39 | clean |
| chr1 | ind6 | 37 | clean |
| chr1 | ind13 | 36 | clean |
| chr2 | ind12 | 11838 | **CONTAM sev=5×** |
| chr2 | ind1 | 8911 | **CONTAM sev=2×** |
| chr2 | ind10 | 40 | clean |
| chr2 | ind3 | 34 | clean |
| chr2 | ind14 | 33 | clean |
| chr2 | ind4 | 33 | clean |
| chr3 | ind12 | 11391 | **CONTAM sev=5×** |
| chr3 | ind1 | 8943 | **CONTAM sev=2×** |
| chr3 | ind2 | 37 | clean |
| chr3 | ind9 | 34 | clean |
| chr3 | ind3 | 32 | clean |
| chr3 | ind15 | 31 | clean |

## 3. mutation_masked windows vs prune intervals

When the per-window outlier-fraction exceeds `--mutload_fraction`, the window is written to `mutation_masked.bed` rather than the per-individual outlier BED. The dominant cause is per-window sample pruning: dropped samples report load = 0 vs expected > 0, which flags them on the low end of the cutoff. So `mutation_masked.bed` should overlap the ground-truth `prune_intervals`.

| chrom | prune truth (Mb) | called (Mb) | overlap (Mb) | precision | recall |
|---|---:|---:|---:|---:|---:|
| chr1 | 5.00 | 2.56 | 1.40 | 0.548 | 0.281 |
| chr2 | 5.00 | 2.50 | 1.36 | 0.544 | 0.271 |
| chr3 | 5.00 | 2.49 | 1.46 | 0.585 | 0.291 |

Modest recall: with `--prune-frac-samples 0.25`, exactly 25% of samples are dropped in each pruned window — sitting right at the default `--mutload_fraction 0.2` boundary. Poisson variation pushes some windows over and some under. Tightening `mutload_fraction` to 0.15 would lift recall.

## 4. trim_samples recovery (step5_trimmed_samples)

For each ground-truth `(window, drop_sample_ids)` entry, count how many isolated samples the post-trim ts has at the window's midpoint. Cap at `len(drop_sample_ids)` per entry.

| chrom | expected isolations | recovered | rate |
|---|---:|---:|---:|
| chr1 | 1600 | 1600 | 1.000 |
| chr2 | 1600 | 1600 | 1.000 |
| chr3 | 1600 | 1600 | 1.000 |

## 5. Mutload outlier per-(window, individual) FP rate

Each row in `outliers.bed` flags one or more individuals in one window. Treating each (window, individual) flag as one draw:

- **TP**: flag on a contaminated individual.
- **Spurious**: flag on a non-contaminated individual (any cause).
- **Prune-explained**: subset of spurious flags on samples pruned in that window (load = 0 vs expected > 0 always trips the low end of the cutoff band).
- **Net FP**: spurious flags that are NOT prune-explained — the actual statistical-noise floor of the outlier test.

| chrom | windows | TP | spurious | prune-expl | net FP | spurious rate | net FP rate |
|---|---:|---:|---:|---:|---:|---:|---:|
| chr1 | 11714 | 19917 | 460 | 222 | 238 | 0.28% | 0.18% |
| chr2 | 11839 | 20749 | 371 | 159 | 212 | 0.22% | 0.16% |
| chr3 | 11747 | 20334 | 399 | 186 | 213 | 0.24% | 0.17% |
| **total** | **35300** | **61000** | **1230** | **567** | **663** | **0.25%** | **0.17%** |

**TP recall**: 61000 / 70600 expected = 86.40% — i.e., when an outlier window exists, the contaminated individuals are flagged in it ~86% of the time.

**Denominators**:
- spurious rate = spurious / (windows × 14 non-contam inds)
- net FP rate = net_FP / (windows × non-contam inds NOT pruned in that window)

Net FP rate ~0.2% means ~1 in 500 (window, non-contam, non-pruned) cells is flagged by chance. This is the actual statistical noise floor of the per-window cutoff test at the `--mutload_cutoff` value used in this run; well below the 5-10% you'd see from a coarser window-median scheme.

Note: `outliers.bed` only enumerates windows with at least one flag; windows with zero flags don't appear. Denominators here are conditional on the window having ≥1 flag. The unconditional per-cell FP rate would be lower.
