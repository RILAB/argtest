# Chromosome/Replicate Pipeline Guide

This document describes the end-to-end batch pipeline for this repository when data are organized as:

- one subdirectory per chromosome
- multiple tree sequence files per chromosome
- each file named by replicate (for example: `1.tsz`, `2.tsz`, ..., `50.tsz`)

The pipeline runs steps 1-5 for every treefile, then concatenates chromosomes per replicate.

## 1) Expected Input Layout

### Root directory with chromosome subdirectories

```text
<root>/
  chr1/
    1.tsz
    2.tsz
    ...
  chr2/
    1.tsz
    2.tsz
    ...
  ...
  chr10/
    1.tsz
    2.tsz
    ...
```

Supported tree suffixes:

- `.tree`
- `.trees`
- `.tsz`

Replicate ID is taken from the filename stem (for example `1` from `1.tsz`).

### Required metadata/reference inputs

- HapMap recombination map (`--hapmap`) with columns:
  - `Chromosome`
  - `Position(bp)`
  - `Rate(cM/Mb)`
- FASTA index (`--fai`) with chromosome lengths
- Mutation-rate map(s) (`*.mut_rate.p`) discoverable near each chromosome directory

For replicate-named files like `chr1/1.tsz`, mutation-map discovery checks:

- `<root>/chr1/1.mut_rate.p`
- `<root>/chr1/chr1.mut_rate.p`
- `<root>/1.mut_rate.p`
- `<root>/chr1.mut_rate.p`

Recommended pattern: one mutation map per chromosome at `<root>/<chromosome>.mut_rate.p` (for example `chr1.mut_rate.p`).

## 2) Environment Setup

From the repo root:

```bash
module load conda
conda activate argtest
```

If needed:

```bash
conda env create -f environment.yml
module load conda
conda activate argtest
```

## 3) Run the Full Pipeline (Steps 1-5 + Per-Replicate Concatenation)

Script:

- `scripts/run_steps1_5_and_concat.py`

Example:

```bash
python scripts/run_steps1_5_and_concat.py \
  --root /path/to/root_with_chr_dirs \
  --hapmap /path/to/map.hapmap.tsv \
  --fai /path/to/genome.fai \
  --rec-fraction 0.1 \
  --window-size 50000 \
  --cutoff-bp 2500 \
  --mutload-window-size 50000 \
  --mutload-cutoff 0.25 \
  --mutload-fraction 0.2 \
  --pattern "*.tsz"
```

Key required arguments:

- `--root`: chromosome-subdirectory root
- `--hapmap`, `--fai`
- `--rec-fraction`
- `--window-size`, `--cutoff-bp`
- one of:
  - `--mutload-window-size`
  - `--mutload-snp-window`

Common optional arguments:

- `--pattern` (default `*`)
- `--out-dir` (default `<root>/batch_steps1_5`)
- `--merged-dir` (default `<out-dir>/combined`)
- `--base-name` (default: root directory name; used in merged filenames)
- `--out-suffix` (`.tree`, `.trees`, `.tsz`)
- `--allow-missing-replicates` (by default, each replicate must exist in every chromosome directory)
- `--suffix-to-strip` (default `_anchorwave`)

## 4) What Each Step Does in This Batch Pipeline

1. Low-recombination mask (per chromosome):
   - builds `<chrom>.low_rec.mask.bed` from `--hapmap` + `--fai`
2. Low-accessibility mask (per chromosome):
   - uses the first treefile in that chromosome directory plus inferred `*.mut_rate.p`
3. Mutational-load outliers (per chromosome x replicate):
   - writes per-replicate outlier BED
   - optionally writes mutation-masked BED when `--mutload-fraction` is provided
4. Region trimming (per chromosome x replicate):
   - combines step 1 + step 2 + optional step 3 mutation-masked BED
   - trims those regions from each treefile
5. Sample trimming (per chromosome x replicate):
   - removes outlier individuals using the step-3 outlier BED
6. Concatenate per replicate across chromosomes:
   - output name: `<base-name>.combined.<replicate><suffix>`

Chromosomes are concatenated in natural sort order (`chr1`, `chr2`, ..., `chr10`).

## 5) Output Structure

By default (`--out-dir` omitted), outputs are under:

```text
<root>/batch_steps1_5/
  step1_low_rec/
  step2_low_access/<chrom>/
  step3_mutload/<chrom>/
  step4_masks/<chrom>/
  step4_trimmed_regions/<chrom>/
  step5_trimmed_samples/<chrom>/
  combined/
  run_summary.log
```

Final concatenated files are in:

- `<root>/batch_steps1_5/combined/`

## 6) Merge-Only Mode (No Steps 1-5)

If you already have per-chromosome treefiles and only want per-replicate concatenation:

```bash
python scripts/merge_treefiles_by_replicate.py \
  --ts-dir /path/to/root_with_chr_dirs \
  --layout nested \
  --base-name myrun
```

`merge_treefiles_by_replicate.py` supports:

- `--layout nested`: `<chromosome>/<replicate><suffix>`
- `--layout flat`: `<base>.<chromosome>.<replicate><suffix>` in one directory
- `--layout auto`: detect either layout automatically

## 7) Practical Notes

- Chromosome directory names must match chromosome names in both HapMap and FAI.
- Replicate IDs are matched by exact filename stem (for example `1`, `2`, ...).
- Without `--allow-missing-replicates`, concatenation fails if a replicate is missing in any chromosome.
- For stable behavior, keep all replicates for a chromosome in the same suffix format.
- If mutation-map inference is ambiguous (multiple matching `*.mut_rate.p`), rename files to remove ambiguity.
