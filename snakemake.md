# Snakemake Workflow

This repository now includes a rule-based Snakemake pipeline in `Snakefile` for the layout:

```text
<root>/
  chr1/1.tsz  chr1/2.tsz ...
  chr2/1.tsz  chr2/2.tsz ...
  ...
```

It runs:

1. low-recombination masking
2. low-accessibility masking
3. mutational-load outlier masking
4. region trimming
5. sample trimming
6. per-replicate concatenation across chromosomes

## Setup

```bash
module load conda
conda activate argtest
```

## Configure

Copy and edit:

```bash
cp config/snakemake.example.yaml config/snakemake.yaml
```

Set at least:

- `root_dir`
- `hapmap`
- `fai`

And choose exactly one:

- `mutload_window_size`
- `mutload_snp_window`

## Run

Dry run:

```bash
snakemake -n -p --configfile config/snakemake.yaml
```

Execute with parallel jobs:

```bash
snakemake --cores 16 --rerun-incomplete --keep-going --configfile config/snakemake.yaml
```

## Outputs

Outputs are written under `out_dir` (from config), with these key directories:

- `step1_low_rec/`
- `step2_low_access/`
- `step3_mutload/`
- `step4_masks/`
- `step4_trimmed_regions/`
- `step5_trimmed_samples/`
- `combined/` (final concatenated outputs)

Final concatenated files are named:

- `<base_name>.combined.<replicate>.<ext>`

## Notes

- By default, `allow_missing_replicates: false` enforces complete chromosome sets for each replicate.
- Set `allow_missing_replicates: true` to concatenate partial replicate sets.
- Mutation maps (`*.mut_rate.p`) are inferred from tree filenames/directories using existing repo logic.
