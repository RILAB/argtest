# Project notes for Claude

## Pipeline overview

Snakemake pipeline for analyzing tree sequences from amaranth data (16 chromosomes, 150 replicates each). Steps:
1. `step1_low_rec_masks` — hapmap-based low-recombination masks (`hapmap_low_rec_mask.py`)
2. `step2_low_access` — low-accessibility window masks (`find_low_access_regions.py`)
3. `step3_mutload_masks` — mutation load outlier masks (`mutload_masks.py`)
4. `step4_combine_remove_masks` / `step4_trim_regions_single` — merge masks and trim tree sequences
5. `step5_trim_samples_single` — trim samples

Input data lives in `amaranth/` (one subdir per chromosome, e.g. `amaranth/amaranth.1/`). Results go to `amaranth_results/` (path set by `out_dir` in `config/snakemake.yaml`, tilde-expanded).

## Known issues diagnosed (2026-04-09)

- **Tilde in OUT_DIR**: `out_dir` in config was not tilde-expanded, producing `~/src/argtest/...` paths that Python could not open. Fixed in commit `b9a4722`.
- **Hapmap repeated headers**: `hapmap_low_rec_mask.py` failed on repeated header lines and cross-naming-convention chroms. Fixed in commit `371b526`.
- **`find_low_access_regions.py` interface change**: script argument was renamed from `--ts-dir` (directory + `--pattern`) to `--ts` (single file path). Snakefile already uses `--ts`; uncommitted working-tree changes to the script bring it in sync.
- **Missing `amaranth.16.mut_rate.p`**: mutation rate file was absent for chromosome 16, causing `infer_mu_path` to raise "Ambiguous mutation maps" (it fell back to a glob on the broad base `"amaranth"`, matching all 15 other chroms). File was added manually to `amaranth/` to fix.

## `infer_mu_path` known fragility

`argtest_common::infer_mu_path` (moved from `find_low_access_regions.py` 2026-04-16) derives candidate base names by stripping numeric suffixes from the ts stem and parent dir. For chroms named `amaranth.N`, this adds `"amaranth"` as a fallback base, which matches every `*.mut_rate.p` in the directory. If the exact file (`amaranth/amaranth.N.mut_rate.p`) is missing, the error message will say "Ambiguous" rather than "not found" — check for a missing `mut_rate.p` first.

## Upstream reference

This pipeline is a fork/adaptation of https://github.com/nspope/singer-snakemake. When diagnosing differences in plot output or statistics, compare against that repo's `workflow/scripts/diagnostics.py` and `workflow/scripts/tree_statistics.py`. Key difference: the original simulates mutations via `msprime.sim_mutations` and uses site-mode stats; our scripts currently use branch-mode stats directly (see project memory for details).

## Conda environment

Environment name is in `environment.yml`. Activate with:
```
module load conda && conda activate <name-from-environment.yml>
```
