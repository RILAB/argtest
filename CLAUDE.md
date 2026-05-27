# Project notes for Claude

## Pipeline overview

Snakemake pipeline for analyzing tree sequences from amaranth data (16 chromosomes, 150 replicates each). Steps:
1. `step1_low_rec_masks` — hapmap-based low-recombination masks (`hapmap_low_rec_mask.py`)
2. `step2_low_access` — low-accessibility window masks (`find_low_access_regions.py`)
3. `step3_mutload_masks` — mutation load outlier masks (`mutload_masks.py`)
4. `step4_combine_remove_masks` / `step4_trim_regions_single` — merge masks and trim tree sequences
5. `step5_trim_samples_single` — trim samples

Input data lives in `amaranth/` (one subdir per chromosome, e.g. `amaranth/amaranth.1/`). Results go to `amaranth_results/` (path set by `out_dir` in `config/snakemake.yaml`, tilde-expanded).

## Merge step memory and cluster execution

`merge_replicates` (`merge_treefiles_by_replicate.py`) loads ALL chromosomes of a replicate into memory at once and `.concatenate()`s them, so its peak RAM is far higher than any other step (steps 1-5 each hold a single tree sequence). Empirically, on the `admix/` data (10 chroms, 24 samples, 2.13 Gbp, ~101M edges per replicate): whole-genome tables ≈ 8 GB, but the merge transient peak is ~52 GB and it **OOM-killed at `--mem=64G`** — budget ~96-128 GB per merge (more for larger genomes / higher sample counts).

**Running plain `snakemake --cores N` (the README command) is OOM-prone at the merge stage.** The Snakefile declares no `resources: mem_mb`, so `--cores N` schedules up to N jobs concurrently based only on CPU. Once step5 finishes, every per-replicate merge becomes ready and Snakemake will launch N of them at once — N × ~50-128 GB, which exceeds even a 1 TB node. Mitigations: (a) declare `resources: mem_mb` per rule and run `snakemake --cores N --resources mem_mb=<node_RAM_in_mb>` so Snakemake serializes the memory-heavy jobs to fit the node; or (b) submit via a SLURM profile so each rule is its own right-sized job. Both are tracked in project memory (SLURM-ification todo).

`measure_merge_mem.py` (repo root) measures the actual merge peak RSS for one replicate without writing anything to disk: `python measure_merge_mem.py [REP] [GLOB]` (GLOB has a single `{rep}` field; default is the `admix/` layout). Use it to size `mem_mb` for a new dataset before a full run.

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
