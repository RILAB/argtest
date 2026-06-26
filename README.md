# ARG Tree Sequence Utilities and Validation Plotting

Snakemake pipeline for post-processing, QC, and basic visualization of ARG tree sequences (`.ts`, `.trees`, `.tsz`). Written with the aid of [Codex](https://openai.com/codex/) and [Claude](https://claude.ai/).

If you use, please cite:

Ross-Ibarra, J. 2026. ARGtest: tools for QC and validation of ancestral recombination graphs. [doi: 10.5281/zenodo.19698118](https://doi.org/10.5281/zenodo.19698118)

## Contents

- [Install](#install)
- [Quick start](#quick-start)
- [Suggested Workflow](#suggested-workflow)
- [Snakemake Pipeline](#snakemake-pipeline)
- [Scripts](#scripts)
- [Auxiliary scripts](#auxiliary-scripts)
- [Inputs, formats, defaults & logs](#inputs-formats-defaults--logs)
- [Repository notes](#repository-notes)
- [Changelog](CHANGELOG.md)

See [NOTES.md](NOTES.md) for the shared module layout and `trim_samples.py` sample-ID matching rules. Run any script with `--help` for its full argument list.

## Install

```bash
conda env create -f environment.yml
conda activate argtest
```

**Which version to use:** work from the most recent tagged release rather than an
arbitrary commit on `main`. To check out the latest tag:

```bash
git fetch --tags
git checkout "$(git tag -l 'v*' | sort -V | tail -1)"
```

See [CHANGELOG.md](CHANGELOG.md) for a per-version breakdown of changes.

Core dependencies are in `environment.yml` (`numpy`, `matplotlib`, `tskit`, `tszip`, `msprime`, `snakemake`, `pytest`).

## Quick start

New here? Run the whole pipeline end-to-end on the bundled example dataset:

```bash
# 1. Create and activate the environment
conda env create -f environment.yml && conda activate argtest

# 2. Generate the example tree sequences. These are not committed; they
#    regenerate deterministically from the seed pinned in ground_truth.json
#    (a few minutes of simulation).
python scripts/make_realistic_example.py --out-dir argtest-realistic-example \
    --n-chrom 3 --n-reps 8 --n-samples 16 --seq-length 10000000

# 3. Run the pipeline. The committed config/snakemake.yaml already points at
#    the example dataset above.
snakemake --cores 16 --configfile config/snakemake.yaml
```

Results land in `argtest-realistic-example-out/`. For a dry-run preview, a
step-by-step explanation, and cluster execution, see
[How To Run](#how-to-run) below.

## Suggested Workflow

One reasonable post-processing workflow for ARG tree sequences in this repo is:

1. **Find low rec regions** Because regions of low recombination are more affected by linked selection, for analyses assuming the neutrality of the ARG it may be a good idea to remove low recombination regions ahead of time. Find windows in the genetic map in the bottom `X` percentile of `cM/Mb` using [scripts/hapmap_low_rec_mask.py](scripts/hapmap_low_rec_mask.py). This turns a HapMap-style recombination map plus a `.fai` into per-chromosome BED masks for very low-recombination regions. The `rec_fraction` cutoff ranks the intervals *between map markers* by recombination rate and masks the lowest fraction of those **intervals** — not the lowest fraction of base pairs. Because low-recombination intervals tend to be long, masking e.g. the bottom 5 % of intervals can remove substantially more than 5 % of the genome by bp. If your ARG was created with [singer-snakemake](https://github.com/nspope/singer-snakemake/tree/main) and the `paircoal-reweight: true` flag, depending on your intereste you might not need to run this this step (set `rec_fraction` to 0). 
2. **Find regions of poor alignment** Find windows of `size` kb where more than `X`% of bp are masked using [scripts/find_low_access_regions.py](scripts/find_low_access_regions.py). This inspects the inferred mutation map for a tree sequence and writes low-accessibility windows to a BED file.
3. **Find aberrant mutational load** All samples in a tree should have similar derived-mutation loads after accounting for local mutation-rate variation. In windows of `number` SNPs or fixed bp size, simulate mutations once on the input ARG to get a per-window per-individual expected load. [scripts/mutload_masks.py](scripts/mutload_masks.py) flags any (window, individual) pair whose observed load differs from its simulated expectation by more than `X`%, writing one BED row per outlier window with the affected individuals comma-joined in column 4 for [scripts/trim_samples.py](scripts/trim_samples.py) to prune those samples from just those windows. The comparison is per individual: a sample's per-window load is summed across all of that individual's sample nodes before the test, so for haploid data (one sample node per individual, e.g. the `admix` dataset) the comparison is effectively per sample, while for diploid/polyploid data the haploid copies are pooled into a single per-individual test. This is what the Snakemake workflow calls. [scripts/mutload_summary.py](scripts/mutload_summary.py) provides an HTML diagnostic that prunes the same (window, individual) pairs, plots each individual's residual observed load as an ASCII bar chart, and highlights in red anyone still outside the cutoff band after pruning as a candidate for manual whole-sample removal. A zero simulated expectation in a window flags any observed-positive load as a high-load outlier.
4. **Remove affected regions** For each chromosome, combine the BED files from steps 1-3 (<chr>.low_rec.mask.bed, low_access.ws<window>.accbp<cutoff>.bed, and <ts_stem>_mutation_masked.bed), then remove those genomic regions from a directory of tree sequences with [scripts/trim_regions.py](scripts/trim_regions.py). This script applies a supplied BED mask and writes trimmed tree sequences.
5. **Remove affected samples** In many cases, only a few samples within a window will be problematic. They could have evidence of introgression (identified using e.g. [TRACE](https://github.com/YulinZhang9806/trace)) or odd patterns of derived mutations (see step 3). Using a bedfile specifying regions and individuals, prune those individuals from the trees with [scripts/trim_samples.py](scripts/trim_samples.py). When run through the Snakemake pipeline this step automatically prunes the step-3 mutload outliers; to also remove your own hand-curated individuals (independent of step 3), set `trim_individuals` (genome-wide IDs) and/or `trim_remove_bed` (per-individual interval BED files) in the config — both are merged with the mutload outliers. Running [scripts/trim_samples.py](scripts/trim_samples.py) directly (`--individuals` / `--remove`) gives the same control outside the pipeline.
6. **Validate** Run the validation plots with [scripts/validation_plots_from_ts.py](scripts/validation_plots_from_ts.py) to get a sense of the cleaned ARG. This gives a compact set of QC plots for mutational load, diversity, Tajima's D, and the site-frequency spectrum (both folded and unfolded). Compare these to the same plots run on the original tree sequences.
7. If satisfied, merge chromosomes for each replicate for downstream analysis using [scripts/merge_treefiles_by_replicate.py](scripts/merge_treefiles_by_replicate.py). This concatenates chromosome-specific tree files into one combined tree sequence per replicate. Mutation-rate ratemaps embedded in each chromosome's metadata are merged and carried forward in the combined output.
8. **Summarise** Generate a self-contained HTML pipeline summary with [scripts/pipeline_summary.py](scripts/pipeline_summary.py). This collects genome retention statistics across all pipeline steps (mean ± SD across replicates for per-replicate steps), per-individual outlier counts from step 5, and embeds the validation plots for every chromosome in a single HTML file.

## Snakemake Pipeline

For a rule-based version of the workflow above, this repository now includes a Snakemake pipeline in [Snakefile](Snakefile).

The Snakemake workflow is designed for a directory layout with one subdirectory per chromosome and one treefile per replicate inside each chromosome directory:

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
```

The Snakemake workflow discovers treefiles with suffixes `.ts`, `.trees`, and `.tsz`. Replicate IDs are taken from the filename stem, so `chr1/1.tsz` is replicate `1` for chromosome `chr1`.

If the treefiles live in a subdirectory of each chromosome directory rather than directly inside it — for example SINGER output where each `chrN/trees/` holds the replicates — set `tree_subdir` to that subdirectory name and discovery looks there instead:

```text
<root>/
  chr1/
    trees/
      chr1.1.tsz
      chr1.2.tsz
      ...
```

The chromosome name still comes from the chromosome directory (`chr1`), and a leading `chrN.` prefix on the filename is stripped to get the replicate ID, so `chr1/trees/chr1.2.tsz` is replicate `2`.

### Required Inputs

The Snakemake config is in [config/snakemake.yaml](config/snakemake.yaml). Edit it for your project and supply it with `--configfile`. The file has an inline comment for every option.

Required keys:

- `root_dir`: path to the chromosome-subdirectory root
- `hapmap`: single HapMap recombination map covering **all** chromosomes (one combined file, not one per chromosome — rows are grouped by the `Chromosome` column), used for step 1
- `fai`: FASTA index used for chromosome lengths
- `rec_fraction`: fraction of recombination-rate **intervals** (ranked by `Rate(cM/Mb)`, lowest first) to include in the low-recombination mask; e.g. `0.1` masks the bottom 10 % of intervals, while `0.0` writes empty low-recombination masks. Note this is a fraction of the *number of intervals between map markers*, not of base pairs — because the lowest-recombination intervals tend to be the longest, masking the bottom 5 % of intervals can remove well over 5 % of the genome by bp
- `low_access_window_size`: window size in bp for step 2
- `low_access_cutoff_bp`: minimum accessible bp per window for step 2
- exactly one of `mutload_window_size` or `mutload_snp_window` for step 3

Optional keys (all have sensible defaults):

- `tree_pattern`: glob for treefiles within each chromosome directory (default: `"*"`), for example `"*.trees"` or `"*.tsz"`
- `tree_subdir`: optional subdirectory within each chromosome directory that holds the treefiles (default: unset → files live directly in the chromosome dir); e.g. `"trees"` for SINGER-style `chrN/trees/` layouts
- `mutload_cutoff`: outlier cutoff fraction for step 3 (default: `0.5`)
- `mutation_rate`: shared scalar mutation rate (per bp per generation) used as the fallback for **both** step 3 and step 6 when no embedded or sibling ratemap is available; set this once instead of the two per-step keys below
- `mutload_mutation_rate`: per-step override of `mutation_rate` for step 3 only (defaults to `mutation_rate` if unset)
- `mutload_random_seed`: base seed for the per-replicate mutation simulation in step 3 (default: `1`)
- `mutload_fraction`: fraction threshold for writing mutation-masked BED rows in step 3
- `suffix_to_strip`: suffix removed from sample IDs before matching in step 3 and step 5 (default: `"_anchorwave"`)
- `trim_individuals`: extra individual IDs removed genome-wide in step 5, **in addition** to the step-3 mutload outliers (e.g. introgressed samples). Comma-separated string (`"id1,id2"`) or a YAML list; matched with `suffix_to_strip` like step 3. Default: unset (trim only mutload outliers)
- `trim_remove_bed`: extra BED file(s) of per-individual intervals removed in step 5, in addition to the mutload outliers. Column 4 holds comma-separated sample IDs (or the filename stem if absent). Single path or a YAML list of paths; applied identically to every (chrom, rep). Default: unset
- `allow_missing_replicates`: set to `true` to concatenate partial replicate sets (default: `false`)
- `burnin`: number of leading discovered replicates to discard before concatenation (default: `0`); must be smaller than the number of replicates discovered after applying `tree_pattern`
- `base_name`: prefix used for merged outputs (default: name of `root_dir`)
- `merged_out_suffix`: force a specific output suffix for merged files (`.ts`, `.trees`, `.tsz`); default is to inherit the suffix of the first input
- `out_dir`: output root for Snakemake products (default: `snakemake_out`; tilde is expanded)
- `validation_mutation_rate`: per-step override of `mutation_rate` for step 6 validation plots (defaults to `mutation_rate` if unset); omit/leave both unset or set to `null` to skip step 6
- `validation_first_chrom_only`: run step 6 only on the first chromosome (default: `true`)
- `validation_sim_branch`: simulate site mutations on each ARG replicate with msprime for a posterior-predictive check instead of scaling branch statistics (default: `false`)
- `emit_vcf`: if `true`, export one `.vcf.gz` per (chromosome, replicate) from the trimmed step-5 tree sequences into `<out_dir>/vcf/` (default: `false`); see [VCF export](#vcf-export) below
- `vcf_reps`: restrict VCF output to specific replicate IDs (a subset of the post-`burnin` replicates); leave unset/null to emit every post-`burnin` replicate

**Where the mutation rate comes from.** Steps 3 and 6 resolve a per-bp mutation rate in this order: (1) a ratemap embedded in the tree-sequence metadata, (2) a sibling `*.mut_rate.p` file near the treefile, (3) the scalar `mutation_rate` (or its per-step override). SINGER output normally provides (1) or (2), so a `*.mut_rate.p` file is *not* required if a ratemap is embedded. **For ARGs produced by non-SINGER software** — which typically carry no embedded ratemap and ship no `*.mut_rate.p` — just set the scalar `mutation_rate` and the pipeline runs without any ratemap file. One caveat: step 3's outlier test is designed to correct for *local* mutation-rate variation, so a flat scalar reduces it to a uniform-rate expectation (no spatial correction); prefer an embedded or sibling ratemap for step 3 when you have one. Step 6's diversity scaling is unaffected by this and works fine with a scalar.

### File naming and what must match

The pipeline derives the chromosome label, the replicate ID, and (optionally) the mutation map from your directory layout and filenames. Two of these must line up with the *contents* of your input files; getting them wrong is the most common setup failure (`KeyError: Chromosome '<label>' not found in <file>` at step 1, or a "could not infer mutation map" error at step 3/6).

- **Chromosome label** — the name of each chromosome subdirectory directly under `root_dir` (e.g. `chr1` in `<root>/chr1/...`; see the layout diagrams above). This one string is written as the BED chromosome column **and** is the key looked up in your hapmap and `.fai`. Keep it short and chromosome-like (`1`, `chr1`); **do not embed run descriptions in the directory name** — a directory called `chr.10.combined.snp.te.sorted` becomes the label verbatim and will not match a normal hapmap. (If your treefiles sit one level deeper, e.g. `chrN/trees/`, set `tree_subdir` rather than pushing `root_dir` down or lengthening the chromosome directory name.)

- **Replicate ID** — the treefile stem with a leading `<chromosome-label>.` prefix stripped, so `chr1/trees/chr1.2.tsz` is replicate `2`. One treefile per replicate per chromosome.

- **Hapmap `Chromosome` column** — the hapmap is one combined file given by the `hapmap:` path; **its filename is irrelevant**. For each chromosome the pipeline looks up the chromosome label in the `Chromosome` column, trying, in order: the label as-is, the substring after its first `.`, then `chr<that>` and `chr_<that>`. So a label `amaranth.1` matches a column value of `1`, `chr1`, or `chr_1`; a bare label `chr1` matches only `chr1` (there is no dot to strip, so `1` alone would *not* match). The **same matching is applied to the first column of the `.fai`.** A label like `chr.10.combined.snp.te.sorted` only ever reduces to `10.combined.snp.te.sorted`, which is why such names fail to match.

- **`*.mut_rate.p` files** (optional — see "Where the mutation rate comes from" above) — unlike the hapmap, these **are** discovered by name, from the treefile path, independent of `root_dir`. For a treefile `<stem>.<rep>.<suffix>` the pipeline searches the treefile's own directory **and one level up** for `<base>.mut_rate.p`, where `<base>` is the stem with the trailing `.<rep>` removed (the full stem is also accepted). Example: for `.../trees/foo.bar.12.tsz`, name the map `foo.bar.mut_rate.p` (one per chromosome) or `foo.bar.12.mut_rate.p` (one per replicate) and place it in `trees/` or its parent. If two different `*.mut_rate.p` files could match, the run aborts with an "Ambiguous mutation maps" error.

### Workflow Steps

The Snakemake workflow runs the same logical steps as the manual workflow:

1. Build per-chromosome low-recombination BED masks
2. Build per-chromosome low-accessibility BED masks
3. Build per-chromosome, per-replicate mutational-load outlier BEDs and optional mutation-masked BEDs
4. Combine the BED masks and trim affected regions from each treefile; the mutation-rate ratemap is embedded in each trimmed treefile's metadata
5. Trim the affected samples from each trimmed treefile
6. Concatenate chromosomes for each replicate into one combined tree sequence; ratemaps are merged across chromosomes and carried forward
7. (optional) Run validation plots on original and cleaned tree sequences and generate a self-contained HTML pipeline summary
8. (optional) Export a per-(chromosome, replicate) VCF from the trimmed tree sequences (`emit_vcf`); see [VCF export](#vcf-export) below

The final merged outputs are named:

```text
<base_name>.combined.<replicate>.<suffix>
```

and are written under the configured `out_dir` in a `combined/` directory.

### VCF export

Set `emit_vcf: true` to write one bgzip-less `.vcf.gz` per `(chromosome, replicate)` under `<out_dir>/vcf/<chrom>/<rep>.vcf.gz`, produced by [`scripts/export_vcf.py`](scripts/export_vcf.py) from the trimmed step-5 tree sequences (so coordinates are real per-chromosome positions, not the concatenated coordinates of the merged ARG). Notes:

- **Variable sites only** — the records are the sites carried on the ARG; the pipeline does not synthesize invariant/monomorphic positions.
- **Pruned samples are missing, not dropped.** Because `trim_samples` leaves a pruned sample *isolated* over its intervals, that sample is written as a missing genotype (`.`) at any site inside those intervals (via tskit's `isolated_as_missing`), while remaining present elsewhere. The sample/site is never globally removed.
- **Ploidy-aware genotypes** — a haploid individual (one sample node, e.g. the `admix` data) is coded as a single allele `0`/`1`; a diploid individual as `0|1`-style. ARG genotypes are phased (`|`).
- **One VCF per replicate** — each replicate is a distinct ARG with its own topology, sites, and per-replicate trimming, so VCFs differ across replicates. The replicate set follows the pipeline `burnin` (leading replicates already dropped); restrict further with `vcf_reps`. To get one genome-wide VCF per replicate, `bcftools concat` the per-chromosome files.

### How To Run

From the repo root:

```bash
module load conda
conda activate argtest
```

Edit `config/snakemake.yaml` to point at your data (at minimum set `root_dir`, `hapmap`, and `fai`), then run.

#### Walk-through using the example dataset

The committed [config/snakemake.yaml](config/snakemake.yaml) points at the realistic example dataset under [argtest-realistic-example/](argtest-realistic-example/) — a deliberately-flawed dataset (contaminated individuals, per-window sample pruning, and an accessibility mask) with a ground-truth scorecard at [argtest-realistic-example-out/scoring_report.md](argtest-realistic-example-out/scoring_report.md) (precision/recall of the pipeline's masks against the injected flaws).

The `.trees` inputs are **not** committed (they regenerate deterministically), so generate them first with [`scripts/make_realistic_example.py`](scripts/make_realistic_example.py) (seed pinned in `ground_truth.json`); see [MAKE_REALISTIC_EXAMPLE.md](MAKE_REALISTIC_EXAMPLE.md) for the generator's CLI and the ground-truth schema:

```bash
python scripts/make_realistic_example.py --out-dir argtest-realistic-example \
    --n-chrom 3 --n-reps 8 --n-samples 16 --seq-length 10000000
```

(These non-default flags reproduce the committed dataset — 3 chromosomes × 8 replicates × 16 diploids × 10 Mb — that `config/snakemake.yaml` expects; the seed and rates match `ground_truth.json`.)

Dry-run first to preview the jobs:

```bash
snakemake -n -p --configfile config/snakemake.yaml
```

Then run the workflow for real:

```bash
snakemake --cores 16 --rerun-incomplete --keep-going --configfile config/snakemake.yaml
```

This is fine for the small example dataset. **For real datasets, prefer the [SLURM route](#running-on-a-slurm-cluster) below.** The `merge_replicates` step loads all chromosomes of a replicate into memory at once and can need ~50–128 GB each; with plain `--cores N`, Snakemake may launch several merges concurrently and OOM the machine. If you must run locally, cap memory with `--resources mem_mb=<node_RAM_in_mb>` so the merges serialize to fit.

> **Sandboxed / read-only `~/.cache` environments** (some HPC or container setups) may need cache and temp-dir redirects — see [Sandboxed environments](NOTES.md#sandboxed-environments-read-only-cache) in NOTES.md.

#### Running on a SLURM cluster

Add `--profile profiles/slurm` to submit each job to SLURM instead of running locally:

```bash
snakemake --profile profiles/slurm --configfile config/snakemake.yaml
```

This is the recommended way to run real datasets: besides parallelism, it fans the ~per-(chromosome, replicate) jobs out across the cluster, sending light steps to one partition and the memory-heavy `merge_replicates` / `step6_validation_plots` steps to a big-mem partition — so each merge gets its own right-sized node and avoids the local-run OOM noted above. **All cluster knobs live in your configfile** — `slurm_account`, `slurm_partition`, and the per-rule `resources:` block (mem/time/threads/partition). The profile file itself ([profiles/slurm/config.yaml](profiles/slurm/config.yaml)) is set-and-forget; you do not edit it.

Notes:
- **Partition names are cluster-specific.** The shipped config uses `"low"` and `"high"`, which are particular to our cluster; change `slurm_partition` (and any per-rule `partition` in the `resources:` block) in your configfile to match the partitions on your own cluster.
- **`out_dir` must be on a shared filesystem** (not node-local `/tmp`) — each job runs on a different node, so a `/tmp` output path silently loses results.
- The snakemake process here is just a controller (it submits jobs and polls), but it runs for the whole pipeline — launch it inside `tmux`/`screen` or a small `srun` so it survives disconnects.
- Per-job SLURM logs are written to `logs/slurm/`.
- Plain `snakemake --cores N …` (above) still runs everything locally; SLURM-only settings such as account, partition, memory, and walltime are ignored, while per-rule `threads` still affects local scheduling.

### Output Layout

By default, Snakemake writes outputs beneath `out_dir` with subdirectories for each stage:

```text
<out_dir>/
  step1_low_rec/
  step2_low_access/
  step3_mutload/
  step4_masks/
  step4_trimmed_regions/
  step5_trimmed_samples/
  combined/
  step6_validation/    # step 6 validation plots (original and cleaned), if configured
  vcf/                 # per-(chromosome, replicate) VCFs, if emit_vcf: true
  pipeline_summary.html
  logs/
```

Intermediate filenames include both chromosome and replicate information so they stay unique across the full workflow.

## Scripts

Pipeline scripts (called by the Snakefile). Run any with `--help` for arguments, defaults, and examples.

- [`hapmap_low_rec_mask.py`](scripts/hapmap_low_rec_mask.py) — per-chromosome BED of the bottom `--rec-fraction` of HapMap recombination-rate intervals.
- [`find_low_access_regions.py`](scripts/find_low_access_regions.py) — BED of low-accessibility windows, computed from a tree sequence's inferred mutation map.
- [`mutload_summary.py`](scripts/mutload_summary.py) — interactive HTML diagnostic for the mutload step: per-individual residual load after window-level pruning (ASCII bar chart, red highlight on individuals still outside the cutoff band, lineage table of flagged counts).
- [`mutload_masks.py`](scripts/mutload_masks.py) — outlier and mutation-masked BED files for one tree sequence (pipeline step 3).
- [`combine_remove_masks.py`](scripts/combine_remove_masks.py) — merge the step 1–3 BED masks into a single combined BED per chromosome.
- [`trim_regions.py`](scripts/trim_regions.py) / [`trim_regions_single.py`](scripts/trim_regions_single.py) — apply a BED mask to a directory (or single file) of tree sequences and write trimmed outputs with compacted coordinates.
- [`trim_samples.py`](scripts/trim_samples.py) — remove individuals genome-wide (`--individuals`) or over BED intervals (`--remove`). See [NOTES.md](NOTES.md) for the exact sample-ID matching rules.
- [`validation_plots_from_ts.py`](scripts/validation_plots_from_ts.py) — SINGER-style QC plots (mutational load, diversity, Tajima's D, folded/unfolded SFS) across TS replicates; optional observed-vs-simulated overlays.
- [`merge_treefiles_by_replicate.py`](scripts/merge_treefiles_by_replicate.py) — concatenate chromosome-specific tree-sequence files by replicate; embedded mutation-rate ratemaps are merged and carried forward.
- [`export_vcf.py`](scripts/export_vcf.py) — export a `.vcf`/`.vcf.gz` from a (filtered) tree sequence: variable sites only, ploidy-aware genotypes, and samples pruned by `trim_samples` written as missing (`.`) via `isolated_as_missing`. Driven by `emit_vcf` (pipeline step 8); see [VCF export](#vcf-export).
- [`pipeline_summary.py`](scripts/pipeline_summary.py) — self-contained HTML report of genome retention, per-individual outlier counts, and embedded validation plots.

### ⚠️ Warning: summary statistics after sample pruning

Diversity (π), Tajima's D, and SFS in this pipeline are computed with [tskit](https://tskit.dev/)'s built-in methods (`ts.diversity`, `ts.Tajimas_D`, `ts.allele_frequency_spectrum`). These normalize by a **constant** `n · (n − 1) / 2` based on the sample set passed in, so when the sample size varies across regions — as it does after `trim_samples.py` removes individuals over BED intervals — the per-window statistics are **not correctly normalized** for the locally retained sample count. Treat post-pruning stats with caution, and prefer comparisons on replicates where sample membership is uniform across the genome.

## Auxiliary scripts

Scripts not called by the Snakemake pipeline.

- [`coalescence_ne_plots_from_ts.py`](scripts/coalescence_ne_plots_from_ts.py) — pair-coalescence and Ne plots from TS replicates. Choose the time grid with either `--time-bins-file` (explicit edges) or `--num-bins N` (equal-coalescence-event bins derived from `pair_coalescence_quantiles` averaged over post-burnin replicates). Optional Demes-based coalescent simulations (`--sim N`) produce window-stat and SFS TSVs for observed-vs-sim density plots in `validation_plots_from_ts.py`.

> ⚠️ **Per-window caveat (pair-coalescence rescale):** `compute_pair_coal` rescales tskit's `pair_coalescence_counts(pair_normalise=True)` PDF by a single scalar (`n_pairs · seq_length / connected_pair_span`) to fix the over-counted denominator when partial trees leave some samples isolated. This is the same bug fixed upstream in [nspope/tskit `nsp-paircoal-partial-missing`](https://github.com/nspope/tskit/tree/nsp-paircoal-partial-missing) (commit `be9f67f2`). The scalar rescale is mathematically exact for the script's current global, single-sample-set, time-binned output, but would **not** be correct if the script were extended to emit spatial-`windows=` rates or multiple-`sample_sets` output — those need per-window or per-index-pair denominators rather than a single scalar. The rescale will be retired in favor of the in-engine fix once it ships in a released `tskit`.
- [`compare_trees_html.py`](scripts/compare_trees_html.py) — render one tree index from each of two tree sequences side-by-side into a single HTML file.
- [`trees_gallery_html.py`](scripts/trees_gallery_html.py) — scrollable HTML gallery of all trees from two tree sequences, useful for quick before/after inspection.
- [`simulate_two_bottleneck_demography.py`](scripts/simulate_two_bottleneck_demography.py) — simulate replicate ARGs under a fixed two-bottleneck demography (35 ka + 9 ka bottlenecks, present-day expansion) for known-truth pipeline tests.
- [`make_realistic_example.py`](scripts/make_realistic_example.py) — generate a realistic synthetic example dataset (ARGs from the two-bottleneck model + three injected flaws: contaminated individuals, per-window sample pruning, and a `mut_rate.p` accessibility mask) for end-to-end pipeline testing. Emits a `ground_truth.json` for scoring the pipeline's masks. See [MAKE_REALISTIC_EXAMPLE.md](MAKE_REALISTIC_EXAMPLE.md) for details, CLI options, and the ground-truth schema.

## Inputs, formats, defaults & logs

- **Tree-sequence files:** scripts accept `.ts`, `.trees`, and `.tsz` files. Loading/writing `.tsz` requires `tszip` to be installed; scripts will raise a clear error if `tszip` is missing when a `.tsz` is used.
- **BED files:** expected as whitespace-separated lines `chrom  start  end  [name]`. `start` and `end` are numeric (half-open intervals `[start, end)`). If a fourth column `name` is present it may list one or more comma-separated sample IDs; if omitted the BED filename stem is used as the sample name. Lines starting with `#` and blank lines are ignored.
- **HapMap recombination maps:** when required (e.g. `hapmap_low_rec_mask.py`), the script expects the HapMap format used by [`msprime.RateMap.read_hapmap`](https://tskit.dev/msprime/docs/stable/api.html#msprime.RateMap.read_hapmap).
- **Glob `--pattern`:** arguments named `--pattern` accept shell-style glob patterns (for example "*.tsz") and are matched against filenames in the supplied directory.
- **Defaults & output locations:** many scripts write to a `results/` directory or to an output directory under the input tree-directory when `--out`/`--out-dir` are not provided. Examples:
  - `trim_samples.py`: default output is `results/<ts_stem>_trimmed.tsz` when `--out` is not given.
  - `trim_regions.py`: default `--out-dir` is `<ts-dir>/trimmed` and default log is `<out-dir>/logs/trim_regions.log`.
  - `mutload_summary.py` writes `results/<name>.html` and `logs/<name>.log`; no BED files are written (use `mutload_masks.py` for BED output).
  - Several plotting scripts write PNG files into `results/` by default; most have `--out` or `--out-dir` flags to override this.
- **Logging & errors:** scripts write summary logs into `logs/` or into the chosen `--out-dir` (see each script's `--log`/`--out-dir` options). Common failure modes include missing `tszip` for `.tsz` I/O, mismatched sequence lengths across chromosome files (checked by `trim_regions.py`), and invalid BED line formats (the loader will raise `ValueError` when a BED line has fewer than 3 columns). When a script prints an `ERROR:` or raises an exception, check the corresponding `logs/` or `<out-dir>/` log file for the detailed run record.

## Repository notes

- Generated `logs/` and `results/` are git-ignored.
- `.DS_Store` is git-ignored.
- Per-release changes are tracked in [CHANGELOG.md](CHANGELOG.md), keyed to the
  `v1.x` git tags.

## Acknowledgements

None of this would be possible without the patient help and advice of Nate Pope. Any errors, bad code, or poor interpretations, however, are my responsbility alone. This repo also uses code from Nate Pope's [singer-snakemake](https://github.com/nspope/singer-snakemake).
