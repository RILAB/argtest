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

`tskit` is not installed from conda-forge: `environment.yml` pins it to
[nspope/tskit commit `73d8cd9`](https://github.com/nspope/tskit/commit/73d8cd922482475020ae01180cae95bf5abbf067)
and installs it with pip, so the build needs pip, git, and access to GitHub. The
pinned build is what provides the native partial-missing-data pair-coalescence
normalization (see [Auxiliary scripts](#auxiliary-scripts)) and it requires
Python ≥ 3.11 and NumPy ≥ 2. When upgrading from an environment built before
v1.9, recreate it rather than updating in place, so the conda-forge `tskit` is
replaced:

```bash
conda env remove -n argtest
conda env create -f environment.yml
```

**Which version to use:** work from the most recent tagged release rather than an
arbitrary commit on `main`. To check out the latest tag:

```bash
git fetch --tags
git checkout "$(git tag -l 'v*' | sort -V | tail -1)"
```

See [CHANGELOG.md](CHANGELOG.md) for a per-version breakdown of changes.

## Quick start

New here? After [installing](#install), run the whole pipeline end-to-end on the bundled example dataset:

```bash
# 1. Generate the example tree sequences. These are not committed; they
#    regenerate deterministically from the seed pinned in ground_truth.json
#    (3 chromosomes × 8 replicates × 16 diploids × 10 Mb; a few minutes).
python scripts/make_realistic_example.py --out-dir argtest-realistic-example \
    --n-chrom 3 --n-reps 8 --n-samples 16 --seq-length 10000000

# 2. Run the pipeline. The committed config/snakemake.yaml already points at
#    the example dataset above.
snakemake --cores 16 --configfile config/snakemake.yaml
```

Results land in `argtest-realistic-example-out/`. The example is a deliberately-flawed dataset (contaminated individuals, per-window sample pruning, an accessibility mask) with a ground-truth scorecard at [scoring_report.md](argtest-realistic-example-out/scoring_report.md) grading the pipeline's masks against the injected flaws; see [MAKE_REALISTIC_EXAMPLE.md](MAKE_REALISTIC_EXAMPLE.md) for the generator's CLI and ground-truth schema. For dry-run preview, run flags, and cluster execution, see [How To Run](#how-to-run) below.

## Suggested Workflow

One reasonable post-processing workflow for ARG tree sequences is below. The [Snakemake pipeline](#snakemake-pipeline) automates exactly these steps; each can also be run by hand via the linked script (`--help` for full arguments).

1. **Find low-recombination regions.** Regions of low recombination are more affected by linked selection, so for analyses assuming a neutral ARG it can help to remove them first. [scripts/hapmap_low_rec_mask.py](scripts/hapmap_low_rec_mask.py) turns a HapMap recombination map plus a `.fai` into per-chromosome BED masks of the lowest-recombination intervals (set by `rec_fraction`; see [Required Inputs](#required-inputs) for the intervals-vs-base-pairs caveat). If your ARG came from [singer-snakemake](https://github.com/nspope/singer-snakemake) with `paircoal-reweight: true`, you may be able to skip this step (set `rec_fraction` to 0).
2. **Find regions of poor alignment.** [scripts/find_low_access_regions.py](scripts/find_low_access_regions.py) inspects a tree sequence's inferred mutation map and writes windows where more than a cutoff fraction of base pairs are masked to a BED file.
3. **Find aberrant mutational load.** All samples should have similar derived-mutation loads after accounting for local mutation-rate variation. [scripts/mutload_masks.py](scripts/mutload_masks.py) simulates mutations once on the input ARG to get a per-window, per-individual expected load, then flags any (window, individual) pair whose observed load deviates from expectation by more than `mutload_cutoff`, writing one BED row per outlier window with the affected individuals in column 4. The test is per individual: a diploid/polyploid individual's haploid copies are pooled before testing, so for haploid data it is effectively per sample. [scripts/mutload_summary.py](scripts/mutload_summary.py) gives an HTML diagnostic of residual load after pruning, flagging candidates for whole-sample removal.
4. **Remove affected regions.** Combine the step 1–3 BED masks per chromosome and trim those regions from each tree sequence with [scripts/trim_regions.py](scripts/trim_regions.py). The mutation-rate ratemap is embedded in each trimmed tree sequence's metadata.
5. **Remove affected samples.** Often only a few samples in a window are problematic — e.g. introgressed (cf. [TRACE](https://github.com/YulinZhang9806/trace)) or with odd derived-mutation patterns (step 3). [scripts/trim_samples.py](scripts/trim_samples.py) prunes given individuals over given intervals. The pipeline automatically prunes the step-3 outliers; add your own with `trim_individuals` (genome-wide IDs) and/or `trim_remove_bed` (per-individual BEDs), or run the script directly with `--individuals` / `--remove`.
6. **Validate** (optional). [scripts/validation_plots_from_ts.py](scripts/validation_plots_from_ts.py) produces compact QC plots (mutational load, diversity, Tajima's D, folded/unfolded SFS) for the cleaned ARG; compare against the same plots on the originals.
7. **Merge chromosomes per replicate.** [scripts/merge_treefiles_by_replicate.py](scripts/merge_treefiles_by_replicate.py) concatenates the per-chromosome tree files into one tree sequence per replicate, merging the embedded mutation-rate ratemaps across chromosomes.
8. **Summarise** (optional). [scripts/pipeline_summary.py](scripts/pipeline_summary.py) collects genome-retention statistics across all steps (mean ± SD across replicates), per-individual outlier counts, and the validation plots into one self-contained HTML file. Final retained Mb are measured directly from the filtered ARGs as positive-rate accessible sequence with a non-empty local genealogy, so they exclude both regions masked in the input ARG and regions removed by this pipeline. Retention percentages use the full FASTA-index reference length as the denominator.

Optionally, set `emit_vcf: true` to export a per-(chromosome, replicate) VCF from the filtered per-chromosome tree sequences; see [VCF export](#vcf-export).

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
- `mutation_rate`: single scalar mutation rate (per bp per generation), the shared fallback for **both** step 3 and step 6 when no embedded or sibling ratemap is available
- `mutload_random_seed`: base seed for the per-replicate mutation simulation in step 3 (default: `1`)
- `mutload_fraction`: fraction threshold for writing mutation-masked BED rows in step 3
- `suffix_to_strip`: suffix removed from sample IDs before matching in step 3 and step 5 (default: `"_anchorwave"`)
- `trim_individuals`: extra individual IDs removed genome-wide in step 5, **in addition** to the step-3 mutload outliers (e.g. introgressed samples). Comma-separated string (`"id1,id2"`) or a YAML list; matched with `suffix_to_strip` like step 3. Default: unset (trim only mutload outliers)
- `trim_remove_bed`: extra BED file(s) of per-individual intervals removed in step 5, in addition to the mutload outliers. Column 4 holds comma-separated sample IDs (or the filename stem if absent). Single path or a YAML list of paths; applied identically to every (chrom, rep). Default: unset
- `min_samples`: minimum number of non-isolated **retained sample nodes** (haploids, not individuals) a local tree must have; intervals below it are dropped in an optional **step 5b** that runs after step 5 and before the merge. Dropping uses `delete_intervals`, so sequence coordinates are **preserved** (the removed spans become empty gaps — there is no coordinate compaction) and the `kept_intervals` metadata is intersected with the surviving spans so downstream accessibility is not overestimated. Each `(chrom, rep)` also gets a diagnostic BED of the dropped intervals (`chrom start end retained_samples min_samples`) under `<out_dir>/step5b_min_samples/`. When set, the merge, VCF export, and validation steps automatically consume the step-5b filtered tree sequences instead of the raw step-5 output. Default: unset/null (skip step 5b entirely, so existing configs are unaffected). Sample pruning in step 5 is what creates the low-sample intervals this filter targets
- `allow_missing_replicates`: set to `true` to concatenate partial replicate sets (default: `false`)
- `burnin`: number of leading discovered replicates to discard before concatenation (default: `0`); must be smaller than the number of replicates discovered after applying `tree_pattern`
- `base_name`: prefix used for merged outputs (default: name of `root_dir`)
- `merged_out_suffix`: force a specific output suffix for merged files (`.ts`, `.trees`, `.tsz`); default is to inherit the suffix of the first input
- `out_dir`: output root for Snakemake products (default: `snakemake_out`; tilde is expanded)
- `run_validation`: master switch for step 6 (default: `true`); set `false` to skip the validation plots while keeping `mutation_rate` set for step 3. Step 6 also auto-skips when no rate source is available (neither `mutation_rate` nor `validation_sim_branch`)
- `validation_first_chrom_only`: run step 6 only on the first chromosome (default: `true`)
- `validation_window_size`: window size in bp for step-6 diversity, Tajima's D, and segregating-sites validation plots (default: `100000`); larger windows run faster and use less memory on large ARGs but give coarser QC curves
- `validation_sim_branch`: simulate site mutations on each ARG replicate with msprime for a posterior-predictive check instead of scaling branch statistics (default: `false`); can run without a scalar `mutation_rate` when every validated tree sequence has an embedded/sibling ratemap
- `emit_vcf`: if `true`, export one `.vcf.gz` per (chromosome, replicate) from the filtered per-chromosome tree sequences into `<out_dir>/vcf/` (default: `false`); see [VCF export](#vcf-export) below
- `vcf_reps`: restrict VCF output to specific replicate IDs (a subset of the post-`burnin` replicates); leave unset/null to emit every post-`burnin` replicate

**Where the mutation rate comes from.** Steps 3 and 6 resolve a per-bp mutation rate in this order: (1) a ratemap embedded in the tree-sequence metadata, (2) a sibling `*.mut_rate.p` file near the treefile, (3) the scalar `mutation_rate`. SINGER output normally provides (1) or (2), so a `*.mut_rate.p` file is *not* required if a ratemap is embedded. **For ARGs produced by non-SINGER software** — which typically carry no embedded ratemap and ship no `*.mut_rate.p` — just set the scalar `mutation_rate` and the pipeline runs without any ratemap file. One caveat: step 3's outlier test is designed to correct for *local* mutation-rate variation, so a flat scalar reduces it to a uniform-rate expectation (no spatial correction); prefer an embedded or sibling ratemap for step 3 when you have one. Step 6 can run in `validation_sim_branch` mode with ratemaps alone; otherwise it needs a scalar rate for simulation or branch scaling.

### File naming and what must match

The pipeline derives the chromosome label, the replicate ID, and (optionally) the mutation map from your directory layout and filenames. Two of these must line up with the *contents* of your input files; getting them wrong is a common setup failure.

**Startup chromosome validation.** During Snakefile evaluation, before the DAG is built or any local/SLURM job is submitted, the pipeline checks every chromosome label discovered from the ARG directory layout against both the HapMap `Chromosome` column and the first column of the `.fai`. If a label cannot be resolved, the workflow stops once with a `Chromosome naming mismatch detected before job execution` error that lists all unmatched pipeline chromosomes and the available reference names. For example, ARG directories named `chrom1` through `chrom10` do not match HapMap/FAI names `1` through `10` under the aliases described below, so the complete mismatch is reported at startup instead of failing later in ten separate step-1 jobs.

The workflow has no input VCF whose chromosome names need an independent check. When `emit_vcf: true`, each output VCF's `CHROM` value and contig header are generated from the already-validated ARG chromosome-directory label, so the ARG-derived output VCF label agrees by construction. The startup check therefore covers **ARG directory labels ↔ HapMap ↔ FAI**; it does not inspect chromosome metadata embedded inside a tree sequence.

- **Chromosome label** — the name of each chromosome subdirectory directly under `root_dir` (e.g. `chr1` in `<root>/chr1/...`; see the layout diagrams above). This one string is written as the BED chromosome column **and** is the key looked up in your hapmap and `.fai`. Keep it short and chromosome-like (`1`, `chr1`); **do not embed run descriptions in the directory name** — a directory called `chr.10.combined.snp.te.sorted` becomes the label verbatim and will not match a normal hapmap. (If your treefiles sit one level deeper, e.g. `chrN/trees/`, set `tree_subdir` rather than pushing `root_dir` down or lengthening the chromosome directory name.)

- **Replicate ID** — the treefile stem with a leading `<chromosome-label>.` prefix stripped, so `chr1/trees/chr1.2.tsz` is replicate `2`. One treefile per replicate per chromosome.

- **Hapmap `Chromosome` column** — the hapmap is one combined file given by the `hapmap:` path; **its filename is irrelevant**. For each chromosome the pipeline looks up the chromosome label in the `Chromosome` column, trying, in order: the label as-is, the substring after its first `.`, then `chr<that>` and `chr_<that>`. So a label `amaranth.1` matches a column value of `1`, `chr1`, or `chr_1`; a bare label `chr1` matches only `chr1` (there is no dot to strip, so `1` alone would *not* match). The **same matching is applied to the first column of the `.fai`.** A label like `chr.10.combined.snp.te.sorted` only ever reduces to `10.combined.snp.te.sorted`, which is why such names fail to match.

- **`*.mut_rate.p` files** (optional — see "Where the mutation rate comes from" above) — unlike the hapmap, these **are** discovered by name, from the treefile path, independent of `root_dir`. For a treefile `<stem>.<rep>.<suffix>` the pipeline searches the treefile's own directory **and one level up** for `<base>.mut_rate.p`, where `<base>` is the stem with the trailing `.<rep>` removed (the full stem is also accepted). Example: for `.../trees/foo.bar.12.tsz`, name the map `foo.bar.mut_rate.p` (one per chromosome) or `foo.bar.12.mut_rate.p` (one per replicate) and place it in `trees/` or its parent. If two different `*.mut_rate.p` files could match, the run aborts with an "Ambiguous mutation maps" error.

### Output naming

The Snakemake workflow automates the [Suggested Workflow](#suggested-workflow) steps above. The final merged outputs are named:

```text
<base_name>.combined.<replicate>.<suffix>
```

and are written under the configured `out_dir` in a `combined/` directory.

### Locating a tree by chromosome and position

The merge lays chromosomes end-to-end along a single coordinate axis, so a within-chromosome position (e.g. position 1234 on chromosome 8) sits at `chromosome_offset + position` in the merged sequence. To spare you that arithmetic, the merge records a `chrom_offsets` table (`[{chrom, offset, length}, ...]`) in the merged tree sequence's top-level metadata.

Use [`scripts/locate_tree.py`](scripts/locate_tree.py) to find the local tree at a `(chromosome, position)`:

```bash
python scripts/locate_tree.py --ts <out_dir>/combined/<base>.combined.<rep>.tsz --chrom 8 --position 1234
```

It prints the genome coordinate, the covering tree's index/interval, and flags when the position falls in a masked/trimmed region (a tree with no topology, `num_edges == 0`). The same mapping is available programmatically in [`scripts/argtest_common.py`](scripts/argtest_common.py):

```python
from argtest_common import load_ts, tree_at_chrom_position, genome_position, chrom_position_from_genome
ts = load_ts("combined/run.combined.101.tsz")
tree = tree_at_chrom_position(ts, "8", 1234)   # local tree at chr8:1234
gpos = genome_position(ts, "8", 1234)          # concatenated coordinate
chrom, pos = chrom_position_from_genome(ts, gpos)  # inverse mapping
```

**Merged files produced before this feature** lack `chrom_offsets` and will raise a `KeyError` asking you to re-merge. For those, either read the local tree straight from the per-chromosome file (native coordinates, no offset needed) — `tszip.decompress(".../step5_trimmed_samples/<chrom>/<rep>.tsz").at(position)` — or add the offset by hand. Since trimming preserves each chromosome's length, the offset table is identical across replicates, so you only need to compute it once (sum the per-chromosome `sequence_length`s in natural chromosome order).

### VCF export

Set `emit_vcf: true` to write one bgzip-less `.vcf.gz` per `(chromosome, replicate)` under `<out_dir>/vcf/<chrom>/<rep>.vcf.gz`, produced by [`scripts/export_vcf.py`](scripts/export_vcf.py) from the filtered per-chromosome tree sequences: step 5 output by default, or step 5b output when `min_samples` is enabled. Coordinates are real per-chromosome positions, not the concatenated coordinates of the merged ARG. Notes:

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

Preview the planned jobs with a dry run, then run for real:

```bash
snakemake -n -p --configfile config/snakemake.yaml
snakemake --cores 16 --rerun-incomplete --keep-going --configfile config/snakemake.yaml
```

`-n -p` prints the planned jobs and their commands without executing; `--rerun-incomplete` re-runs any jobs a previous interrupted run left half-finished; `--keep-going` lets independent jobs continue when one fails.

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
  step5b_min_samples/    # low-sample intervals dropped, if min_samples is set
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
- [`trim_regions.py`](scripts/trim_regions.py) — apply a BED mask to a directory of tree sequences and write trimmed outputs with compacted coordinates.
- [`trim_regions_single.py`](scripts/trim_regions_single.py) — apply a BED mask to one tree sequence while preserving the original coordinate system.
- [`trim_samples.py`](scripts/trim_samples.py) — remove individuals genome-wide (`--individuals`) or over BED intervals (`--remove`). See [NOTES.md](NOTES.md) for the exact sample-ID matching rules.
- [`filter_min_samples.py`](scripts/filter_min_samples.py) — drop intervals whose local trees retain fewer than `--min-samples` non-isolated sample nodes, via `delete_intervals` (coordinates preserved); updates `kept_intervals` metadata and writes a diagnostic BED. Optional pipeline step 5b, driven by `min_samples`.
- [`validation_plots_from_ts.py`](scripts/validation_plots_from_ts.py) — SINGER-style QC plots (mutational load, diversity, Tajima's D, folded/unfolded SFS) across TS replicates; optional observed-vs-simulated overlays.
- [`merge_treefiles_by_replicate.py`](scripts/merge_treefiles_by_replicate.py) — concatenate chromosome-specific tree-sequence files by replicate; embedded mutation-rate ratemaps are merged and carried forward.
- [`export_vcf.py`](scripts/export_vcf.py) — export a `.vcf`/`.vcf.gz` from a (filtered) tree sequence: variable sites only, ploidy-aware genotypes, and samples pruned by `trim_samples` written as missing (`.`) via `isolated_as_missing`. Driven by `emit_vcf`; see [VCF export](#vcf-export).
- [`pipeline_summary.py`](scripts/pipeline_summary.py) — self-contained HTML report of genome retention, per-individual outlier counts, and embedded validation plots. Requires `--filtered-ts` (the final per-chromosome tree sequences, in the `<chrom>/<rep>.<suffix>` layout produced by steps 5/5b), which it loads to measure retained sequence directly from each ARG.

### ⚠️ Warning: summary statistics after sample pruning

Diversity (π), Tajima's D, and SFS in this pipeline are computed with [tskit](https://tskit.dev/)'s built-in methods (`ts.diversity`, `ts.Tajimas_D`, `ts.allele_frequency_spectrum`). These normalize by a **constant** `n · (n − 1) / 2` based on the sample set passed in, so when the sample size varies across regions — as it does after `trim_samples.py` removes individuals over BED intervals — the per-window statistics are **not correctly normalized** for the locally retained sample count. Treat post-pruning stats with caution, and prefer comparisons on replicates where sample membership is uniform across the genome.

## Auxiliary scripts

Scripts not called by the Snakemake pipeline.

- [`coalescence_ne_plots_from_ts.py`](scripts/coalescence_ne_plots_from_ts.py) — pair-coalescence and Ne plots from TS replicates. Choose the time grid with one of `--time-bins-file` (explicit edges), `--num-quantiles N` (equal-coalescence-event bins derived from `pair_coalescence_quantiles` averaged over post-burnin replicates), or `--num-bins N` (uniform log-spaced bins across the observed coalescence-time range). Note that `--num-bins` meant equal-coalescence-event bins up to v1.8; that mode is now `--num-quantiles`, so pre-v1.9 commands need renaming to reproduce their old time grid. Optional Demes-based coalescent simulations (`--sim N`) produce window-stat and SFS TSVs for observed-vs-sim density plots in `validation_plots_from_ts.py`.

  The environment pins [nspope/tskit commit `73d8cd9`](https://github.com/nspope/tskit/commit/73d8cd922482475020ae01180cae95bf5abbf067), whose pair-coalescence quantiles, counts, and rates normalize over locally non-missing pair-spans. Consequently, sample-isolated intervals are handled natively for global, spatial-window, and multiple-sample-set calculations; no script-level scalar correction is applied.
- [`compare_trees_html.py`](scripts/compare_trees_html.py) — render one tree index from each of two tree sequences side-by-side into a single HTML file.
- [`locate_tree.py`](scripts/locate_tree.py) — find the local tree at a `(chromosome, position)` in a merged tree sequence, using the `chrom_offsets` metadata to map the within-chromosome coordinate to the concatenated axis. See [Locating a tree by chromosome and position](#locating-a-tree-by-chromosome-and-position).
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

None of this would be possible without the patient help and advice of Nate Pope. Any errors, bad code, or poor interpretations, however, are my responsibility alone. This repo also uses code from Nate Pope's [singer-snakemake](https://github.com/nspope/singer-snakemake).
