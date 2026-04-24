# ARG Tree Sequence Utilities and Validation Plotting

Snakemake pipeline for post-processing, QC, and visualization of ARG tree sequences (`.ts`, `.trees`, `.tsz`).

If you use, please cite:

Ross-Ibarra, J. 2026. ARGtest: tools for QC and validation of ancestral recombination graphs. [doi: 10.5281/zenodo.19698118](https://doi.org/10.5281/zenodo.19698118)

## Contents

- [Install](#install)
- [Suggested Workflow](#suggested-workflow)
- [Snakemake Pipeline](#snakemake-pipeline)
- [Scripts](#scripts)
- [Auxiliary scripts](#auxiliary-scripts)
- [Inputs, formats, defaults & logs](#inputs-formats-defaults--logs)
- [Repository notes](#repository-notes)

See [NOTES.md](NOTES.md) for the shared module layout and `trim_samples.py` sample-ID matching rules. Run any script with `--help` for its full argument list.

## Install

```bash
conda env create -f environment.yml
conda activate argtest
```

Core dependencies are in `environment.yml` (`numpy`, `matplotlib`, `tskit`, `tszip`, `msprime`, `snakemake`, `pytest`).

## Suggested Workflow

One reasonable post-processing workflow for ARG tree sequences in this repo is:

1. **Find low rec regions** Because regions of low recombination are more affected by linked selection, for analyses assuming the neutrality of the ARG it may be a good idea to remove low recombination regions ahead of time. Find windows in the genetic map in the bottom `X` percentile of `cM/Mb` using [scripts/hapmap_low_rec_mask.py](scripts/hapmap_low_rec_mask.py). This turns a HapMap-style recombination map plus a `.fai` into per-chromosome BED masks for very low-recombination regions.
2. **Find regions of poor alignment** Find windows of `size` kb where more than `X`% of bp are masked using [scripts/find_low_access_regions.py](scripts/find_low_access_regions.py). This inspects the inferred mutation map for a tree sequence and writes low-accessibility windows to a BED file.
3. **Find windows with aberrant mutational load** All samples in a tree should have the same number of derived mutations, since all have the same distance from root. In windows of `number` SNPs, identify individuals with `X`% more or fewer derived mutations than the window median. Use [scripts/mutload_summary.py](scripts/mutload_summary.py) for an interactive HTML summary (ASCII bar charts with outlier individuals highlighted in red, plus a lineage table). Use [scripts/mutload_masks.py](scripts/mutload_masks.py) to write the outlier and mutation-masked BED files needed for the pipeline; this is what the Snakemake workflow calls.
4. **Remove affected regions** For each chromosome, combine the BED files from steps 1-3 (<chr>.low_rec.mask.bed, low_access.ws<window>.accbp<cutoff>.bed, and <ts_stem>_mutation_masked.bed), then remove those genomic regions from a directory of tree sequences with [scripts/trim_regions.py](scripts/trim_regions.py). This script applies a supplied BED mask and writes trimmed tree sequences.
5. **Remove affected samples** In many cases, only a few samples within a window will be problematic. They could have evidence of introgression (identified using e.g. [TRACE](https://github.com/YulinZhang9806/trace)) or odd patterns of derived mutations (see step 3). Using a bedfile specifying regions and individuals, prune those individuals from the trees with [scripts/trim_samples.py](scripts/trim_samples.py).
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

### Required Inputs

The Snakemake config is in [config/snakemake.yaml](config/snakemake.yaml). Edit it for your project and supply it with `--configfile`. The file has an inline comment for every option.

Required keys:

- `root_dir`: path to the chromosome-subdirectory root
- `hapmap`: HapMap recombination map used for step 1
- `fai`: FASTA index used for chromosome lengths
- `rec_fraction`: fraction of recombination-rate intervals (ranked by `Rate(cM/Mb)`, lowest first) to include in the low-recombination mask; e.g. `0.1` masks the bottom 10 % of intervals, while `0.0` writes empty low-recombination masks
- `low_access_window_size`: window size in bp for step 2
- `low_access_cutoff_bp`: minimum accessible bp per window for step 2
- exactly one of `mutload_window_size` or `mutload_snp_window` for step 3

Optional keys (all have sensible defaults):

- `tree_pattern`: glob for treefiles within each chromosome directory (default: `"*"`), for example `"*.trees"` or `"*.tsz"`
- `mutload_cutoff`: outlier cutoff fraction for step 3 (default: `0.25`)
- `mutload_fraction`: fraction threshold for writing mutation-masked BED rows in step 3
- `suffix_to_strip`: suffix removed from sample IDs before matching in step 3 and step 5 (default: `"_anchorwave"`)
- `allow_missing_replicates`: set to `true` to concatenate partial replicate sets (default: `false`)
- `burnin`: number of leading discovered replicates to discard before concatenation (default: `0`); must be smaller than the number of replicates discovered after applying `tree_pattern`
- `base_name`: prefix used for merged outputs (default: name of `root_dir`)
- `merged_out_suffix`: force a specific output suffix for merged files (`.ts`, `.trees`, `.tsz`); default is to inherit the suffix of the first input
- `out_dir`: output root for Snakemake products (default: `snakemake_out`; tilde is expanded)
- `validation_mutation_rate`: mutation rate used by step 6 validation plots; omit or set to `null` to skip step 6
- `validation_first_chrom_only`: run step 6 only on the first chromosome (default: `true`)
- `validation_sim_branch`: simulate site mutations on each ARG replicate with msprime for a posterior-predictive check instead of scaling branch statistics (default: `false`)

Mutation-rate maps are inferred from the treefile location using the same logic as the scripts, so the usual `*.mut_rate.p` files should be available near each chromosome directory.

### Workflow Steps

The Snakemake workflow runs the same logical steps as the manual workflow:

1. Build per-chromosome low-recombination BED masks
2. Build per-chromosome low-accessibility BED masks
3. Build per-chromosome, per-replicate mutational-load outlier BEDs and optional mutation-masked BEDs
4. Combine the BED masks and trim affected regions from each treefile; the mutation-rate ratemap is embedded in each trimmed treefile's metadata
5. Trim the affected samples from each trimmed treefile
6. Concatenate chromosomes for each replicate into one combined tree sequence; ratemaps are merged across chromosomes and carried forward
7. (optional) Run validation plots on original and cleaned tree sequences and generate a self-contained HTML pipeline summary

The final merged outputs are named:

```text
<base_name>.combined.<replicate>.<suffix>
```

and are written under the configured `out_dir` in a `combined/` directory.

### How To Run

From the repo root:

```bash
module load conda
conda activate argtest
```

Edit `config/snakemake.yaml` to point at your data (at minimum set `root_dir`, `hapmap`, and `fai`), then run.

#### Walk-through using bundled example data

The committed [config/snakemake.yaml](config/snakemake.yaml) already points at the example dataset at [example_data/sim_2chr_5rep_clean](example_data/sim_2chr_5rep_clean) and can be used as-is for a test run. It uses `rec_fraction: 0.0` to write empty low-recombination masks and `burnin: 0` because the bundled example has 5 discovered replicates; `burnin` must always be smaller than the discovered replicate count.

In sandboxed environments where `~/.cache` is read-only, set cache/temp dirs to `/tmp` when running Snakemake:

```bash
XDG_CACHE_HOME=/tmp/argtest-xdg-cache TMPDIR=/tmp/argtest-tmp \
python -m snakemake -n -p --configfile config/snakemake.yaml
```

Run the workflow for real:

```bash
XDG_CACHE_HOME=/tmp/argtest-xdg-cache TMPDIR=/tmp/argtest-tmp \
python -m snakemake --cores 16 --rerun-incomplete --keep-going --configfile config/snakemake.yaml
```

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
  pipeline_summary.html
  logs/
```

Intermediate filenames include both chromosome and replicate information so they stay unique across the full workflow.

## Scripts

Pipeline scripts (called by the Snakefile). Run any with `--help` for arguments, defaults, and examples.

- [`hapmap_low_rec_mask.py`](scripts/hapmap_low_rec_mask.py) — per-chromosome BED of the bottom `--rec-fraction` of HapMap recombination-rate intervals.
- [`find_low_access_regions.py`](scripts/find_low_access_regions.py) — BED of low-accessibility windows, computed from a tree sequence's inferred mutation map.
- [`mutload_summary.py`](scripts/mutload_summary.py) — interactive HTML summary of per-individual mutational-load outliers (ASCII bar charts, outlier histogram, lineage table).
- [`mutload_masks.py`](scripts/mutload_masks.py) — outlier and mutation-masked BED files for one tree sequence (pipeline step 3).
- [`combine_remove_masks.py`](scripts/combine_remove_masks.py) — merge the step 1–3 BED masks into a single combined BED per chromosome.
- [`trim_regions.py`](scripts/trim_regions.py) / [`trim_regions_single.py`](scripts/trim_regions_single.py) — apply a BED mask to a directory (or single file) of tree sequences and write trimmed outputs with compacted coordinates.
- [`trim_samples.py`](scripts/trim_samples.py) — remove individuals genome-wide (`--individuals`) or over BED intervals (`--remove`). See [NOTES.md](NOTES.md) for the exact sample-ID matching rules.
- [`validation_plots_from_ts.py`](scripts/validation_plots_from_ts.py) — SINGER-style QC plots (mutational load, diversity, Tajima's D, folded/unfolded SFS) across TS replicates; optional observed-vs-simulated overlays.
- [`merge_treefiles_by_replicate.py`](scripts/merge_treefiles_by_replicate.py) — concatenate chromosome-specific tree-sequence files by replicate; embedded mutation-rate ratemaps are merged and carried forward.
- [`pipeline_summary.py`](scripts/pipeline_summary.py) — self-contained HTML report of genome retention, per-individual outlier counts, and embedded validation plots.

### ⚠️ Warning: summary statistics after sample pruning

Diversity (π), Tajima's D, and SFS in this pipeline are computed with [tskit](https://tskit.dev/)'s built-in methods (`ts.diversity`, `ts.Tajimas_D`, `ts.allele_frequency_spectrum`). These normalize by a **constant** `n · (n − 1) / 2` based on the sample set passed in, so when the sample size varies across regions — as it does after `trim_samples.py` removes individuals over BED intervals — the per-window statistics are **not correctly normalized** for the locally retained sample count. Treat post-pruning stats with caution, and prefer comparisons on replicates where sample membership is uniform across the genome.

## Auxiliary scripts

Scripts not called by the Snakemake pipeline.

- [`coalescence_ne_plots_from_ts.py`](scripts/coalescence_ne_plots_from_ts.py) — pair-coalescence and Ne plots from TS replicates with explicit time bins; optional Demes-based coalescent simulations (`--sim N`) that produce window-stat and SFS TSVs for observed-vs-sim density plots in `validation_plots_from_ts.py`.
- [`compare_trees_html.py`](scripts/compare_trees_html.py) — render one tree index from each of two tree sequences side-by-side into a single HTML file.
- [`trees_gallery_html.py`](scripts/trees_gallery_html.py) — scrollable HTML gallery of all trees from two tree sequences, useful for quick before/after inspection.

> **TODO:** Running `coalescence_ne_plots_from_ts.py --sim N` before and after the pipeline (on the raw input tree sequences and on the step 5 output) and comparing the resulting Ne trajectories and simulation TSVs would be a useful formal QC step. Consider making this a standard part of the pipeline alongside `validation_plots_from_ts.py --compare`.

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

## Acknowledgements

None of this would be possible without the patient help and advice of Nate Pope. Any errors, bad code, or poor interpretations, however, are my responsbility alone. This repo also uses code from Nate Pope's [singer-snakemake](https://github.com/nspope/singer-snakemake).
