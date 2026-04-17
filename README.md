# ARG Tree Sequence Utilities and Validation Plotting

Standalone scripts for post-processing, QC, and visualization of ARG tree sequences (`.ts`, `.trees`, `.tsz`).

## Contents

- [Install](#install)
- [Suggested Workflow](#suggested-workflow)
- [Snakemake Pipeline](#snakemake-pipeline)
- [Scripts](#scripts): `hapmap_low_rec_mask` · `find_low_access_regions` · `mutload_summary` · `mutload_masks` · `combine_remove_masks` · `trim_regions` · `trim_regions_single` · `trim_samples` · `validation_plots_from_ts` · `merge_treefiles_by_replicate`
- [Auxiliary scripts](#auxiliary-scripts): `run_steps1_5_and_concat` · `coalescence_ne_plots_from_ts` · `compare_trees_html` · `trees_gallery_html`
- [Shared module](#shared-module)
- [Inputs, formats, defaults & logs](#inputs-formats-defaults--logs)
- [Sample ID matching](#sample-id-matching-trim_samplespy)
- [Repository notes](#repository-notes)

## Install

```bash
conda env create -f environment.yml
conda activate argtest
```

Core dependencies are in `environment.yml` (`numpy`, `matplotlib`, `tskit`, `tszip`, `msprime`).

## Suggested Workflow

One reasonable post-processing workflow for ARG tree sequences in this repo is:

1. **Find low rec regions** Because regions of low recombination are more affected by linked selection, for analyses assuming the neutrality of the ARG it may be a good idea to remove low recombination regions ahead of time. Find windows in the genetic map in the bottom `X` percentile of `cM/Mb` using [scripts/hapmap_low_rec_mask.py](scripts/hapmap_low_rec_mask.py). This turns a HapMap-style recombination map plus a `.fai` into per-chromosome BED masks for very low-recombination regions.
2. **Find regions of poor alignment** Find windows of `size` kb where more than `X`% of bp are masked using [scripts/find_low_access_regions.py](scripts/find_low_access_regions.py). This inspects the inferred mutation map for a tree sequence and writes low-accessibility windows to a BED file.
3. **Find windows with aberrant mutational load** All samples in a tree should have the same number of derived mutations, since all have the same distance from root. In windows of `number` SNPs, identify individuals with `X`% more or fewer derived mutations than the window median. Use [scripts/mutload_summary.py](scripts/mutload_summary.py) for an interactive HTML summary (ASCII bar charts with outlier individuals highlighted in red, plus a lineage table). Use [scripts/mutload_masks.py](scripts/mutload_masks.py) to write the outlier and mutation-masked BED files needed for the pipeline; this is what the Snakemake workflow calls.
4. **Remove affected regions** For each chromosome, combine the BED files from steps 1-3 (<chr>.low_rec.mask.bed, low_access.ws<window>.accbp<cutoff>.bed, and <ts_stem>_mutation_masked.bed), then remove those genomic regions from a directory of tree sequences with [scripts/trim_regions.py](scripts/trim_regions.py). This script applies a supplied BED mask and writes trimmed tree sequences.
5. **Remove affected samples** In many cases, only a few samples within a window will be problematic. They could have evidence of introgression (identified using e.g. [TRACE](https://github.com/YulinZhang9806/trace)) or odd patterns of derived mutations (see step 3). Using a bedfile specifying regions and individuals, prune those individuals from the trees with [scripts/trim_samples.py](scripts/trim_samples.py).
6. **Validate** Run the validation plots with [scripts/validation_plots_from_ts.py](scripts/validation_plots_from_ts.py) to get a sense of the cleaned ARG. This gives a compact set of QC plots for mutational load, diversity, Tajima's D, and related summaries. Compare these to the same plots run on the original treesequences.
7. If satisfied, merge chromosomes for each replicate for downstream analysis using [scripts/merge_treefiles_by_replicate.py](scripts/merge_treefiles_by_replicate.py). This concatenates chromosome-specific tree files into one combined tree sequence per replicate.

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

Supported treefile suffixes are `.ts`, `.trees`, and `.tsz`. Replicate IDs are taken from the filename stem, so `chr1/1.tsz` is replicate `1` for chromosome `chr1`.

### Required Inputs

The Snakemake config expects these keys in a YAML file such as [config/snakemake.example.yaml](config/snakemake.example.yaml):

- `root_dir`: path to the chromosome-subdirectory root
- `hapmap`: HapMap recombination map used for step 1
- `fai`: FASTA index used for chromosome lengths
- `tree_pattern`: glob for treefiles within each chromosome directory, for example `"*.tsz"`
- `rec_fraction`: fraction of recombination-rate intervals (ranked by `Rate(cM/Mb)`, lowest first) to include in the low-recombination mask; e.g. `0.1` masks the bottom 10 % of intervals
- `low_access_window_size`: window size in bp for step 2
- `low_access_cutoff_bp`: minimum accessible bp per window for step 2
- exactly one of `mutload_window_size` or `mutload_snp_window` for step 3
- `mutload_cutoff`: outlier cutoff fraction for step 3
- `mutload_fraction`: optional fraction threshold for writing mutation-masked BED rows in step 3
- `suffix_to_strip`: suffix removed from sample IDs before matching in step 3 and step 5
- `allow_missing_replicates`: set to `true` if you want to concatenate partial replicate sets
- `base_name`: prefix used for merged outputs
- `out_dir`: output root for Snakemake products

Mutation-rate maps are inferred from the treefile location using the same logic as the scripts, so the usual `*.mut_rate.p` files should be available near each chromosome directory.

### Workflow Steps

The Snakemake workflow runs the same logical steps as the manual workflow:

1. Build per-chromosome low-recombination BED masks
2. Build per-chromosome low-accessibility BED masks
3. Build per-chromosome, per-replicate mutational-load outlier BEDs and optional mutation-masked BEDs
4. Combine the BED masks and trim affected regions from each treefile
5. Trim the affected samples from each trimmed treefile
6. Concatenate chromosomes for each replicate into one combined tree sequence

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
cp config/snakemake.example.yaml config/snakemake.yaml
```

#### Walk-through using bundled example data

Use the committed example dataset at [example_data/sim_2chr_5rep_clean](example_data/sim_2chr_5rep_clean):

```yaml
root_dir: example_data/sim_2chr_5rep_clean/trees
hapmap: example_data/sim_2chr_5rep_clean/hapmap/all.hapmap.tsv
fai: example_data/sim_2chr_5rep_clean/sim.fai
tree_pattern: "*.trees"
rec_fraction: 0.1
low_access_window_size: 10000
low_access_cutoff_bp: 9000
mutload_window_size: 10000
mutload_cutoff: 0.25
mutload_fraction: 0.8
suffix_to_strip: ""
allow_missing_replicates: false
base_name: "sim2chr"
out_dir: "/tmp/argtest-snakemake-example-out"
```

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
  logs/
```

Intermediate filenames include both chromosome and replicate information so they stay unique across the full workflow.

For a more detailed walkthrough of the Snakemake inputs and outputs, see [snakemake.md](snakemake.md).

## Scripts

### `scripts/hapmap_low_rec_mask.py`
Builds one BED per chromosome marking the bottom `--rec-fraction` of HapMap recombination-rate intervals, ranked by `Rate(cM/Mb)` separately within each chromosome.

The HapMap input should follow the format described in the msprime docs: <https://tskit.dev/msprime/docs/stable/api.html#msprime.RateMap.read_hapmap>.

Behavior:
- requires `--hapmap` and `--fai`
- uses the first row's `Rate(cM/Mb)` from position `0` up to that row's `Position(bp)`
- uses each row's `Rate(cM/Mb)` from its position up to the next position on the same chromosome
- uses the last row's `Rate(cM/Mb)` from its position to the chromosome end from the `.fai`
- writes one BED per chromosome named `<chr>.low_rec.mask.bed` to `--out-dir` (default: current directory)

Example:
```bash
python scripts/hapmap_low_rec_mask.py \
  --hapmap maize.hapmap.tsv \
  --fai maize.fa.fai \
  --rec-fraction 0.1 \
  --out-dir low_rec_masks
```

Options:

| Flag | Description |
|------|-------------|
| `--hapmap PATH` | *(required)* HapMap recombination map file. |
| `--fai PATH` | *(required)* FASTA index (`.fai`) supplying chromosome lengths. |
| `--rec-fraction FLOAT` | *(required)* Fraction of intervals (ranked by `Rate(cM/Mb)`) to mark as low-recombination, per chromosome. E.g. `0.1` marks the bottom 10 %. |
| `--out-dir PATH` | Directory for output BED files (default: current directory). |
| `--chrom NAME` | If given, process only this chromosome; otherwise all chromosomes in the HapMap are processed. |
| `--log PATH` | Log file path (default: `<out-dir>/logs/hapmap_low_rec_mask.log`). |

### `scripts/find_low_access_regions.py`
Finds low-accessibility windows from the inferred mutation map for a single tree sequence file and writes those windows as a BED file.

Behavior:
- infers the mutation map (`*.mut_rate.p`) from the path of the `--ts` file
- computes accessible bp in windows of `--window-size`
- writes windows with accessible bp below `--cutoff-bp` to a BED file
- writes to `--out` or to `<ts-parent>/low_access.ws<window>.accbp<cutoff>.bed` by default

Inputs:
- one tree sequence file
- inferred mutation-rate map file (must be co-located with the tree sequence)

Key outputs:
- one BED file of low-accessibility windows

Example:
```bash
python scripts/find_low_access_regions.py \
  --ts /path/to/trees/chr1/1.tsz \
  --window-size 50000 \
  --cutoff-bp 2500
```

Options:

| Flag | Description |
|------|-------------|
| `--ts PATH` | *(required)* Tree sequence file (`.tsz`, `.ts`, or `.trees`). The mutation-rate map is inferred from this path. |
| `--window-size FLOAT` | *(required)* Window size in bp for accessibility calculation. |
| `--cutoff-bp FLOAT` | *(required)* Windows with fewer accessible bp than this value are written to the BED output. |
| `--out PATH` | Output BED path (default: `<ts-parent>/low_access.ws<window>.accbp<cutoff>.bed`). |
| `--log PATH` | Log file path (default: same base as `--out` with `.log` extension). |

### `scripts/mutload_summary.py`
Builds a self-contained HTML summary of per-individual mutational load. When `--window-size` or `--snp-window` is used, each individual is compared to that window's median mutational load and is marked as an outlier if its load is greater than `(1 + cutoff) * median` or less than `(1 - cutoff) * median`; windows with median zero are skipped. The HTML contains ASCII bar charts (outlier individuals shown in red), a per-individual outlier-window-count histogram, and a lineage summary table. This script is intended for interactive review; for pipeline BED output see `mutload_masks.py`.

Inputs:
- one tree sequence file

Key outputs:
- `results/<name>.html`
- `logs/<name>.log`

Example:
```bash
python scripts/mutload_summary.py example_data/maize.tsz --window-size 50000 --cutoff 0.25 --out mutload.html
```

Example with fixed-variant windows:
```bash
python scripts/mutload_summary.py example_data/maize.tsz --snp-window 500 --cutoff 0.25 --out mutload_snp.html
```

Options:

| Flag | Description |
|------|-------------|
| `ts` | *(positional, required)* Tree sequence file (`.ts`, `.trees`, or `.tsz`). |
| `--window-size FLOAT` | Window size in bp. Mutually exclusive with `--snp-window`. |
| `--snp-window INT` | Window size as a fixed number of variants. Mutually exclusive with `--window-size`. |
| `--cutoff FLOAT` | Outlier cutoff as a fraction of the window median (default: `0.25`). An individual is an outlier if its load is outside `[(1−cutoff)×median, (1+cutoff)×median]`. |
| `--out PATH` | Output filename (default: `mutational_load_summary.html`). Only the filename part is used; the file is always written to `<repo-root>/results/<filename>`. |
| `--suffix-to-strip STR` | Suffix stripped from sample names before display (default: `""`). |

### `scripts/mutload_masks.py`
Writes outlier and mutation-masked BED files for one tree sequence. This is the script called by the Snakemake pipeline for step 3; use `mutload_summary.py` for a human-readable HTML report.

Behavior:
- identifies per-window outlier individuals using the same median-based cutoff as `mutload_summary.py`
- writes outlier windows (with outlier sample IDs and load values) to `--outlier-bed`
- if `--fraction` is set, windows where the outlier fraction exceeds that threshold are written to `--masked-bed` instead and excluded from `--outlier-bed`; otherwise `--masked-bed` is created empty
- chromosome label written to BED column 1 is set by `--chrom`

Inputs:
- one tree sequence file (`--ts`)

Key outputs:
- `--outlier-bed`: per-window outlier intervals with sample IDs
- `--masked-bed`: windows with too many outliers (empty if `--fraction` is not set)
- log file (defaults to `<outlier-bed-parent>/logs/<chrom>.<ts_stem>.mutload.log`)

Example:
```bash
python scripts/mutload_masks.py \
  --ts /path/to/trees/chr1/1.tsz \
  --chrom chr1 \
  --window-size 50000 \
  --cutoff 0.25 \
  --fraction 0.2 \
  --outlier-bed step3/chr1/1.outliers.bed \
  --masked-bed step3/chr1/1.mutation_masked.bed
```

Options:

| Flag | Description |
|------|-------------|
| `--ts PATH` | *(required)* Input tree sequence file. |
| `--chrom NAME` | *(required)* Chromosome label written in BED column 1. |
| `--window-size FLOAT` | *(required, or use `--snp-window`)* Window size in bp. |
| `--snp-window INT` | *(required, or use `--window-size`)* Window size as a fixed number of variants. |
| `--cutoff FLOAT` | Outlier cutoff as a fraction of the window median (default: `0.25`). |
| `--fraction FLOAT` | If provided, windows where the outlier fraction exceeds this threshold are written to `--masked-bed` and excluded from `--outlier-bed`; if omitted, `--masked-bed` is created empty. |
| `--outlier-bed PATH` | *(required)* Output BED of per-window outlier intervals, with sample IDs and load values in extra columns. |
| `--masked-bed PATH` | *(required)* Output BED of windows with too many outliers (empty when `--fraction` is not set). |
| `--suffix-to-strip STR` | Suffix removed from sample names before grouping to individuals (default: `""`). |
| `--log PATH` | Log file path (default: `<outlier-bed-parent>/logs/<chrom>.<ts_stem>.mutload.log`). |

### `scripts/combine_remove_masks.py`
Merges multiple BED mask files (low-rec, low-access, mutation-masked) into a single combined BED. Missing input files are silently skipped. This is called by the Snakemake pipeline between steps 3 and 4.

Example:
```bash
python scripts/combine_remove_masks.py \
  --chrom chr1 \
  --out step4_masks/chr1.combined.bed \
  --inputs step1/chr1.low_rec.mask.bed step2/low_access.ws50000.accbp2500.bed step3/chr1/1_mutation_masked.bed
```

Options:

| Flag | Description |
|------|-------------|
| `--chrom NAME` | *(required)* Chromosome name written in BED column 1 of the output. |
| `--out PATH` | *(required)* Output combined BED file. |
| `--inputs PATH [PATH ...]` | *(required)* One or more input BED files to merge. Missing files are ignored. |
| `--log PATH` | Log file path (default: `<out-parent>/logs/<chrom>.combine_masks.log`). |

### `scripts/trim_regions.py`
Applies a BED mask to a directory of tree sequences and writes trimmed tree files with compacted coordinates.

Inputs:
- directory of tree sequences
- one BED file of regions to remove

Key outputs:
- trimmed tree sequences in output directory
- one summary log (`<out-dir>/logs/trim_regions.log` by default)

Example:
```bash
python scripts/trim_regions.py \
  --ts-dir /path/to/trees \
  --remove low_access.bed \
  --pattern "*.tsz" \
  --out-dir /path/to/trimmed
```

Options:

| Flag | Description |
|------|-------------|
| `--ts-dir PATH` | *(required)* Directory containing tree sequence files. |
| `--remove PATH` | *(required)* BED file of regions to remove. |
| `--out-dir PATH` | Output directory for trimmed tree files (default: `<ts-dir>/trimmed`). |
| `--pattern GLOB` | Glob pattern to filter input filenames (default: `*`). |
| `--log PATH` | Log file path (default: `<out-dir>/logs/trim_regions.log`). |

### `scripts/trim_regions_single.py`
Applies a BED mask to a **single** tree sequence file and writes a trimmed output. This is the per-file variant called by the Snakemake pipeline (step 4); use `trim_regions.py` to batch-process a whole directory.

Example:
```bash
python scripts/trim_regions_single.py \
  --ts /path/to/trees/chr1/1.tsz \
  --remove combined_masks.bed \
  --out trimmed/chr1/1.tsz
```

Options:

| Flag | Description |
|------|-------------|
| `--ts PATH` | *(required)* Input tree sequence file. |
| `--remove PATH` | *(required)* BED file of regions to remove. |
| `--out PATH` | *(required)* Output tree sequence file. |
| `--log PATH` | Log file path (default: `<out-parent>/logs/<ts_stem>.trim_regions.log`). |

### `scripts/trim_samples.py`
Removes selected individuals either genome-wide (`--individuals`) or over BED intervals (`--remove`).

For `--remove`, you can use one BED for all samples or multiple BED files. Each BED line should contain `chrom  start  end  sample_id` in column 4, and column 4 may also list multiple comma-separated sample IDs to apply the same interval to several individuals. If column 4 is omitted, the BED filename stem is used as the sample name.

Inputs:
- one tree sequence file
- optional BED(s) with per-individual intervals

Key output:
- trimmed tree sequence (`--out`, or `results/<ts_stem>_trimmed.tsz`)

Example:
```bash
python scripts/trim_samples.py example_data/maize.tsz --individuals B73,Mo17 --out results/maize_trimmed.tsz
```

Example BED lines:
```text
chr1    1000    5000    A
chr1    8000    9000    B,C
```

Options:

| Flag | Description |
|------|-------------|
| `ts` | *(positional, required)* Tree sequence file (`.ts`, `.trees`, or `.tsz`). |
| `--individuals IDs` | Comma-separated individual IDs to remove across the entire sequence. |
| `--remove PATH` | BED file of per-individual intervals to remove. Column 4 should contain the sample ID (comma-separated for multiple IDs sharing the same interval); if column 4 is absent the BED filename stem is used. Can be repeated or given as a comma-separated list to supply multiple BED files. |
| `--out PATH` | Output tree sequence path (default: `results/<ts_stem>_trimmed.tsz`). |
| `--suffix-to-strip STR` | Suffix removed from sample names before matching (default: `""`). |
| `--log PATH` | Log file path (default: `<out-parent>/logs/<ts_stem>_trim_samples.log`). |

### `scripts/validation_plots_from_ts.py`
Generates SINGER-style validation/diagnostic plots (excluding coalescence/Ne curves) directly from a set of TS replicates.

Plots produced:
- `mutational-load.png`
- `mutational-load-trace.png`
- `diversity-scatter.png`
- `diversity-skyline.png`
- `diversity-trace.png`
- `tajima-d-scatter.png`
- `tajima-d-skyline.png`
- `tajima-d-trace.png`
- `frequency-spectrum.png`
- optional observed-vs-sim density plots (when `--sim` TSV is provided):
  - `diversity-density-vs-sim.png`
  - `tajima-d-density-vs-sim.png`
- `summary.txt`

Notes:
- branch diversity is scaled by `--mutation-rate` for site-vs-branch comparison
- trace plots are branch-only MCMC outcomes
- for observed-vs-simulated plots, first run `coalescence_ne_plots_from_ts.py --sim N ...` to generate the simulation TSVs, then pass them via `--sim` (window stats) and `--sim-sfs` (SFS); the simulations can be reused across multiple validation runs without re-computing them
- **`*.mut_rate.p` file required for original (untrimmed) tree sequences.** When running on the original SINGER output (before any masking), the script auto-detects the nearest `*.mut_rate.p` pickle file (searched in the ts directory and its parent, same logic as `find_low_access_regions.py`) and uses it to normalize diversity per accessible bp rather than per total window span. If no `*.mut_rate.p` is found, the script falls back to total window span and prints no warning. Trimmed tree sequences (step 4/5 output) carry `kept_intervals` metadata and do not need the file.

Example:
```bash
python scripts/validation_plots_from_ts.py \
  --ts-dir /path/to/treefiles \
  --pattern "*.tsz" \
  --window-size 100000 \
  --mutation-rate 3.3e-8 \
  --burnin-frac 0.5 \
  --out-dir results/validation_plots
```

Example with simulation comparison:
```bash
python scripts/validation_plots_from_ts.py \
  --ts-dir /path/to/treefiles \
  --pattern "*.tsz" \
  --window-size 100000 \
  --mutation-rate 3.3e-8 \
  --burnin-frac 0.5 \
  --sim results/coalescence_ne_plots/sim-window-stats.tsv \
  --out-dir results/validation_plots
```

Options:

| Flag | Description |
|------|-------------|
| `--ts-dir PATH` | *(required)* Directory containing tree sequence files. |
| `--mutation-rate FLOAT` | *(required)* Mutation rate used to scale branch diversity for comparison to site diversity. |
| `--pattern GLOB` | Glob pattern for input trees (default: `*.tsz`). |
| `--window-size FLOAT` | Window size in bp for diversity and Tajima's D plots (default: `50000`). |
| `--burnin-frac FLOAT` | Fraction of initial trees to discard as burn-in when computing posterior means (default: `0.5`). |
| `--folded` | Plot folded SFS (minor-allele frequency) instead of polarised derived-frequency SFS. |
| `--sim PATH` | Optional TSV of simulated window statistics (from `coalescence_ne_plots_from_ts.py --sim`) for observed-vs-simulated density plots. |
| `--sim-sfs PATH` | Optional TSV of simulated site frequency spectra (from `coalescence_ne_plots_from_ts.py`) for an observed-vs-simulated SFS plot. |
| `--compare PATH` | Optional second tree-sequence directory to overlay on all plots for comparison (e.g. pre- vs post-pipeline). Uses the same `--pattern`, `--window-size`, `--burnin-frac`, and `--mutation-rate` as the primary directory. Each dataset is labelled by its directory name. |
| `--out-dir PATH` | Output directory for plots (default: `results/validation_plots`). |
| `--prefix STR` | Optional prefix prepended to output plot filenames. |

### `scripts/merge_treefiles_by_replicate.py`
Merges chromosome-specific tree sequence files by replicate. Input files must be named like `<base>.<chromosome>.<replicate>` with suffix `.tree`, `.trees`, or `.tsz`.

Behavior:
- groups files by `<base>` and `<replicate>`
- sorts chromosomes in natural order within each replicate group
- concatenates all chromosomes for a replicate into one combined tree sequence
- writes outputs named `<base>.combined.<replicate><suffix>`
- writes to `--out-dir` or to `<ts-dir>/combined` by default
- preserves the suffix of the first file in each group unless `--out-suffix` is provided

Inputs:
- directory of chromosome-specific tree sequence files

Key outputs:
- one combined tree sequence per replicate

Example:
```bash
python scripts/merge_treefiles_by_replicate.py \
  --ts-dir /path/to/treefiles
```

Options:

| Flag | Description |
|------|-------------|
| `--ts-dir PATH` | *(required)* Directory containing chromosome-specific tree sequence files. |
| `--layout {auto,flat,nested}` | Input layout. `flat` expects `<base>.<chromosome>.<replicate><suffix>` files directly in `--ts-dir`; `nested` expects `--ts-dir/<chromosome>/<replicate><suffix>`; `auto` detects either (default: `auto`). |
| `--pattern GLOB` | Glob pattern to filter input filenames (default: `*`). |
| `--base-name STR` | Base name for nested-layout outputs (default: the `--ts-dir` directory name). Ignored for flat layout. |
| `--out-dir PATH` | Output directory for merged tree sequences (default: `<ts-dir>/combined`). |
| `--out-suffix SUFFIX` | Output file suffix (`.tree`, `.trees`, or `.tsz`; default: suffix of the first file in each group). |
| `--replicate ID` | If set, only write the merged output for this replicate ID. |

## Auxiliary scripts

Scripts not called by the Snakemake pipeline.

### `scripts/run_steps1_5_and_concat.py`
Convenience script that runs all five pipeline steps (low-rec mask → low-access mask → mutational-load masks → trim regions → trim samples) and the final per-replicate merge in a single invocation without Snakemake. Useful for small runs or testing.

Example:
```bash
python scripts/run_steps1_5_and_concat.py \
  --root /path/to/trees \
  --hapmap maize.hapmap.tsv \
  --fai maize.fa.fai \
  --rec-fraction 0.1 \
  --window-size 50000 \
  --cutoff-bp 2500 \
  --mutload-window-size 50000 \
  --mutload-cutoff 0.25 \
  --pattern "*.tsz" \
  --out-dir results/batch
```

Options:

| Flag | Description |
|------|-------------|
| `--root PATH` | *(required)* Root directory containing one subdirectory per chromosome. |
| `--hapmap PATH` | *(required)* HapMap recombination map for step 1. |
| `--fai PATH` | *(required)* FASTA index for chromosome lengths (step 1). |
| `--rec-fraction FLOAT` | *(required)* Low-recombination fraction for step 1. |
| `--window-size FLOAT` | *(required)* Window size in bp for low-accessibility step 2. |
| `--cutoff-bp FLOAT` | *(required)* Accessibility cutoff in bp for step 2. |
| `--mutload-window-size FLOAT` | *(required, or use `--mutload-snp-window`)* Mutational-load window size in bp for step 3. |
| `--mutload-snp-window INT` | *(required, or use `--mutload-window-size`)* Mutational-load window size in variants for step 3. |
| `--mutload-cutoff FLOAT` | Outlier cutoff fraction for step 3 (default: `0.25`). |
| `--mutload-fraction FLOAT` | If provided, windows with outlier fraction above this threshold are written to the mutation-masked BED in step 3. |
| `--pattern GLOB` | Glob pattern for tree files inside each chromosome directory (default: `*`). |
| `--suffix-to-strip STR` | Suffix stripped from individual IDs before name matching (default: `""`). |
| `--out-dir PATH` | Output directory (default: `<root>/batch_steps1_5`). |
| `--merged-dir PATH` | Directory for per-replicate merged outputs (default: `<out-dir>/combined`). |
| `--base-name STR` | Base name for merged output filenames (default: root directory name). |
| `--out-suffix SUFFIX` | Output suffix for merged files (`.tree`, `.trees`, or `.tsz`; default: suffix of first input file). |
| `--allow-missing-replicates` | Allow concatenation when a replicate is absent in some chromosome directories. |

### `scripts/coalescence_ne_plots_from_ts.py`
Generates pair coalescence and effective population size plots from a set of TS replicates using explicit time bins.

Plots produced:
- `pair-coalescence-pdf.png`
- `pair-coalescence-rates.png`
- `effective-pop-size.png` (`Ne = 1 / (2 * coal_rate)`)
- optional simulation summary table: `sim-window-stats.tsv` (when `--sim > 0`)
- `summary.txt`

Notes:
- time bins come from `--time-bins-file` (explicit bin edges)
- `--time-adjust` rescales plotted x-axis values by dividing time by a factor
- `--year` adds a red dashed vertical marker to the Ne plot
- optional `--sim X` mode converts the inferred piecewise-constant Ne trajectory into a one-deme Demes model and runs `X` coalescent simulations
- simulations are `1 Mb`, use the same sample size as the observed TS, and use `--mu` for both mutation and recombination rates
- simulated nucleotide diversity and Tajima's D are computed in `50 Kb` windows and written as a TSV

Example:
```bash
python scripts/coalescence_ne_plots_from_ts.py \
  --ts-dir /path/to/treefiles \
  --pattern "*.tsz" \
  --time-bins-file /path/to/time_bins.txt \
  --burnin-frac 0.5 \
  --time-adjust 6.19476 \
  --year 534 \
  --log-rates \
  --out-dir results/coalescence_ne_plots
```

Example with simulations:
```bash
python scripts/coalescence_ne_plots_from_ts.py \
  --ts-dir /path/to/treefiles \
  --pattern "*.tsz" \
  --time-bins-file /path/to/time_bins.txt \
  --burnin-frac 0.5 \
  --time-adjust 6.19476 \
  --year 534 \
  --sim 100 \
  --mu 3.3e-8 \
  --sim-outfile results/coalescence_ne_plots/sim-window-stats.tsv \
  --out-dir results/coalescence_ne_plots \
  --log-rates
```

Options:

| Flag | Description |
|------|-------------|
| `--ts-dir PATH` | *(required)* Directory containing tree sequence files. |
| `--time-bins-file PATH` | *(required)* File of explicit time-bin edges (one per line) defining the coalescence time grid. |
| `--pattern GLOB` | Glob pattern for input trees (default: `*.tsz`). |
| `--burnin-frac FLOAT` | Fraction of initial trees to discard as burn-in (default: `0.5`). |
| `--tail-cutoff FLOAT` | Minimum probability mass threshold for pair-coalescence tail trimming (default: `1e-12`). |
| `--time-adjust FLOAT` | Divide plotted time-axis values by this factor (e.g. generations per year) to convert to calendar time (default: `1.0`). |
| `--year FLOAT` | If provided, draw a red dashed vertical marker at this x-position on the Ne plot (e.g. a known bottleneck year). |
| `--log-rates` | Use log y-axis for pair-coalescence-rates and Ne plots. |
| `--sim INT` | Number of 1 Mb coalescent simulations to run under a Demes model built from the inferred Ne trajectory (default: `0` = no simulations). |
| `--mu FLOAT` | Mutation rate for simulations; recombination rate is set equal to this value. Required when `--sim > 0`. |
| `--sim-outfile PATH` | Output TSV of simulated 50 Kb window statistics (default: `<out-dir>/<prefix>sim-window-stats.tsv`). |
| `--sim-sfs-outfile PATH` | Output TSV of simulated site frequency spectra across simulations (default: `<out-dir>/<prefix>sim-sfs.tsv`). |
| `--sim-length FLOAT` | Sequence length in bp for each simulation (default: `1000000`). |
| `--sim-window-size FLOAT` | Window size in bp for simulated diversity/Tajima's D TSV (default: `50000`). |
| `--out-dir PATH` | Output directory for plots (default: `results/coalescence_ne_plots`). |
| `--prefix STR` | Optional prefix prepended to output filenames. |

> **TODO:** Running `coalescence_ne_plots_from_ts.py --sim N` before and after the pipeline (on the raw input tree sequences and on the step 5 output) and comparing the resulting Ne trajectories and simulation TSVs would be a useful formal QC step. Consider making this a standard part of the pipeline alongside `validation_plots_from_ts.py --compare`.

### `scripts/compare_trees_html.py`
Renders one tree index from each of two tree sequences side-by-side into a single HTML file for visual comparison.

Example:
```bash
python scripts/compare_trees_html.py a.tsz 9 b.tsz 9 --out tree_compare.html
```

Options:

| Flag | Description |
|------|-------------|
| `ts_a` | *(positional, required)* First tree sequence file (`.ts`, `.trees`, or `.tsz`). |
| `index_a` | *(positional, required)* 0-based tree index to render from `ts_a`. |
| `ts_b` | *(positional, required)* Second tree sequence file (`.ts`, `.trees`, or `.tsz`). |
| `index_b` | *(positional, required)* 0-based tree index to render from `ts_b`. |
| `--out PATH` | Output HTML file (default: `tree_compare.html`). |

### `scripts/trees_gallery_html.py`
Renders all trees from two tree sequences as a top/bottom-row scrollable HTML gallery, useful for quickly comparing before/after trimming.

Example:
```bash
python scripts/trees_gallery_html.py a.tsz b.tsz --out trees_gallery.html
```

Options:

| Flag | Description |
|------|-------------|
| `ts_top` | *(positional, required)* Tree sequence file rendered in the top row (`.ts`, `.trees`, or `.tsz`). |
| `ts_bottom` | *(positional, required)* Tree sequence file rendered in the bottom row. |
| `--out PATH` | Output HTML file (default: `trees_gallery.html`). |

## Shared module

`scripts/argtest_common.py` contains shared tree-sequence helpers used by multiple scripts:
- TS I/O (`load_ts`, `dump_ts`)
- mutational load/stat helpers
- trimming and masking helpers

Use this module for internal script imports.

## Inputs, formats, defaults & logs

- **Tree-sequence files:** scripts accept `.ts`, `.trees`, and `.tsz` files. Loading/writing `.tsz` requires `tszip` to be installed; scripts will raise a clear error if `tszip` is missing when a `.tsz` is used.
- **BED files:** expected as whitespace-separated lines `chrom  start  end  [name]`. `start` and `end` are numeric (half-open intervals `[start, end)`). If a fourth column `name` is present it may list one or more comma-separated sample IDs; if omitted the BED filename stem is used as the sample name. Lines starting with `#` and blank lines are ignored.
- **HapMap recombination maps:** when required (e.g. `hapmap_low_rec_mask.py`), the script expects the HapMap format used by `msprime.RateMap.read_hapmap`.
- **Glob `--pattern`:** arguments named `--pattern` accept shell-style glob patterns (for example "*.tsz") and are matched against filenames in the supplied directory.
- **Defaults & output locations:** many scripts write to a `results/` directory or to an output directory under the input tree-directory when `--out`/`--out-dir` are not provided. Examples:
  - `trim_samples.py`: default output is `results/<ts_stem>_trimmed.tsz` when `--out` is not given.
  - `trim_regions.py`: default `--out-dir` is `<ts-dir>/trimmed` and default log is `<out-dir>/logs/trim_regions.log`.
  - `mutload_summary.py` writes `results/<name>.html` and `logs/<name>.log`; no BED files are written (use `mutload_masks.py` for BED output).
  - Several plotting scripts write PNG files into `results/` by default; most have `--out` or `--out-dir` flags to override this.
- **Logging & errors:** scripts write summary logs into `logs/` or into the chosen `--out-dir` (see each script's `--log`/`--out-dir` options). Common failure modes include missing `tszip` for `.tsz` I/O, mismatched sequence lengths across chromosome files (checked by `trim_regions.py`), and invalid BED line formats (the loader will raise `ValueError` when a BED line has fewer than 3 columns). When a script prints an `ERROR:` or raises an exception, check the corresponding `logs/` or `<out-dir>/` log file for the detailed run record.

## Sample ID matching (trim_samples.py)

- `trim_samples.py` matches sample/individual IDs exactly against the tree sequence internal individual names as produced by `scripts/argtest_common.py`: `get_individual_name()` (this prefers `individual.metadata['id']` when present, otherwise a synthetic `ind<id>` name is used). Matching is exact and case-sensitive. The `--suffix-to-strip` option (default `""`) is applied by name lookup and is removed via simple string replacement before matching; provide names that match the post-stripping individual names.

## Repository notes

- Generated `logs/` and `results/` are git-ignored.
- `.DS_Store` is git-ignored.
