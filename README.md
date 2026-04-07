# ARG Tree Sequence Utilities and Validation Plotting

Standalone scripts for post-processing, QC, and visualization of ARG tree sequences (`.ts`, `.trees`, `.tsz`).

## Install

```bash
conda env create -f environment.yml
conda activate argtest
```

Core dependencies are in `environment.yml` (`numpy`, `matplotlib`, `tskit`, `tszip`, `msprime`).

## Suggested Workflow

One reasonable post-processing workflow for ARG tree sequences in this repo is:

1. **Find low rec regions** Because regions of low recombination are more affected by linked selection, for analyses assuming the neutrality of the ARG it may be a good idea to remove low recombination regions ahead of time. Find windows in the genetic map in the bottom `X` percentile of `cM/Mb` using [scripts/hapmap_low_rec_mask.py](/home/jri/src/argtest/scripts/hapmap_low_rec_mask.py). This turns a HapMap-style recombination map plus a `.fai` into per-chromosome BED masks for very low-recombination regions.
2. **Find regions of poor alignment** Find windows of `size` kb where more than `X`% of bp are masked using [scripts/find_low_access_regions.py](/home/jri/src/argtest/scripts/find_low_access_regions.py). This inspects the inferred mutation map for a tree sequence and writes low-accessibility windows to a BED file.
3. **Find windows with aberrant mutational load** All samples in a tree should have the same number of derived mutations, since all have the same distance from root. In windows of `number` SNPs, identify individuals with `X`% more or fewer derived mutations than the window median using [scripts/mutload_summary.py](/home/jri/src/argtest/scripts/mutload_summary.py). This reports per-window outliers, can write mutation-masked BED intervals when too many individuals are outliers, and highlights outlier individuals in red in the HTML summary.
4. **Remove affected regions** For each chromosome, combine the BED files from steps 1-3 (<chr>.low_rec.mask.bed, low_access.ws<window>.accbp<cutoff>.bed, and <ts_stem>_mutation_masked.bed), then remove those genomic regions from a directory of tree sequences with [scripts/trim_regions.py](/home/jri/src/argtest/scripts/trim_regions.py). This script applies a supplied BED mask and writes trimmed tree sequences.
5. **Remove affected samples** In many cases, only a few samples within a window will be problematic. They could have evidence of introgression (identified using e.g. [TRACE](https://github.com/YulinZhang9806/trace)) or odd patterns of derived mutations (see step 3). Using a bedfile specifying regions and individuals, prune those individuals from the trees with [scripts/trim_samples.py](/home/jri/src/argtest/scripts/trim_samples.py). 
6. **Validate** Run the validation plots with [scripts/validation_plots_from_ts.py](/home/jri/src/argtest/scripts/validation_plots_from_ts.py) to get a sense of the cleaned ARG. This gives a compact set of QC plots for mutational load, diversity, Tajima's D, and related summaries. Compare these to the same plots run on the original treesequences.
7. If satisfied, merge chromosomes for each replicate for downstream analysis using [scripts/merge_treefiles_by_replicate.py](/home/jri/src/argtest/scripts/merge_treefiles_by_replicate.py). This concatenates chromosome-specific tree files into one combined tree sequence per replicate.

## Scripts

### `scripts/hapmap_low_rec_mask.py`
Builds one BED per chromosome marking the bottom `--rec-fraction` of HapMap recombination-rate intervals, ranked by `Rate(cM/Mb)` separately within each chromosome.

The HapMap input should follow the format described in the msprime docs: <https://tskit.dev/msprime/docs/stable/api.html#msprime.RateMap.read_hapmap>.

Behavior:
- requires `--hapmap` and `--fai`
- uses the first row's `Rate(cM/Mb)` from position `0` up to that row's `Position(bp)`
- uses each row's `Rate(cM/Mb)` from its position up to the next position on the same chromosome
- uses the last row's `Rate(cM/Mb)` from its position to the chromosome end from the `.fai`
- writes `results` as one BED per chromosome named `<chr>.low_rec.mask.bed` in `--out-dir` (or the current directory by default)

Example:
```bash
python scripts/hapmap_low_rec_mask.py \
  --hapmap maize.hapmap.tsv \
  --fai maize.fa.fai \
  --rec-fraction 0.1 \
  --out-dir low_rec_masks
```

Main options:
- `--hapmap`
- `--fai`
- `--rec-fraction`
- `--out-dir`

### `scripts/find_low_access_regions.py`
Finds low-accessibility windows from the inferred mutation map for a representative tree sequence and writes those windows as a BED file.

Behavior:
- infers the mutation map (`*.mut_rate.p`) from the first tree sequence in `--ts-dir`
- computes accessible bp in windows of `--window-size`
- writes windows with accessible bp below `--cutoff-bp` to a BED file
- writes to `--out` or to `<ts-dir>/low_access.ws<window>.accbp<cutoff>.bed` by default

Inputs:
- directory of tree sequence files
- inferred mutation-rate map file(s)

Key outputs:
- one BED file of low-accessibility windows

Example:
```bash
python scripts/find_low_access_regions.py \
  --ts-dir /path/to/trees \
  --window-size 50000 \
  --cutoff-bp 2500
```

Main options:
- `--ts-dir`
- `--window-size`
- `--cutoff-bp`
- `--pattern`
- `--out`

### `scripts/mutload_summary.py`
Builds an HTML summary of per-individual mutational load and writes outlier windows as BED when windowing is enabled. When `--window-size` or `--snp-window` is used, each individual is compared to that window's median mutational load and is marked as an outlier if its load is greater than `(1 + cutoff) * median` or less than `(1 - cutoff) * median`; windows with median zero are skipped. In the per-window plots, outlier individuals are shown with red bars. If `--fraction` is provided, windows with an outlier fraction greater than that threshold are written to `results/<ts_stem>_mutation_masked.bed` and are excluded from `results/<ts_stem>_outliers.bed`.

Inputs:
- one tree sequence file

Key outputs:
- `results/<name>.html`
- `results/<ts_stem>_outliers.bed` (if `--window-size` or `--snp-window` is used; excludes windows written to the mutation-masked BED when `--fraction` is set)
- `results/<ts_stem>_mutation_masked.bed` (if `--fraction` is used)
- `logs/<name>.log`

Example:
```bash
python scripts/mutload_summary.py example_data/maize.tsz --window-size 50000 --cutoff 0.25 --out mutload.html
```

Example with fixed-variant windows:
```bash
python scripts/mutload_summary.py example_data/maize.tsz --snp-window 500 --cutoff 0.25 --out mutload_snp.html
```

Main options:
- `--window-size`
- `--snp-window`
- `--cutoff`
- `--fraction`
- `--out`
- `--suffix-to-strip`

### `scripts/trim_regions.py`
Applies a BED mask to a directory of tree sequences and writes trimmed tree files with compacted coordinates.

Inputs:
- directory of tree sequences
- one BED file of regions to remove

Key outputs:
- trimmed tree sequences in output directory
- one summary log (`trim_regions_log.txt` by default)

Example:
```bash
python scripts/trim_regions.py \
  --ts-dir /path/to/trees \
  --remove low_access.bed \
  --pattern "*.tsz" \
  --out-dir /path/to/trimmed
```

Main options:
- `--ts-dir`
- `--remove`
- `--out-dir`
- `--pattern`
- `--log`

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

Main options:
- `--individuals`
- `--remove`
- `--out`
- `--suffix-to-strip`

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
- `--sim` can read the simulation TSV from `scripts/coalescence_ne_plots_from_ts.py --sim` and compare observed vs simulated windowed `pi` and Tajima's D using overlapping semi-transparent density plots

Example:
```bash
python scripts/validation_plots_from_ts.py \
  --ts-dir ~/crud/collapsed \
  --pattern "*.tsz" \
  --window-size 100000 \
  --mutation-rate 3.3e-8 \
  --burnin-frac 0.5 \
  --out-dir results/validation_plots
```

Example with simulation comparison:
```bash
python scripts/validation_plots_from_ts.py \
  --ts-dir ~/crud/collapsed \
  --pattern "*.tsz" \
  --window-size 100000 \
  --mutation-rate 3.3e-8 \
  --burnin-frac 0.5 \
  --sim results/coalescence_ne_plots/sim-window-stats.tsv \
  --out-dir results/validation_plots
```

Main options:
- `--ts-dir`
- `--pattern`
- `--mutation-rate`
- `--burnin-frac`
- `--window-size`
- `--folded`
- `--sim`
- `--out-dir`
- `--prefix`

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

Main options:
- `--ts-dir`
- `--out-dir`
- `--pattern`
- `--out-suffix`

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
  --ts-dir ~/crud/collapsed \
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
  --ts-dir ~/crud/collapsed \
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

Main options:
- `--ts-dir`
- `--pattern`
- `--time-bins-file`
- `--burnin-frac`
- `--tail-cutoff`
- `--time-adjust`
- `--year`
- `--sim`
- `--mu`
- `--sim-outfile`
- `--log-rates`
- `--out-dir`
- `--prefix`

### `scripts/compare_trees_html.py`
Renders one tree index from each of two tree sequences into a single HTML comparison.

Example:
```bash
python scripts/compare_trees_html.py a.tsz 9 b.tsz 9 --out tree_compare.html
```

### `scripts/trees_gallery_html.py`
Renders all trees from two tree sequences as a top/bottom-row scrollable HTML gallery.

Example:
```bash
python scripts/trees_gallery_html.py a.tsz b.tsz --out trees_gallery.html
```

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
  - `trim_regions.py`: default `--out-dir` is `<ts-dir>/trimmed` and default log is `<out-dir>/trim_regions_log.txt`.
  - `mutload_summary.py` and several plotting scripts write HTML/PNG files into `results/` by default; many have `--out` or `--out-dir` flags to override this.
- **Logging & errors:** scripts write summary logs into `logs/` or into the chosen `--out-dir` (see each script's `--log`/`--out-dir` options). Common failure modes include missing `tszip` for `.tsz` I/O, mismatched sequence lengths across chromosome files (checked by `trim_regions.py`), and invalid BED line formats (the loader will raise `ValueError` when a BED line has fewer than 3 columns). When a script prints an `ERROR:` or raises an exception, check the corresponding `logs/` or `<out-dir>/` log file for the detailed run record.

## Sample ID matching (trim_samples.py)

- `trim_samples.py` matches sample/individual IDs exactly against the tree sequence internal individual names as produced by `scripts/argtest_common.py`: `get_individual_name()` (this prefers `individual.metadata['id']` when present, otherwise a synthetic `ind<id>` name is used). Matching is exact and case-sensitive. The `--suffix-to-strip` option (default `_anchorwave`) is applied by name lookup and is removed via simple string replacement before matching; provide names that match the post-stripping individual names.

## Repository notes

- Generated `logs/` and `results/` are git-ignored.
- `.DS_Store` is git-ignored.
