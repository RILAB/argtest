# ARG Tree Sequence Utilities and Validation Plotting

Standalone scripts for post-processing, QC, and visualization of ARG tree sequences (`.ts`, `.trees`, `.tsz`).

## Install

```bash
conda env create -f environment.yml
conda activate argtest
```

Core dependencies are in `environment.yml` (`numpy`, `matplotlib`, `tskit`, `tszip`, `msprime`).

## Scripts

### `scripts/mutload_summary.py`
Builds an HTML summary of per-individual mutational load and writes outlier windows as BED when windowing is enabled.

Inputs:
- one tree sequence file

Key outputs:
- `results/<name>.html`
- `results/<ts_stem>_outliers.bed` (if `--window-size` is used)
- `logs/<name>.log`

Example:
```bash
python scripts/mutload_summary.py example_data/maize.tsz --window-size 50000 --cutoff 0.25 --out mutload.html
```

Main options:
- `--window-size`
- `--cutoff`
- `--out`
- `--suffix-to-strip`

### `scripts/trim_samples.py`
Removes selected individuals either genome-wide (`--individuals`) or over BED intervals (`--remove`).

Inputs:
- one tree sequence file
- optional BED(s) with per-individual intervals

Key output:
- trimmed tree sequence (`--out`, or `results/<ts_stem>_trimmed.tsz`)

Example:
```bash
python scripts/trim_samples.py example_data/maize.tsz --individuals B73,Mo17 --out results/maize_trimmed.tsz
```

Main options:
- `--individuals`
- `--remove`
- `--out`
- `--suffix-to-strip`

### `scripts/trim_regions.py`
Applies a shared accessibility-based mask across a directory of tree sequences by:
1. inferring mutation map (`*.mut_rate.p`) from TS names,
2. keeping windows with accessible bp >= `--cutoff-bp`,
3. collapsing masked intervals and writing renamed output TS files.

Inputs:
- directory of tree sequences
- inferred mutation-rate map file(s)

Key outputs:
- collapsed tree sequences in output directory
- one summary log (`collapse_log.txt` by default)

Example:
```bash
python scripts/trim_regions.py \
  --ts-dir /path/to/trees \
  --window-size 50000 \
  --cutoff-bp 2500 \
  --pattern "*.tsz"
```

Main options:
- `--ts-dir`
- `--window-size`
- `--cutoff-bp`
- `--out-dir`
- `--pattern`
- `--log`

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

## Repository notes

- Generated `logs/` and `results/` are git-ignored.
- `.DS_Store` is git-ignored.
