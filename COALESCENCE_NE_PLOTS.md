# Coalescence and Ne plots

Full option reference for [`scripts/coalescence_ne_plots_from_ts.py`](scripts/coalescence_ne_plots_from_ts.py),
an auxiliary script that is not part of the Snakemake pipeline — see
[README.md](README.md#auxiliary-scripts) for the others.

It reads a directory of tree-sequence replicates (treated as MCMC draws, ordered naturally by trailing replicate number), computes pair-coalescence mass and rates on a shared time grid, and writes plots plus the underlying numbers. Times are in generations throughout.

## Inputs and outputs

| Option | Default | Meaning |
| --- | --- | --- |
| `--ts-dir DIR` | *required* | Directory of tree-sequence replicates (`.tsz`, `.ts`, `.trees`). |
| `--pattern GLOB` | `*.tsz` | Glob matched against filenames in `--ts-dir`. |
| `--out-dir DIR` | `results/coalescence_ne_plots` | Output directory; created if missing. |
| `--prefix STR` | none | Prefix prepended to every output filename. |

## Time grid

Supply exactly one of these three; zero or two or more is an error.

| Option | Default | Meaning |
| --- | --- | --- |
| `--time-bins-file FILE` | — | Explicit bin edges, one per line. |
| `--num-quantiles N` | — | `N` equal-coalescence-event bins (N ≥ 2), from connected-pair coalescence quantiles averaged across post-burnin replicates. This mode was called `--num-bins` up to v1.8. |
| `--num-bins N` | — | `N` bins of equal width in log-time (N ≥ 2), spanning the youngest-to-oldest node time observed across post-burnin replicates. |

With `--num-quantiles` / `--num-bins`, `0` and `inf` are padded onto the grid so tskit's rate calculation accepts it; those two padded intervals are dropped from the plots and the TSV, leaving exactly `N` plotted bins.

## Estimation and plotting

| Option | Default | Meaning |
| --- | --- | --- |
| `--burnin-frac F` | `0.5` | Fraction of leading replicates discarded as burn-in. Only post-burnin replicates are plotted and averaged into the posterior mean. |
| `--tail-cutoff X` | `1e-12` | Minimum probability mass retained when trimming the pair-coalescence tail. |
| `--log-rates` | off | Log y-axis on the pair-coalescence-rate and Ne plots (the time axis is always log-scaled). |
| `--year X` | none | Draw a red dashed vertical marker at generation `X` on the Ne plot. |

## Optional Demes simulations

A Demes model is built from the posterior-mean Ne trajectory and simulated with msprime; the resulting TSVs feed the observed-vs-simulated overlays in `validation_plots_from_ts.py`.

| Option | Default | Meaning |
| --- | --- | --- |
| `--sim N` | `0` | Number of independent simulations (`0` skips simulation entirely). |
| `--mu X` | none | Mutation rate for the simulations; recombination rate is set to the same value. Required when `--sim > 0`. |
| `--sim-length BP` | `1000000` | Sequence length per simulation. |
| `--sim-window-size BP` | `50000` | Window size for the simulated diversity / Tajima's D table. |
| `--sim-outfile FILE` | `<out-dir>/<prefix>sim-window-stats.tsv` | Simulated window-statistics TSV. |
| `--sim-sfs-outfile FILE` | `<out-dir>/<prefix>sim-sfs.tsv` | Simulated site-frequency-spectra TSV. |

## Outputs

Written to `--out-dir`, each name prefixed by `--prefix`:

- `pair-coalescence-pdf.png`, `pair-coalescence-rates.png`, `effective-pop-size.png` — post-burnin replicate trajectories (thin lines) with the posterior mean overlaid.
- `coalescence-ne-estimates.tsv` — the plotted numbers, with columns `series` (`replicate` or `posterior_mean`), `replicate_index`, `tree_file`, `bin_index`, `time_left`, `time_right`, `pair_coalescence_mass`, `pair_coalescence_log_density`, `pair_coalescence_rate`, `effective_population_size`.
- `summary.txt` — the run's resolved settings (inputs, burn-in index, time-window edges, simulation parameters, sample count, sequence-length range) and the paths of everything written.
- `sim-window-stats.tsv`, `sim-sfs.tsv` — only when `--sim > 0`.

## Examples

```bash
# 20 equal-coalescence-event bins, log-scaled rate axes
python scripts/coalescence_ne_plots_from_ts.py \
  --ts-dir <out_dir>/combined --pattern "*.tsz" \
  --num-quantiles 20 --burnin-frac 0.5 --log-rates \
  --out-dir results/coalescence_ne_plots

# same grid, plus 100 simulations under the inferred Ne trajectory
python scripts/coalescence_ne_plots_from_ts.py \
  --ts-dir <out_dir>/combined --num-quantiles 20 \
  --sim 100 --mu 1e-8 --prefix cleaned- \
  --out-dir results/coalescence_ne_plots
```

## Pinned tskit fork

This script requires the pinned `nspope/tskit` fork from `environment.yml`. It probes the behaviour at
startup and exits with recreate-the-environment instructions when stock tskit is installed.

The environment pins [nspope/tskit commit `73d8cd9`](https://github.com/nspope/tskit/commit/73d8cd922482475020ae01180cae95bf5abbf067),
whose pair-coalescence quantiles, counts, and rates normalize over locally non-missing pair-spans.
Consequently, sample-isolated intervals are handled natively for global, spatial-window, and
multiple-sample-set calculations; no script-level scalar correction is applied.
