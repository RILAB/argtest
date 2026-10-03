# Scripts

← [Back to README](../README.md)

## Pipeline scripts

Pipeline scripts (called by the Snakefile). Run any with `--help` for arguments, defaults, and examples.

- [`hapmap_low_rec_mask.py`](../scripts/hapmap_low_rec_mask.py) — per-chromosome BED of the bottom `--rec-fraction` of HapMap recombination-rate intervals. Chromosome names are matched symmetrically across the two common conventions: `1`, `chr1`, `chr_1` and `Chr-1` all resolve to each other in either direction, and a pipeline base-name prefix (`combined.1`, `amaranth.16`) is stripped from the query first. An exact match always wins, and a file holding two spellings of the same chromosome is an error rather than a coin flip.
- [`find_low_access_regions.py`](../scripts/find_low_access_regions.py) — BED of low-accessibility windows, computed from a tree sequence's inferred mutation map.
- [`mutload_summary.py`](../scripts/mutload_summary.py) — interactive HTML diagnostic for the mutload step: per-individual residual load after window-level pruning (ASCII bar chart, red highlight on individuals still outside the cutoff band, lineage table of flagged counts).
- [`mutload_masks.py`](../scripts/mutload_masks.py) — outlier and mutation-masked BED files for one tree sequence (pipeline step 3).
- [`combine_remove_masks.py`](../scripts/combine_remove_masks.py) — merge the step 1–3 BED masks into a single combined BED per chromosome.
- [`trim_regions_single.py`](../scripts/trim_regions_single.py) — the pipeline step-4 command; apply a BED mask to one tree sequence while preserving the original coordinate system and embedding the resolved mutation map.
- [`trim_samples.py`](../scripts/trim_samples.py) — remove individuals genome-wide (`--individuals`) or over BED intervals (`--remove`). See [Sample ID matching](#sample-id-matching-trim_samplespy) for the exact sample-ID matching rules.
- [`filter_min_samples.py`](../scripts/filter_min_samples.py) — drop intervals whose local trees retain fewer than `--min-samples` non-isolated sample nodes, via `delete_intervals` (coordinates preserved); updates `kept_intervals` metadata and writes a diagnostic BED. Optional pipeline step 5b, driven by `min_samples`.
- [`validation_plots_from_ts.py`](../scripts/validation_plots_from_ts.py) — SINGER-style QC plots (mutational load, diversity, Tajima's D, folded/unfolded SFS) across TS replicates; optional observed-vs-simulated overlays. Also dumps the plotted values to `windows.tsv` / `samples.tsv` / `sfs.tsv`, one table per plot axis, with a `dataset` column so a `--compare` run keeps both series in the same file.
- [`genomewide_expected_vs_observed.py`](../scripts/genomewide_expected_vs_observed.py) — pool the per-chromosome `windows.tsv` files into one expected-vs-observed scatter per statistic (π and Tajima's D), coloured by each window's length-weighted mean recombination rate. Pass `--all-chroms` to have each figure state how much of the genome it actually covers; it prints a warning and titles the plot `PARTIAL GENOME` when chromosomes are missing. Windows the map does not cover are drawn grey rather than as zero-rate; when the positive rates span more than a decade the colour axis is log-scaled and zero-rate windows are clipped to the smallest positive rate (stated on the colourbar).
- [`merge_treefiles_by_replicate.py`](../scripts/merge_treefiles_by_replicate.py) — concatenate chromosome-specific tree-sequence files by replicate; embedded mutation-rate ratemaps are merged and carried forward.
- [`export_vcf.py`](../scripts/export_vcf.py) — export a `.vcf`/`.vcf.gz` from a (filtered) tree sequence: variable sites only, ploidy-aware genotypes, and samples pruned by `trim_samples` written as missing (`.`) via `isolated_as_missing`. Driven by `emit_vcf`; see [VCF export](outputs.md#vcf-export).
- [`pipeline_summary.py`](../scripts/pipeline_summary.py) — self-contained HTML report of genome retention, per-individual outlier counts, the genome-wide expected-vs-observed panels from step 6b, and the embedded per-chromosome validation plots. Requires `--filtered-ts` (the final per-chromosome tree sequences, in the `<chrom>/<rep>.<suffix>` layout produced by steps 5/5b), which it loads to measure retained sequence directly from each ARG.

## Auxiliary scripts

Scripts not called by the Snakemake pipeline. Run any with `--help` for its full
argument list, defaults, and examples.

- [`coalescence_ne_plots_from_ts.py`](../scripts/coalescence_ne_plots_from_ts.py) — pair-coalescence and Ne plots from TS replicates. Choose the time grid with one of `--time-bins-file` (explicit edges), `--num-quantiles N` (equal-coalescence-event bins derived from `pair_coalescence_quantiles` averaged over post-burnin replicates), or `--num-bins N` (uniform log-spaced bins across the observed coalescence-time range). Note that `--num-bins` meant equal-coalescence-event bins up to v1.8; that mode is now `--num-quantiles`, so pre-v1.9 commands need renaming to reproduce their old time grid. Optional Demes-based coalescent simulations (`--sim N`) produce window-stat and SFS TSVs for observed-vs-sim density plots in `validation_plots_from_ts.py`. Full option list and outputs: [coalescence_ne_plots.md](coalescence_ne_plots.md).
- [`compare_trees_html.py`](../scripts/compare_trees_html.py) — render one tree index from each of two tree sequences side-by-side into a single HTML file.
- [`locate_tree.py`](../scripts/locate_tree.py) — find the local tree at a `(chromosome, position)` in a merged tree sequence, using the `chrom_offsets` metadata to map the within-chromosome coordinate to the concatenated axis. See [Locating a tree by chromosome and position](outputs.md#locating-a-tree-by-chromosome-and-position).
- [`trees_gallery_html.py`](../scripts/trees_gallery_html.py) — scrollable HTML gallery of all trees from two tree sequences, useful for quick before/after inspection.
- [`simulate_two_bottleneck_demography.py`](../scripts/simulate_two_bottleneck_demography.py) — simulate replicate ARGs under a fixed two-bottleneck demography (35 ka + 9 ka bottlenecks, present-day expansion) for known-truth pipeline tests.
- [`make_realistic_example.py`](../scripts/make_realistic_example.py) — generate a realistic synthetic example dataset (ARGs from the two-bottleneck model + three injected flaws: contaminated individuals, per-window sample pruning, and a `mut_rate.p` accessibility mask) for end-to-end pipeline testing. Emits a `ground_truth.json` for scoring the pipeline's masks. See [realistic_example.md](realistic_example.md) for details, CLI options, and the ground-truth schema.

## Inputs, formats, defaults & logs

- **Tree-sequence files:** scripts accept `.ts`, `.trees`, and `.tsz` files. Loading/writing `.tsz` requires `tszip` to be installed; scripts will raise a clear error if `tszip` is missing when a `.tsz` is used.
- **BED files:** expected as whitespace-separated lines `chrom  start  end  [name]`. `start` and `end` are numeric (half-open intervals `[start, end)`). If a fourth column `name` is present it may list one or more comma-separated sample IDs; if omitted the BED filename stem is used as the sample name. Lines starting with `#` and blank lines are ignored.
- **HapMap recombination maps:** when required (e.g. `hapmap_low_rec_mask.py`), the script expects the HapMap format used by [`msprime.RateMap.read_hapmap`](https://tskit.dev/msprime/docs/stable/api.html#msprime.RateMap.read_hapmap).
- **Glob `--pattern`:** arguments named `--pattern` accept shell-style glob patterns (for example "*.tsz") and are matched against filenames in the supplied directory.
- **Defaults & output locations:** many scripts write to a `results/` directory or to an output directory under the input tree-directory when `--out`/`--out-dir` are not provided. Examples:
  - `trim_samples.py`: default output is `<ts_parent>/trimmed/<ts_stem>_trimmed.tsz` when `--out` is not given.
  - `mutload_summary.py` writes `results/<name>.html` and `logs/<name>.log`; no BED files are written (use `mutload_masks.py` for BED output).
  - Several plotting scripts write PNG files into `results/` by default; most have `--out` or `--out-dir` flags to override this.
- **Logging & errors:** Snakemake captures each rule's stdout/stderr under `logs/slurm/` when using the supplied SLURM profile. These scheduler/job logs are distinct from optional per-script summary logs written via a script's `--log` or `--out-dir` option. A script exception is therefore diagnosed first from its Snakemake job log; a summary log records completed-script results and is not guaranteed to exist after an early failure. Common failures include missing `tszip` for `.tsz` I/O and malformed BED records; BED parsing requires finite numeric coordinates, `start < end`, and rows within the sequence bounds.

## Shared module

`scripts/argtest_common.py` contains shared tree-sequence helpers used by multiple scripts:
- TS I/O (`load_ts`, `dump_ts`)
- mutational load / stat helpers
- trimming and masking helpers

Use this module for internal script imports rather than duplicating the helpers.

## Sample ID matching (`trim_samples.py`)

`trim_samples.py` matches sample/individual IDs exactly against the tree sequence's internal individual names as produced by `argtest_common.get_individual_name()` (prefers `individual.metadata['id']` when present, otherwise a synthetic `ind<id>` name).

- Matching is **exact** and **case-sensitive**.
- `--name-substring-to-remove` (default `""`) is removed via global string replacement before matching; provide names that match the normalized individual names. Despite replacing the former `--suffix-to-strip` option, the operation is not suffix-specific: every occurrence is removed.

Per-individual analyses may pool multiple sample nodes for one individual (for example, two nodes for a diploid). Uniform ploidy is assessed only among represented individuals; individual-table rows with no sample nodes are ignored. The stricter uniform-ploidy and leaf-sample contract remains audit/evidence-gated until it has been checked on the real input corpus.
