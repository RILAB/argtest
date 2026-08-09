# Changelog

All notable changes to this project are documented here. Versions correspond to
the annotated git tags (`git tag -l`). Dates are the tag dates.

## [Unreleased]

### Breaking

- Rename the sample-name normalization setting from `suffix_to_strip` to
  `name_substring_to_remove`, and the CLI option from `--suffix-to-strip` to
  `--name-substring-to-remove`. The behavior remains global `str.replace`
  removal rather than terminal-suffix removal. Existing configs must rename the
  key; the old key produces a migration error instead of silently using the
  default.

### Changed

- Document `trim_regions_single.py` as the coordinate-preserving pipeline
  step-4 filter; the obsolete coordinate-compacting directory CLI is removed.
- Mutation-load documentation now calls the reference a simulation-based
  expected load and states that it is estimated by one reproducible, seeded
  mutation simulation, so it retains simulation variance.
- Accessibility reporting consistently prefers `kept_intervals`, then the
  embedded positive-rate mutation map, then the documented whole-sequence
  fallback.
- Documentation distinguishes Snakemake-captured job stdout/stderr from
  per-script completed-run summaries and corrects `trim_samples.py`'s default
  output to `<ts_parent>/trimmed/<ts_stem>_trimmed.tsz`.
- The supported individual model allows multiple sample nodes per individual
  and assesses uniform ploidy only among represented individuals. Strict
  uniform-ploidy and leaf-sample enforcement remains contingent on a warn-only
  audit of the real ARG corpus.

## [v1.9] — 2026-07-27 — native missing-pair coalescence + final-ARG retention

### Breaking
- `coalescence_ne_plots_from_ts.py --num-bins N` has changed meaning. Through
  v1.8 it selected equal-coalescence-mass bins derived from
  `pair_coalescence_quantiles`; that behaviour is now `--num-quantiles N`, and
  `--num-bins N` selects uniform log-spaced bins across the observed
  coalescence-time range. Existing commands using `--num-bins` still run without
  error but produce a different time grid — rename them to `--num-quantiles` to
  keep the previous output.
- `pipeline_summary.py` gained a required `--filtered-ts` argument listing the
  final per-chromosome tree sequences (in the `<chrom>/<rep>.<suffix>` layout of
  step 5/5b), since retained sequence is now measured from the ARGs themselves.
  Standalone invocations must pass it; the Snakefile does so automatically.
  `step7_summary` consequently declares every filtered per-chromosome tree
  sequence as an input and loads each one, so it re-runs on existing pipeline
  output and its runtime and memory now scale with the ARGs rather than with the
  BED masks alone.
- `tskit` is no longer installed from conda-forge. `environment.yml` now
  installs it via pip from pinned nspope/tskit commit
  `73d8cd922482475020ae01180cae95bf5abbf067`, so building the environment
  requires pip, git, and access to GitHub. Upgrading an existing environment in
  place is unreliable — recreate it with `conda env remove -n argtest` followed
  by `conda env create -f environment.yml`.

### Added
- `coalescence_ne_plots_from_ts.py` now offers distinct automatic time-grid
  modes: `--num-quantiles N` for equal connected-pair coalescence mass and
  `--num-bins N` for uniform log-spaced bins across the observed time range.
- Coalescence/Ne runs export the plotted post-burn-in replicate trajectories
  and posterior means to `{prefix}coalescence-ne-estimates.tsv`, including
  pair-coalescence mass, log-density, rate, and effective population size.
- Pipeline summaries now include the exact argtest git version in their footer
  and show mask/retention values in Mb and as percentages of the corresponding
  full `.fai` chromosome length.

### Changed
- Pair-coalescence quantiles, counts, and rates now use the native
  partial-missing-data normalization from pinned nspope/tskit commit
  `73d8cd922482475020ae01180cae95bf5abbf067`. The script-level conditional
  quantile adjustment and `connected_pair_span` scalar rescaling have been
  removed, so locally isolated samples are handled correctly by tskit for
  global, spatial-window, and multiple-sample-set calculations.
- The runtime environment now requires Python >= 3.11 and NumPy >= 2, as
  required by the pinned tskit development build.
- Coalescence plots include only post-burn-in replicate trajectories in the
  displayed posterior ensemble, and input ARG replicates are ordered
  naturally by numeric index.
- The example Snakemake configuration now documents the optional `min_samples`
  step-5b filter.

### Fixed
- Final retained sequence in `pipeline_summary.html` is measured directly from
  each filtered ARG as positive-rate accessible sequence with a non-empty local
  genealogy. It therefore excludes both regions masked in the input ARG and
  regions removed by the pipeline.
- Genome-wide Step 3, pipeline-union, and retained Mb are now summed across
  chromosomes within each replicate before reporting mean +/- SD. Step 3 and
  pipeline-union columns also report percentages of chromosome length.
- `trim_regions_single.py` preserves mutation-rate maps already embedded in
  input tree-sequence metadata, ensuring downstream retention calculations can
  recover the original accessibility mask.

## [v1.8] — 2026-07-13 — chromosome-position tree lookup + pipeline hotspot optimizations

### Added
- Merged tree sequences now record a `chrom_offsets` table
  (`[{chrom, offset, length}, ...]`) in their top-level metadata, capturing
  where each chromosome lives on the concatenated coordinate axis. New
  `argtest_common` helpers map between within-chromosome and merged coordinates
  without re-deriving offsets from the per-chromosome inputs:
  `genome_position`, `chrom_position_from_genome`, `tree_at_chrom_position`, and
  `chrom_offsets_from_metadata`.
- `scripts/locate_tree.py`: CLI to find the local tree at a
  `(chromosome, position)` in a merged tree sequence, reporting the genome
  coordinate, covering tree index/interval, and whether the position falls in a
  masked/trimmed region (`num_edges == 0`). Documented in the README. Merged
  files produced before this release lack `chrom_offsets` and raise a clear
  `KeyError` pointing to a re-merge or the per-chromosome fallback.

### Performance
- `trim_samples` trimmed-sample mutation filtering is vectorized: the per-site,
  per-mutation Python loop that rebuilt the site and mutation tables is replaced
  by a numpy `searchsorted` grouping over drop intervals plus table
  `keep_rows`, dropping the same mutations without materializing every kept row
  in Python.
- `merge_treefiles_by_replicate` narrows file discovery to the requested
  replicate up front (`find_tree_files_for_replicate`) instead of scanning and
  grouping the whole tree directory, and loads/concatenates chromosomes
  incrementally rather than materializing every per-chromosome tree sequence
  before merging.
- `hapmap_low_rec_mask` resolves and loads only the requested chromosome's rows
  (scanning the `Chromosome` column for the key, then filtering while loading)
  instead of parsing the entire hapmap and subsetting afterward.
- `filter_min_samples` computes each dropped interval's minimum retained-sample
  count with a single ordered sweep over the per-tree counts, replacing the
  per-interval rescans; the `kept_intervals` intersection uses an advancing
  two-pointer sweep instead of re-scanning all dropped spans per kept interval.
- `coalescence_ne_plots_from_ts` accumulates the masked PDF/rate arrays per
  replicate as it reads them, rather than stacking the full-length arrays for
  all replicates and slicing afterward, lowering peak memory. Companion
  streamlining in `validation_plots_from_ts`.

### Fixed
- `genome_position` now rejects non-finite positions (`NaN`, `inf`) with the
  same out-of-range `ValueError` instead of letting them slip past the bounds
  check.

## [v1.7] — 2026-07-04 — scalar mutation-rate fallback + trimmed-sample mutation dropping

### Added
- Scalar mutation-rate fallback for `find_low_access_regions.py`
  (`--mutation-rate`, wired through the Snakefile `step2_low_access` rule). The
  accessibility rate is now resolved via `argtest_common.resolve_mu_rate` in
  priority order: tree-sequence metadata ratemap → sibling `*.mut_rate.p` file →
  the scalar fallback. A positive scalar marks the whole sequence accessible; a
  non-positive scalar marks none of it accessible (every window flagged
  low-access). Replaces the previous unconditional `*.mut_rate.p` pickle load.

### Changed
- `trim_samples` now drops mutations carried by trimmed sample nodes within
  their trim intervals, so a removed individual no longer retains phantom
  variants over regions it was trimmed from. Dropping is restricted to leaf
  sample nodes; an internal sample node (one that appears as an edge parent)
  keeps its mutations because untrimmed descendants still inherit them.
- SLURM profile: a separate concurrency cap for the RAM-heavy merge step so it
  no longer contends with lighter per-chromosome jobs.
- Aligned the SLURM controller default configfile and docs with the
  consolidated `config/snakemake.yaml` template.

### Fixed
- `pipeline_summary`: the genome-wide retained-sequence percentage is now
  length-weighted per replicate (summing retained and total bp across
  chromosomes before dividing) instead of an unweighted mean of per-chromosome
  percentages, which skewed genomes with unevenly sized chromosomes. Outlier
  summaries now count replicates in which an individual had zero outliers,
  instead of averaging only over the replicates where it was flagged.

## [v1.6] — 2026-07-01 — minimum-retained-samples filter + configurable validation window

### Added
- Optional minimum-retained-samples filter (`min_samples`): a new step 5b
  (`scripts/filter_min_samples.py`), run between `trim_samples` and the merge,
  that drops local-tree intervals retaining fewer than `min_samples`
  non-isolated sample nodes. Dropping uses `delete_intervals`, so sequence
  coordinates are preserved (removed spans become empty gaps, no compaction),
  and the `kept_intervals` metadata is intersected with the surviving spans so
  downstream accessibility is not overestimated. Each (chromosome, replicate)
  also gets a diagnostic BED of the dropped intervals. When set, the merge, VCF
  export, and validation steps consume the filtered tree sequences. Unset/null
  by default, so existing configs are unaffected. Retained-sample counts are
  computed with a vectorized edge-endpoint coverage sweep rather than a
  per-tree, per-sample loop.
- `validation_window_size` config key: window size in bp for the step-6
  diversity, Tajima's D, and segregating-sites validation plots (default:
  `100000`). Larger windows run faster and use less memory on large ARGs at the
  cost of coarser QC curves. The standalone `validation_plots_from_ts.py`
  default was aligned to 100000 to match (previously 50000).

### Fixed
- `scripts/coalescence_ne_plots_from_ts.py` now raises an informative error when
  the empirical pair-coalescence CDF has too few independent deep coalescence
  events to define the requested equal-mass `--num-bins` (naming the cause and
  directing the user to `--time-bins-file`), instead of the opaque "quantile
  edges are not strictly increasing".

## [v1.5] — 2026-06-26 — VCF export + shared mutation-rate config

### Added
- Optional VCF export (`emit_vcf`, `vcf_reps`): one `.vcf.gz` per
  (chromosome, replicate) from the trimmed tree sequences via
  `scripts/export_vcf.py` — variable sites only, ploidy-aware genotypes
  (haploid coded as a single allele `0`/`1`, diploid as `0|1`), and samples
  pruned by `trim_samples` written as missing (`.`) via `isolated_as_missing`.
- Shared `mutation_rate` config key used as the default for both the step-3
  mutload sim and step-6 validation, with per-step overrides
  (`mutload_mutation_rate`, `validation_mutation_rate`).

### Changed
- Consolidated to a single documented `config/snakemake.yaml` template;
  removed older dataset-specific config templates.
- Removed the superseded `sim_2chr` example data (the realistic example is the
  default).
- Documented file naming / chromosome-label matching (hapmap `Chromosome` and
  `.fai` lookup, `*.mut_rate.p` discovery), and clarified that `rec_fraction`
  masks a fraction of recombination-rate intervals, not base pairs.

### Fixed
- Validation plots no longer crash on an empty/all-zero site spectrum (the
  y-axis falls back to linear with a warning instead of failing at savefig).
- Corrected the false "`validation_sim_branch` requires a `*.mut_rate.p` file"
  note (it falls back to the scalar rate) and the stale "5 replicates" burn-in
  comment (the bundled example has 8).

## [v1.4] — 2026-06-03 — per-sample expected load + trim_samples scaling

### Added
- Validation report now plots a per-sample simulated expected mutational load
  (with a 95% interval) under `--sim-branch`, replacing the single flat
  genome-wide expected line so each sample is compared to its own
  genealogy-based expectation.
- SLURM status command (`profiles/slurm/status-sacct.sh`) wired into the
  cluster-generic profile so the controller detects jobs SLURM kills before any
  output is written (TIMEOUT/OUT_OF_MEMORY/NODE_FAIL) instead of polling forever.

### Changed
- Rewrite step5 `trim_samples` as a single vectorized simplify pass instead of
  one `simplify()` per removal interval per individual, making it tractable on
  large ARGs (the admix ARG previously never finished a single step5 job).
  Output verified equivalent (genealogy + genotypes) by a 12-check regression
  test.
- README: document that mutation-load outlier detection is per individual
  (sample nodes pooled), which equals per sample for haploid data.

### Fixed
- Empty `base_name` now falls back to the `root_dir` name instead of producing
  merge output names that mismatched the Snakefile targets and aborted
  `merge_replicates` with a `MissingOutputException`.

## [v1.3] — 2026-05-27 — SLURM cluster execution + realistic example dataset

### Added
- SLURM execution: all cluster knobs in a single config, validated resource
  config, per-rule time defaults, shared-FS-friendly `out_dir`.
- `tree_subdir` support for nested per-chromosome tree-file layouts.
- `measure_merge_mem.py` to size the memory-heavy merge step for new datasets.
- Realistic example dataset: `make_realistic_example.py` (with per-individual
  contamination hotspots, `all.hapmap.tsv` and `sim.fai` emission),
  `score_realistic_example.py` for pipeline-vs-ground-truth scoring,
  `MAKE_REALISTIC_EXAMPLE.md`, and tracked lightweight artifacts.

### Changed
- Match chr-prefixed hapmap chromosomes to pipeline chrom names.
- Resolve staged ts symlinks to absolute paths.
- Quadruple per-rule time defaults; raise default-resources fallback to 04:00:00.
- Clarify that the input dataset is configurable (not always amaranth).
- README: document `tree_subdir`, clarify run instructions, drop stale QC TODO,
  add workflow details.

### Fixed
- Ne bias on tree sequences with isolated samples.
- SFS double-normalization in validation plots.
- Pair-coalescence rescale caveat documented.

## [v1.2] — 2026-05-14 — per-individual mutload flagging + auto Ne bins

### Added
- Per-window per-individual mutation-load flagging.
- Automatic Ne time bins (`--num-bins` for Ne plots).
- Two-bottleneck demography simulation script.
- Segregating-sites trace in validation plots.

### Changed
- Switch mutload to a simulation-based expectation.
- Tighten mutload outlier semantics and surface mu-map ambiguity.

## [v1.1] — 2026-04-25 — validation-plots cleanup + robustness pass

### Changed
- Drop unsimulated branch-mode stats from validation plots.
- Offset-merge `kept_intervals` across chromosomes.
- Fix treefile suffix set, allow an empty low-rec mask, document stat caveat.
- README: citation, acknowledgements, and current-pipeline updates.

### Added
- Warning that tskit coalescence rates are biased on partial trees.
- Simulation-based mutload plan.

## [v1.0] — 2026-04-22 — initial Snakemake pipeline release

### Added
- Snakemake pipeline (steps 1–7), integration tests, and a simulated example
  dataset (`sim_2chr_5rep`).
- `--sim-branch` mode: posterior predictive check via msprime simulation.
- `step7_summary`: self-contained pipeline summary HTML.
- Always produce folded + unfolded SFS; embed ratemap in ts metadata.
- Burn-in, step6 validation rule, accessible-bp diversity scaling.
- Centralized `infer_mu_path`.

### Fixed
- `OUT_DIR` tilde expansion in config paths.
- Hapmap loading: skip repeated headers, resolve cross-naming-convention chroms.
- `TSK_ERR_BAD_MUTATION_PARENT` in `remove_ancestry`.
- Normalize mutload by tree-covered accessible bp; metadata schema bugs.

[v1.9]: https://github.com/RILAB/argtest/releases/tag/v1.9
[v1.8]: https://github.com/RILAB/argtest/releases/tag/v1.8
[v1.7]: https://github.com/RILAB/argtest/releases/tag/v1.7
[v1.6]: https://github.com/RILAB/argtest/releases/tag/v1.6
[v1.5]: https://github.com/RILAB/argtest/releases/tag/v1.5
[v1.4]: https://github.com/RILAB/argtest/releases/tag/v1.4
[v1.3]: https://github.com/RILAB/argtest/releases/tag/v1.3
[v1.2]: https://github.com/RILAB/argtest/releases/tag/v1.2
[v1.1]: https://github.com/RILAB/argtest/releases/tag/v1.1
[v1.0]: https://github.com/RILAB/argtest/releases/tag/v1.0
