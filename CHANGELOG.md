# Changelog

All notable changes to this project are documented here. Versions correspond to
the annotated git tags (`git tag -l`). Dates are the tag dates.

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

[v1.4]: https://github.com/RILAB/argtest/releases/tag/v1.4
[v1.3]: https://github.com/RILAB/argtest/releases/tag/v1.3
[v1.2]: https://github.com/RILAB/argtest/releases/tag/v1.2
[v1.1]: https://github.com/RILAB/argtest/releases/tag/v1.1
[v1.0]: https://github.com/RILAB/argtest/releases/tag/v1.0
