# Changelog

All notable changes to this project are documented here. Versions correspond to
the annotated git tags (`git tag -l`). Dates are the tag dates.

## [v1.12] — 2026-08-19 — documented coalescence plotting, no `--time-adjust`

### Breaking

- Remove `--time-adjust` from `coalescence_ne_plots_from_ts.py`. Plot time axes
  are always in generations, so the axis labels drop "Adjusted" and the
  `adjusted_time_left` / `adjusted_time_right` columns and the `time_adjust`
  line are gone from `coalescence-ne-estimates.tsv` and `summary.txt`.
  Convert to calendar time downstream if needed.

### Added

- Boundary tests for `mutload_masks.py --fraction` covering an outlier fraction
  below, exactly equal to, and above the threshold. The comparison is strict
  (`fraction > threshold`), so a window whose outlier fraction exactly equals
  `--fraction` is deliberately *not* masked; the equality test fails if that
  `>` is ever changed to `>=`.

### Fixed

- `coalescence_ne_plots_from_ts.py --help` (and argument-error reporting) now
  works without the pinned tskit fork installed: arguments are parsed before
  the fork probe runs, instead of the probe aborting the run first.

### Changed

- `pipeline_summary.py` now streams BED and FAI files line by line through a
  shared `iter_data_lines` helper instead of `read_text().splitlines()`, so
  whole-genome masks are no longer held in memory in full. Blank and `#`
  comment lines are skipped for FAI input too, which previously would have
  raised on a commented FAI.
- Correct two `coalescence_ne_plots_from_ts.py` help strings that hard-coded
  the simulation defaults (1 Mb sequence, 50 Kb windows) as if they were fixed;
  they point at `--sim-length` and `--sim-window-size`.
- Document every `coalescence_ne_plots_from_ts.py` option in the README, in a
  new [Coalescence and Ne plots](README.md#coalescence-and-ne-plots) section
  covering the time-grid modes, estimation and plotting flags, the Demes
  simulation flags, and every output file with the estimates-TSV columns.

## [v1.11] — 2026-08-14 — fail-fast chromosome-name validation

### Added

- Validate every ARG chromosome-directory label against both the HapMap
  `Chromosome` column and the FASTA-index chromosome names while evaluating the
  Snakefile. A mismatch now reports all unmatched and available names before
  the DAG is built or any local/SLURM job is submitted, instead of failing one
  chromosome at a time in step 1.
- Document the supported chromosome aliases and validation scope. The workflow
  has no input VCF; exported VCF contig names are derived from the validated ARG
  chromosome label and therefore agree by construction.

## [v1.10] — 2026-08-09 — simplification overhaul: explicit inputs, fewer moving parts

A broad simplification pass driven by an independent code review: one supported
implementation per task, and fail-fast errors in place of silent fallbacks.

Measured against v1.9 across `scripts/`, the `Snakefile`, and config — excluding
tests and documentation — this is **205 fewer total lines (−2.6%)** and **263
fewer logical lines (−4.4%)**, ignoring blanks and comments. Roughly 500 logical
lines were deleted (two whole modules, eleven dead helpers in `argtest_common`,
and consolidation across seven scripts) and roughly 240 added back, mostly as
new diagnostics: the mutation-map checker, the ARG contract audit, mutation-rate
source reporting, and the pinned-dependency probe. The structural change is
larger than the line delta suggests — no shared helper retains a test-only
caller, and there is now a single supported implementation of region trimming
and of sample trimming.

This is a minor bump carrying breaking changes, following the same convention as
v1.9. **Read the Breaking section before upgrading**: several changes turn
previously-working configurations into hard failures, and one statistic — the
retained-sequence figure in `pipeline_summary.py` — changes value, as a
correction. Existing configs must at minimum rename `suffix_to_strip`.

### Breaking

- Rename the sample-name normalization setting from `suffix_to_strip` to
  `name_substring_to_remove`, and the CLI option from `--suffix-to-strip` to
  `--name-substring-to-remove`. The behavior remains global `str.replace`
  removal rather than terminal-suffix removal. Existing configs must rename the
  key; the old key produces a migration error instead of silently using the
  default.
- **Step 4 now requires a resolvable mutation rate.** `trim_regions_single.py`
  always embeds a rate map, so it fails when none of an embedded map, an exact
  sibling `*.mut_rate.p`, or the new `--mutation-rate` scalar is available.
  Previously the map was optional and the step silently produced output without
  rate metadata. The Snakefile forwards the config `mutation_rate` when set; a
  run with no mutation-rate source anywhere must now set one.
- **Mutation-map discovery accepts only exact paths.** `infer_mu_path` tries
  `<ts_stem>.mut_rate.p` and `<parent_dir_name>.mut_rate.p`, beside the treefile
  and one directory above it. The trailing-replicate/numeric-suffix base
  derivation and the directory-wide `*.mut_rate.p` glob are gone, as is the
  resulting "Ambiguous mutation maps" error — a missing map now reports every
  path tried. Datasets that relied on a broad shared base name (for example a
  single `amaranth.mut_rate.p` serving `amaranth.1`…`amaranth.16`) must supply
  per-chromosome or per-replicate files, or set a scalar `mutation_rate`.
- `resolve_mu_rate` no longer rebuilds a `RateMap` from any pickled object that
  happens to expose `.position` and `.rate`. Pickles must contain an
  `msprime.RateMap`.
- `find_low_access_regions.py` gained a **required** `--chrom`. It previously
  wrote `ts_path.stem` into BED column 1, which is the *replicate* ID under the
  nested `<chrom>/<rep>.<suffix>` layout the pipeline uses.
- `combine_remove_masks.py` now fails on bad input instead of skipping it: a
  missing input BED raises `FileNotFoundError`, and a non-blank, non-comment row
  with fewer than three fields or non-numeric coordinates raises `ValueError`
  with file and line context. Empty BED files remain valid.
- **`pipeline_summary.py` retained-sequence accounting changed.** The accessible
  set is now `kept_intervals`, falling back to the embedded positive-rate
  mutation map and then the whole sequence; previously it always used the
  mutation map. Both are still intersected with tree-covered spans, so trimmed
  regions are excluded either way — the difference is `mu == 0` sequence that
  survived step 4 and still carries genealogy, which v1.9 excluded and v1.10
  counts. **Reported retained bp increases** for the same ARGs, by an amount
  bounded by the zero-rate fraction of the mutation map (33% of the sequence in
  the bundled realistic example, though the realized delta is smaller since
  step 2 masks most low-accessibility windows). This is a correction, not a
  redefinition: `kept_intervals` is what the pipeline actually retained, and it
  is what `validation_plots_from_ts.py` already used — the two now agree.
- `coalescence_ne_plots_from_ts.py` refuses to start unless the pinned
  nspope/tskit fork is active, rather than silently producing wrong
  partial-missing-data normalization on stock tskit.
- Removed `scripts/trim_regions.py` (the coordinate-*compacting* directory CLI)
  and `scripts/trim_samples_chunked.py`. `trim_regions_single.py` is the single
  supported, coordinate-preserving step-4 filter; its `load_mask_intervals` and
  `complement_intervals` helpers moved to `argtest_common`.
- `argtest_common` lost eleven helpers with no live callers:
  `validate_trimmed_ts`, `assert_sample_ids_preserved`, `trim_ts_by_intervals`,
  `build_removal_segments`, `build_segments_with_drop_nodes`, `infer_mu_base`,
  `build_shared_mask`, `and_ratemaps_binary`, `collapse_masked_intervals`,
  `collapse_masked_and_low_access_windows`, and `ratemap_from_keep_intervals`.
  `mutational_load` also dropped its unused `remove_intervals` /
  `name_to_nodes` parameters.

### Added

- **The mutation-rate source is now recorded and reported.** Resolution order is
  unchanged (embedded ratemap → exact sibling `*.mut_rate.p` → scalar
  `mutation_rate`), but step 4 stamps which one it used into the output ARG's
  `mu_source` metadata, prints it, and warns on stderr when it fell back to the
  scalar. `pipeline_summary.py` gained a **Mutation-rate source** section that
  counts sources across the run and flags every ARG that used the flat fallback.
  Without this, a rate map that fails discovery is silently replaced by a flat
  rate — which removes the local-rate correction the step-3 outlier test exists
  to apply, changing outlier calls with no error anywhere.
- `scripts/check_mu_paths.py` dry-runs mutation-map discovery over a dataset
  using only `stat` calls (no ARG data loaded, safe on a head node) and exits
  non-zero if any tree file would fail step 4. Use it before a run; use the
  summary section above to audit a run after the fact.
- `scripts/audit_arg_contract.py` and `argtest_common.audit_individual_contract`
  report — **as warnings only** — sample nodes with no individual, duplicate
  normalized individual names, mixed ploidy among represented individuals, and
  sample nodes used as edge parents. Pipeline behavior does not depend on the
  result. This is the evidence-gathering step before any of those conditions
  becomes an error.
- `validation_plots_from_ts.py` accepts explicit `--ts FILE [FILE ...]` and
  `--compare-ts FILE [FILE ...]` file lists, mutually exclusive with the
  directory-plus-glob forms.
- `trim_regions_single.py` gained `--mutation-rate` as a positive scalar
  fallback.
- Tests for the individual/ploidy audit (`tests/test_individual_contract.py`)
  and for validation-plot argument handling (`tests/test_validation_plots.py`).

### Fixed

- Step 2 wrote the replicate ID rather than the chromosome into BED column 1
  under the nested input layout (see the required `--chrom` above).
- Step 4 replaced a tree sequence's top-level metadata wholesale; unrelated
  provenance fields are now preserved and only `kept_intervals` and the rate map
  are overwritten.
- `find_low_access_regions.py` and `trim_samples.py` did not create the output
  parent directory when given an explicit `--out`, so standalone runs into a
  fresh tree failed. (Pipeline runs were unaffected — Snakemake pre-creates
  output parents.)
- `validate_trimmed_ts` never actually validated anything: its `check_index`
  fallback could not obtain an index, returned early, and emitted a misleading
  "skipped check_index" warning on every normal step-5 and step-5b job. It has
  been deleted rather than repaired; malformed tables still fail when
  `tables.tree_sequence()` constructs the result.

### Performance

- `merge_replicates` concatenates all chromosomes of a replicate in a single
  `concatenate(*remaining)` call instead of rebuilding a progressively larger
  intermediate at every chromosome boundary, and releases the merged tree
  sequence before reconstructing metadata. Peak RSS on real data has not yet
  been measured on the pinned environment.
- `mutational_load` iterates the input tree sequence directly. It previously
  called `keep_intervals([(0, sequence_length)], simplify=False)` on the
  supported path, copying and reindexing the entire tree sequence without
  changing it.
- `hapmap_low_rec_mask.py` parses the HapMap once per invocation instead of
  twice, halving I/O for the per-chromosome rule.
- Step 6 invokes the validation script with explicit file arguments, removing
  the `/tmp` staging directories, symlink loops, and inline Python map lookup
  from the rule.
- Plotting scripts close figures with `plt.close(fig)` rather than `plt.clf()`,
  so pyplot no longer retains every cleared figure for the life of the process.

### Changed

- BED interval bounds are rounded outward: `floor(start)` / `ceil(end)` instead
  of truncation with `int()`, so a mask never loses covered sequence to inward
  truncation. In practice this is close to a no-op — `int()` already equals
  `floor()` for non-negative coordinates, so **start positions never change**,
  and ends move by at most 1 bp only where a window edge is non-integral. That
  requires `--snp-window` mode (whose edges are site positions) *and* an ARG with
  non-integral site positions, i.e. a continuous-genome msprime simulation.
  Discrete-genome and SINGER-derived ARGs produce byte-identical output.
- Duplicated helpers are consolidated in `argtest_common`: base-pair window
  construction, SNP window construction, and the seeded simulation-based
  expected-load calculation each have one implementation. The HTML tree-view
  scripts and `score_realistic_example.py` now import the shared `load_ts` and
  `merge_intervals` instead of carrying private copies.
- Mutation-load documentation calls the reference a **simulation-based expected
  load** and states that it is estimated from one reproducible, seeded mutation
  simulation, so it carries simulation variance and is not an analytic
  expectation.
- Documentation distinguishes Snakemake-captured job stdout/stderr from
  per-script completed-run summaries and corrects `trim_samples.py`'s default
  output to `<ts_parent>/trimmed/<ts_stem>_trimmed.tsz`.
- README gained a "Supported input assumptions" section. The supported
  individual model allows multiple sample nodes per individual and assesses
  uniform ploidy only among **represented** individuals (individual-table rows
  with no sample nodes are ignored). Strict uniform-ploidy and leaf-sample
  enforcement, and removal of the VCF one-column-per-node fallback, remain
  contingent on the warn-only audit passing against the real amaranth/admix
  corpus.

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

[v1.11]: https://github.com/RILAB/argtest/releases/tag/v1.11
[v1.10]: https://github.com/RILAB/argtest/releases/tag/v1.10
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
