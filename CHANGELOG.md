# Changelog

All notable changes to this project are documented here. Versions correspond to
the annotated git tags (`git tag -l`). Dates are the tag dates.

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

[v1.3]: https://github.com/RILAB/argtest/releases/tag/v1.3
[v1.2]: https://github.com/RILAB/argtest/releases/tag/v1.2
[v1.1]: https://github.com/RILAB/argtest/releases/tag/v1.1
[v1.0]: https://github.com/RILAB/argtest/releases/tag/v1.0
