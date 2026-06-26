# Upstream vs Fork: per-window / diagnostic statistics comparison

Verification status (2026-05-14): file:line citations and §1–§4 claims
spot-checked against current fork tree. §5 empirical questions resolved as
follows — see updated §5 for details:

- **Q#1 (SFS double-normalization): CONFIRMED BUG.** tskit 1.0.2 defaults
  `span_normalise=True` for site-mode `allele_frequency_spectrum` (verified
  with a 1 Mb msprime sim). Fork divides by sequence length twice — the
  implicit default once, then by `total_accessible` explicitly. Tracked as a
  new todo in `project_todos.md`.
- **Q#2 (kept_intervals & tree-edge gaps): NO.** `trim_regions_single.py:38`
  writes `kept_intervals` as the pure inverse of the mask BED with
  `simplify=False`; no edge-coverage intersection. Per-window pi may be
  biased high inside windows with internal edge gaps; SFS divisor mitigates
  via `tree_covered_accessible_bp` (separate from the double-norm bug).
- **Q#3 (partial-tree coalescence bias): STRUCTURALLY CONFIRMED.**
  `coalescence_ne_plots_from_ts.py` does not call `collapse_masked_intervals`,
  and `trim_samples.py:remove_ancestry` produces partial trees inside
  `kept_intervals`. Magnitude on a real trimmed run still empirical —
  remains under the existing partial-trees Ne fix todo.
- **Q#4 (mutload seed plumbing): OK.** `Snakefile:83 mutload_seed_for`
  derives per-(chrom, rep) seeds via sha1; script's default of 1 only
  affects manual invocations.

Other notable verifications:
- `collapse_masked_intervals` (argtest_common.py:566) is used by the legacy
  multi-mask `trim_regions.py:125`, NOT by the active `trim_regions_single.py`
  step. After `trim_samples_single`, local trees inside kept_intervals are
  partial (some samples excised over per-individual intervals).
- All AFS call sites in the fork are listed in the empirical Q#1 block;
  the asymmetric overlay between observed (double-normalized) and
  `sim-sfs.tsv` (single-normalized) is the visible-on-log-scale risk.

Scope: Comparison between this fork (`/home/jri/src/argtest`) and upstream
[`nspope/singer-snakemake`](https://github.com/nspope/singer-snakemake) at the
level of source code only. No empirical runs were performed. Citations point
to specific lines in both code bases (upstream files were fetched 2026-05-09
into `/tmp/upstream/`).

---

## 1. Stats produced

| Stat | Upstream | Fork | Notes |
|---|---|---|---|
| Per-window nucleotide diversity (pi) — observed | yes (VCF, scikit-allel) | yes (site-mode on TS) | Different data source: VCF vs TS sites. Fork has no VCF stats path. |
| Per-window nucleotide diversity — expected on ARG | yes (sim-on-trees, site-mode) | yes (`--sim-branch`, site-mode on simulated TS) | Same idea, different framing — see below. |
| Per-window Tajima's D — observed | yes (VCF, scikit-allel) | yes (site-mode on TS) | Same caveat. |
| Per-window Tajima's D — expected on ARG | yes (sim-on-trees) | yes with `--sim-branch` | Same idea. |
| Genome-wide SFS — observed | yes (VCF, scikit-allel folded/unfolded) | yes (site-mode on TS) | Fork emits both folded/unfolded in one run; upstream picks per `polarised` config. |
| Genome-wide SFS — expected on ARG | yes (sim-on-trees) | yes with `--sim-branch` | Same idea. |
| Mutational load per sample — observed | yes (`utils.mutational_load` on raw trees) | yes (`argtest_common.mutational_load`) | Fork's helper is richer (windows, intervals, drop sets). |
| Mutational load per sample — expected on ARG | yes (`utils.mutational_load` on sim-on-trees) | partial: validation plots compute a single `load_vals` per replicate (no separate "expected" trace); mutload-mask scripts compute sim expected | Fork drops the upstream observed-vs-expected scatter. |
| Mutational load posterior trace | yes | yes | Fork additionally adds a horizontal mean line and supports `--compare`. |
| Repolarised-fraction trace | yes (`stats['repolarised']`) | no | Fork's TS pipeline doesn't have a `flipped` site-metadata field, so this isn't reproducible without changes. |
| Multimapped-sites mean | yes (`stats['multimapped']`) | no | Same reason. |
| Stratified divergence (between subpops) | yes | no | Fork has no `--stratify` plumbing in the validation/coalescence scripts. |
| Stratified SFS (per-subpop) | yes | no | Same. |
| Cross-population coalescence rates / cross-PDF | yes (`coalescence_rates.py` + `plot_coalescence_rates.py`) | no | Fork's `coalescence_ne_plots_from_ts.py` is single-population only. |
| Pair coalescence rates (single pop) | yes | yes | Different time-grid construction, different masking. |
| Pair coalescence PDF | yes | yes | Fork converts PMF to log-density. Upstream plots the raw bin PDF. |
| Effective population size Ne(t) = 1/(2 rate) | NO direct plot | yes (`--out-dir/effective-pop-size.png`) | Fork-only. |
| Demes-graph from inferred Ne and posterior-predictive simulations | NO | yes (`--sim` flag) | Fork-only. Builds a Demes model from mean Ne, simulates 1 Mb chunks with msprime, dumps `sim-window-stats.tsv` and `sim-sfs.tsv`. |
| Mutational-load outlier BED masks per chromosome | NO | yes (`scripts/mutload_masks.py`) | Fork-only filtering step. |
| Mutational-load HTML summary (per-individual + lineage breakdown) | NO | yes (`scripts/mutload_summary.py`) | Fork-only. |
| Per-replicate "kept-intervals" accessibility metadata | NO | yes (read by validation_plots; populated by `trim_regions_single`) | Fork-only. |
| Genealogical-gap detection from trees | yes (`utils.find_genealogical_gaps`) | partial (tree-coverage via `tree_covered_accessible_bp`) | Fork uses tree coverage but does not re-emit a gap interval list. |

---

## 2. How shared stats are computed (upstream vs fork)

### 2.1 Per-window nucleotide diversity (pi)

- **Upstream observed** (`workflow/scripts/chunk_chromosomes.py:454-460`):
  computed once from the VCF using `allel.windowed_diversity` with explicit
  `is_accessible=statmask`, `fill=0.0`. Windows are `int` bp ranges spanning
  the whole bitmask; sites omitted from dating are excluded via
  `statkeep = statmask[positions - 1]`. Result is masked to NaN at
  `~filter_windows` (`chunk_chromosomes.py:645`). Stored in
  `vcf_stats["diversity"]`.
- **Upstream "expected" / branch-equivalent**
  (`workflow/scripts/tree_statistics.py:68-74`): runs
  `msprime.sim_mutations(trees, rate=adjusted_mu, keep=False)`, then calls
  `ts.diversity(mode='site', windows=windows.position, span_normalise=False)`.
  Per-base normalisation done explicitly via
  `accessible_bp = np.diff(accessible.get_cumulative_mass(windows.position))`,
  and only on windows where `windows.rate == 1.0`; others set to NaN.
  `accessible` is the AND of the inaccessibility ratemap and
  `extract_accessible_ratemap(trees)` (`utils.py:207-217`), which in turn
  zeros every tree interval that has no edges.
- **Fork observed**
  (`scripts/validation_plots_from_ts.py:342, 348`): `ts.diversity(mode="site",
  windows=rep_windows)` directly on the SINGER-output tree sequence, then
  rescales with `div_scale = window_spans / rep_acc` (`:331-333`) where
  `rep_acc` is the per-window accessible bp derived from either
  `kept_intervals` metadata (preferred) or the mu ratemap. Windows with
  zero accessible bp are NaN'd. The "observed" pi here therefore counts
  whatever sites SINGER chose to put on the topology, which under default
  config behaviour mirrors the input VCF.
- **Fork "expected" / sim-branch**
  (`scripts/validation_plots_from_ts.py:362-376`): same pattern as upstream
  — `msprime.sim_mutations(ts, rate=_sim_rate, keep=False, random_seed=...)`
  then `sim_ts.diversity(mode="site", ..., span_normalise=False)` divided by
  `rep_acc`.
- **Diff**: identical underlying algorithm for the sim path. The main
  divergences are (a) upstream uses a fixed pre-computed `windows.rate ∈
  {0,1}` mask while fork rebuilds accessible bp per replicate, (b) fork's
  fallback chain (kept_intervals → mu ratemap → full sequence) is more
  flexible but means two replicates of the same chromosome may have
  slightly different denominators if their kept_intervals differ. **Open
  question**: does upstream's `accessible` differ from fork's intersection
  in practice? See §4.

### 2.2 Per-window Tajima's D

- **Upstream observed**: `allel.windowed_tajima_d` from VCF
  (`chunk_chromosomes.py:461-465`, no `is_accessible` argument).
- **Upstream sim**: `ts.Tajimas_D(mode='site', windows=windows.position)` on
  the simulated TS (`tree_statistics.py:76-80`); NaN'd where
  `windows.rate == 0.0`.
- **Fork observed**: `ts.Tajimas_D(mode="site", windows=rep_windows)`
  (`validation_plots_from_ts.py:343`); NaN'd where `np.isnan(div_scale)`.
- **Fork sim**: identical algorithm to upstream
  (`validation_plots_from_ts.py:368-372`).
- **Diff**: NaN-handling rules differ slightly (fork uses div_scale from
  accessible bp; upstream uses the ratemap window mask). For Tajima's D
  itself there is no per-base normalisation issue — only a per-window
  accept/reject. Upstream's VCF Tajima's D and fork's site-on-TS Tajima's D
  should agree exactly when the SINGER-output TS contains exactly the input
  VCF sites and there is no missing data.

### 2.3 Site frequency spectrum (SFS)

- **Upstream observed**: `allel.sfs` (polarised) or `allel.sfs_folded` from
  VCF, normalised by `np.sum(statmask)` (`chunk_chromosomes.py:466-469`).
  One mode per run, controlled by `polarised` config flag.
- **Upstream sim**: `ts.allele_frequency_spectrum(mode='site',
  span_normalise=False, polarised=...) / accessible_bp[windows.rate==1.0].sum()`
  (`tree_statistics.py:82-86`).
- **Fork observed**: emits BOTH polarised and folded in one run
  (`validation_plots_from_ts.py:347-356`):
  `ts.allele_frequency_spectrum(mode="site", polarised=True/False) /
  total_accessible`, where `total_accessible =
  tree_covered_accessible_bp(ts, acc_intervals)`
  (`argtest_common.py:428-451`).
- **Fork sim**: same divisor, `sim_ts.allele_frequency_spectrum(...) /
  total_accessible` (`validation_plots_from_ts.py:377-384`).
- **Diff**: per-base divisor differs slightly. Upstream divides by the
  cumulative mass of the accessible ratemap inside accepted windows.
  Fork divides by `tree_covered_accessible_bp`, which iterates over
  `ts.trees()` and counts only intervals where `tree.num_edges > 0`,
  intersected with `acc_intervals`. Both compute the "tree-covered &
  accessible" bp; in fully-covered ARGs they should match. They will
  diverge if SINGER ever produces non-edge intervals not already masked
  by the inaccessibility map. Default `span_normalise` is also handled
  differently (upstream sets it explicitly to False, fork omits — tskit
  default is True for `allele_frequency_spectrum`, which means fork's
  numerator is already site-density per bp before dividing again by
  `total_accessible` — see Open Question §4).

### 2.4 Mutational load per sample

- **Upstream** (`utils.py:126-143`, called from
  `tree_statistics.py:57-58, 88-89`): sweeps tree-by-tree, digitises
  `ts.sites_position[ts.mutations_site]` to find a window per mutation,
  uses `tree.samples(m.node)` to find descendant samples, increments per
  `(window, sample)`. **Squeezes single-window output to 1-D.** Final
  normalisation: `observed_load /= np.sum(accessible_bp)` — divides by
  total accessible bp across the whole chromosome, not per window. The
  same is done for `expected_load = mutational_load(sim_ts) / np.sum(...)`
  (the simulated TS).
- **Fork** (`argtest_common.py:107-158`): adds two features —
  (a) `windows` parameter (returns 2-D matrix), (b) `remove_intervals`
  + `name_to_nodes` for per-individual drop segments. The drop logic
  uses `ts.keep_intervals` + `ts.simplify(samples=keep)` then maps node
  ids back via a reverse map to accumulate load on the original sample
  ids. Not enabled by default. When `windows is None`, returns 1-D —
  matches upstream behaviour.
- **Validation plots** (`validation_plots_from_ts.py:338-340`):
  `load = mutational_load(ts) / total_accessible`. **No separate
  "expected" trace from a sim-mutations TS** — the fork dropped that path
  (commit `c2192cf`). Upstream's mutational-load plot (`diagnostics.py:161-172`)
  shows observed-vs-expected; fork's plot shows observed only.
- **Mutload-mask path** (`scripts/mutload_masks.py:138-151`): per-replicate,
  builds `windows` (bp or SNP), computes `load = mutational_load(ts,
  windows=...)`, simulates a single sim TS via `msprime.sim_mutations` with
  deterministic seed, computes `expected = mutational_load(sim_ts,
  windows=...)`, declares an outlier where `load > (1+cutoff)*expected` or
  `load < (1-cutoff)*expected`. Aggregated to individuals via
  `aggregate_by_individual`. Outputs BEDs of (chrom, start, end,
  outlier_names, observed, expected). **This is fork-only.**
- **Diff**: upstream's "expected_load" is one number per sample averaged
  over the whole chromosome; fork's expected is per-window per-individual
  for masking purposes, and fork's validation plots no longer overlay an
  expected trace at all.

### 2.5 Pair coalescence rates and PDF

- **Upstream** (`coalescence_rates.py`): builds time grid
  `np.logspace(log_min, log_max, num_intervals + 1)` with `0` and `inf`
  prepended/appended (`:32-34`). Calls `ts.pair_coalescence_counts(...,
  pair_normalise=True)` and `ts.pair_coalescence_rates(...)`. Tail
  cutoff via cumulative survival. Crucially, runs
  `collapse_masked_intervals(ts, accessible)` first (`:39`) — that
  reindexes the TS to compress masked segments out of coordinates so
  partial trees become full trees on the collapsed sequence.
- **Fork** (`coalescence_ne_plots_from_ts.py`): builds time grid two
  ways — explicit edges from a file (`load_time_windows`,
  `:185-199`), or quantile-based bins from
  `ts.pair_coalescence_quantiles(...)` averaged across post-burnin
  replicates (`compute_quantile_time_windows`, `:202-235`). Pads with 0
  and inf identical to upstream. Tail cutoff identical
  (`compute_pair_coal`, `:248-253`). **Does NOT collapse masked
  intervals** — operates directly on the TS as written by the trim
  step, so partial trees produce biased rates as documented in
  `:1-11` and the warning in commit `c41ccc1`.
- **Plot** of PDF: upstream plots the raw bin PMF (`plot_coalescence_rates.py:85-104`).
  Fork converts PMF to density on the log scale by dividing each bin's
  PMF by `log(right/left)` (`coalescence_ne_plots_from_ts.py:413-415`).
- **Ne plot**: upstream has no Ne plot. Fork plots `1/(2 rate)` and
  optionally a year-marker dashed line (`:447-464`).
- **Diff**: the partial-tree handling is the central methodological
  divergence. On a fully-spanned ARG (no trim step) both should agree
  up to time-grid choice. After fork's `trim_regions_single` /
  `trim_samples_single`, fork's rates are biased (per its own warning);
  upstream's `collapse_masked_intervals` would have rebuilt a fully-spanned
  TS but is not invoked by fork's `coalescence_ne_plots_from_ts.py`.

### 2.6 Window construction

- **Upstream**: windows are picked once in `chunk_chromosomes.py:442-444`
  as `np.linspace(0, bitmask.size, num_stats_windows + 1).astype(int)`,
  pickled to disk, and reused unchanged by `tree_statistics.py` and
  `diagnostics.py`. `windows.rate ∈ {0,1}` carries an accept/reject
  per window from upstream filter logic.
- **Fork** (`validation_plots_from_ts.py:191-205, 308-336`): each
  replicate gets its own windows
  (`np.arange(0, ts.sequence_length + window_size, window_size)` clamped
  at the right edge to satisfy tskit's `windows[-1] == sequence_length`
  requirement), and accessible-bp denominators are recomputed per
  replicate from `kept_intervals` metadata. After collecting,
  `min_n_windows = min(len(a) for a in site_div_vals)` truncates to a
  shared length.
- **Diff**: upstream's window grid is tied to the chunked bitmask
  coordinate system from VCF processing; fork's is purely tree-sequence
  coordinate space. If fork's `trim_regions_single` rewrites coordinates,
  both window grid and `kept_intervals` interpretation depend on what
  step5 leaves behind — see Open Question §4.

### 2.7 Treatment of the mutation rate map

- **Upstream**: pickled `*.mut_rate.p` is the `adjusted_mu` from chunked
  pipeline output and is passed straight into `msprime.sim_mutations`
  in `tree_statistics.py:62-66`.
- **Fork**: `argtest_common.resolve_mu_rate` resolves rate from
  (in order) ts.metadata ratemap, sibling `*.mut_rate.p`, scalar
  fallback (`argtest_common.py:372-399`). Used by both validation_plots
  and mutload_masks. `ratemap_from_metadata` (`:364-369`) lets the
  ratemap travel inside the TS instead of as a sidecar — fork-only.

### 2.8 Treatment of masked / kept intervals

- **Upstream** uses two ratemaps: `inaccessible` (from chunk step) and
  `extract_accessible_ratemap(trees)` (from current trees, derived from
  edge-coverage). The product is the mask used for cumulative-mass per
  window and for `collapse_masked_intervals` in coalescence rates.
- **Fork** uses three sources:
  (1) `ts.metadata["kept_intervals"]` written by `trim_regions_single`
  (the offset-merge fix in commit `781d1c9` is the latest behaviour),
  (2) accessible intervals from the mu ratemap (`accessible_intervals_from_mu`),
  (3) full sequence as fallback. `tree_covered_accessible_bp` further
  intersects with edge-covered tree spans — same idea as upstream's
  `extract_accessible_ratemap` AND step.

---

## 3. Fork-only behaviors

- **Mutload outlier masking pipeline** (`scripts/mutload_masks.py`,
  `scripts/mutload_summary.py`): a per-window per-individual sim-based
  expected load test that emits BEDs and an HTML report. Not present
  upstream.
- **Demes-graph + msprime posterior predictive simulations**
  (`scripts/coalescence_ne_plots_from_ts.py:275-351`): builds a
  piecewise-constant Demes graph from the inferred Ne(t), simulates
  multiple 1 Mb chunks with mutations, and writes per-window pi/Tajima's D
  plus per-sim SFS to TSV. These TSVs feed back into validation_plots
  via `--sim` and `--sim-sfs` to produce density-overlay plots.
- **Quantile-based equal-event time bins** for coalescence rates
  (`compute_quantile_time_windows`, `:202-235`). Upstream uses fixed
  log-spaced bins.
- **Effective-population-size plot with optional year marker**
  (`coalescence_ne_plots_from_ts.py:447-464`).
- **`--compare` overlay** in validation_plots_from_ts.py for pre-vs-post
  pipeline comparison (`:484-493`).
- **Folded AND unfolded SFS in one run**
  (`validation_plots_from_ts.py:347-356, 690-725`).
- **kept_intervals metadata** propagation through `trim_regions_single`
  (`argtest_common.merge_intervals`, `:341-353`; commit `781d1c9`
  offset-fix).
- **mutational_load with windows + remove_intervals + simplify**
  (`argtest_common.py:107-158`) — used to drop individuals on
  per-segment basis without rebuilding the TS.

---

## 4. Upstream-only behaviors

- **VCF-side observed stats**: upstream computes pi, Tajima's D, SFS
  directly from the input VCF via scikit-allel
  (`chunk_chromosomes.py:454-469`). The fork has no VCF input pipeline
  — its inputs are tree sequences from upstream of this fork. Therefore
  fork's "observed" pi/D/SFS comes from `ts.diversity(mode="site")`
  applied to whatever the SINGER-output TS contains. They should
  *typically* agree but are not equivalent: any post-SINGER topology
  manipulation that changes site count (collapse, simplify, mask) will
  alter fork's "observed" stats. Upstream's VCF stats are immutable
  ground truth.
- **Stratified divergence and per-stratum SFS**
  (`tree_statistics.py:103-140`, `chunk_chromosomes.py:655-708`,
  `diagnostics.py:192-294`): fork has no analogue.
- **Cross-coalescence rates / PDF between strata**
  (`coalescence_rates.py:55-95`, `plot_coalescence_rates.py:37-81`):
  fork has no analogue.
- **`repolarised` and `multimapped` summaries**
  (`tree_statistics.py:50-51`): rely on `s.metadata["flipped"]` and
  `s.metadata["omitted"]` populated by `utils.absorb_mutations_above_root`
  / chunking. Fork doesn't carry those fields and doesn't plot them.
- **`collapse_masked_intervals` before coalescence-rate computation**
  (`coalescence_rates.py:39`, `utils.py` upstream version): collapses
  masked sequence out of the TS so trees become fully-spanned, sidestepping
  the partial-tree bias. Fork explicitly omits this step and warns
  about the resulting bias (commit `c41ccc1`,
  `coalescence_ne_plots_from_ts.py:1-11`).
- **`absorb_mutations_above_root`** (`utils.py:79-123`) — preprocessing
  step that records flipped ancestral states. Fork's TS inputs are
  assumed to already be in the right orientation; no analogue.
- **`find_genealogical_gaps`** (`utils.py:158-191`): explicit detection
  of gaps in edge coverage that span uninterrupted intervals. Fork's
  `tree_covered_accessible_bp` only sums covered bp — it does not
  emit a gap-interval list.
- **Site `omitted` filtering before stats**
  (`tree_statistics.py:46-47`): upstream calls
  `trees.delete_sites(np.flatnonzero(omitted))` so dating-omitted sites
  don't contribute to repolarised/multimapped counts. Fork has no such
  step because it doesn't carry the omitted flag.
- **Single-config polarisation choice** (`POLARISED` in upstream Snakefile)
  applies to the entire run; fork emits both folded and unfolded by
  default (no global polarisation setting).

---

## 5. Open questions / things to verify empirically

1. **Default tskit `span_normalise` for `allele_frequency_spectrum`.**
   RESOLVED 2026-05-14 — **fork is double-normalising.** tskit 1.0.2 (the
   active env's installed version) defaults `span_normalise=True` for
   `mode="site"` (verified empirically: a 1 Mb msprime sim returned
   `n_sites/seq_len` from the default call). The earlier speculation that
   site-mode defaults to False was wrong; the default is True for ALL modes
   in tskit ≥ 1.0. Therefore:
   - Fork's `validation_plots_from_ts.py:347-352, 377-384` computes
     `(raw_count / seq_length) / total_accessible` ≈ `raw_count / (seq_length × total_accessible)`.
     For a typical post-trim run where `total_accessible ≈ seq_length`, the
     SFS values are off by a factor of ~`seq_length` (10^6–10^8) relative
     to the intended "variants / base" semantics on the y-axis label.
   - Both observed and sim-branch series in `frequency-spectrum-*.png` are
     double-normalized by the same factor, so they overlay correctly on
     each other but at the wrong absolute scale. The `--sim-sfs` overlay
     path is worse: `coalescence_ne_plots_from_ts.py:330-331` writes
     single-normalized values to `sim-sfs.tsv`, then
     `validation_plots_from_ts.py:783-786` overlays them against
     double-normalized observed — sim and observed will be vertically
     offset by ~`total_accessible` on log scale.
   - **Fix:** add `span_normalise=False` to the four AFS call sites in
     `validation_plots_from_ts.py:348, 351, 378, 382`. The two
     `coalescence_ne_plots_from_ts.py:330, 331` calls in
     `simulate_window_stats_from_ne` are **correct as-is**: the Demes sim
     has `sequence_length=1_000_000` with no masking, so the default
     `span_normalise=True` divides by 1 Mb = the accessible bp of the
     sim, producing "variants per accessible bp" — the same units the
     fixed validation_plots overlay produces. Applied 2026-05-14.

2. **Per-window denominator equivalence.** RESOLVED 2026-05-14 — fork's
   `kept_intervals` does NOT account for tree-edge gaps.
   `trim_regions_single.py:38` writes `kept_intervals` as the pure
   complement of the mask BED (`complement_intervals(masked, seq_length)`)
   with `ts.keep_intervals(keep, simplify=False)`, i.e. it never intersects
   with edge coverage. Consequence: any window whose accepted-region span
   contains a tree-edge-empty subinterval will have a denominator larger
   than the upstream equivalent, biasing per-window pi slightly LOW (more
   bp in denominator than there are edge-covered bp). In practice
   `simplify=False` keeps the kept-interval edges intact, so the effect
   should be small unless SINGER produces internal edge gaps inside the
   kept regions — empirical magnitude not measured here. The SFS divisor
   path is unaffected because `validation_plots_from_ts.py:327` uses
   `tree_covered_accessible_bp`, which DOES intersect with edge coverage.

3. **Behaviour after `collapse_masked_intervals` is *not* called.**
   STRUCTURALLY CONFIRMED 2026-05-14. `collapse_masked_intervals` is
   defined in `argtest_common.py:566` and used by the LEGACY multi-mask
   `trim_regions.py:125`, but the ACTIVE pipeline step
   (`trim_regions_single.py:38`) uses `ts.keep_intervals(keep, simplify=False)`
   without any coordinate collapse. After the subsequent `trim_samples_single`
   step calls `remove_ancestry` (`trim_samples.py:89`), trees inside
   `kept_intervals` lose edges for the excised individuals — i.e. they
   become partial trees within their span. `coalescence_ne_plots_from_ts.py`
   reads this TS as-is, so partial trees ARE present whenever
   `trim_samples_single` removed anyone. Magnitude on a real trimmed run
   is still empirical work; subsumed by the partial-trees Ne-fix todo.

4. **Mutload-mask seed convention.** RESOLVED 2026-05-14 — plumbing
   distinguishes seeds in production. `Snakefile:83 mutload_seed_for`
   computes `(sha1(MUTLOAD_RANDOM_SEED|chrom|rep)[:8] % (2**31-1)) + 1`
   and passes it as `--random-seed {params.seed}` (`Snakefile:325`), so
   each (chrom, rep) pair gets a distinct deterministic seed.
   `mutload_masks.py`'s default of `--random-seed 1` only affects manual
   invocations outside Snakemake — those will collapse all replicates to
   the same draw. Note (independent of upstream comparison): fork's
   `--sim-branch` overlay in `validation_plots_from_ts.py:363` uses
   `random_seed=1 + int(rep_id) * 1000`, matching upstream's convention,
   so the two seed schemes coexist in the same pipeline.

5. **Sample-set / stratum extraction.** Upstream pulls strata from
   `ts.populations()` metadata; fork pulls names from
   `individual.metadata["id"]` and groups by lineage prefix
   (`mutload_summary.py:67-73`). These will produce different groupings
   if the input TS uses populations rather than per-individual id
   metadata.

6. **SFS divisor in `tree_covered_accessible_bp`.** Fork's denominator
   is "accessible bp covered by any non-empty tree". If a tree has
   edges but covers only a strict subset of samples (post-trim), it is
   still counted as fully covered — the SFS bins for that interval
   may be skewed because absent samples don't contribute mutations.
   Compare against upstream's behaviour, which avoids this by
   collapsing masked intervals out.

7. **Window grid alignment under `--compare`.** Validation plots
   truncate to `min_n_windows` across the primary set, then run
   `--compare` with its own collection. The two coordinate systems
   may not be aligned if pre- and post-pipeline TSs have different
   sequence lengths (e.g. after `collapse_masked_intervals` were
   ever applied — currently it is not). Skyline plots use
   `pri["coord"]` and `cmp["coord"]` independently which is fine for
   visualisation but means a side-by-side window-i comparison is
   not coordinate-matched.

---

## 6. File:line reference index

Fork:
- `/home/jri/src/argtest/scripts/validation_plots_from_ts.py:268-455` —
  `collect_stats` core; site-mode pi/D/SFS, sim-branch overlay.
- `/home/jri/src/argtest/scripts/validation_plots_from_ts.py:540-549` —
  mutational-load trace plot (no expected line).
- `/home/jri/src/argtest/scripts/coalescence_ne_plots_from_ts.py:248-253` —
  pair-coal rates with tail cutoff.
- `/home/jri/src/argtest/scripts/coalescence_ne_plots_from_ts.py:275-351` —
  Demes graph + msprime simulations.
- `/home/jri/src/argtest/scripts/coalescence_ne_plots_from_ts.py:447-464` —
  Ne plot.
- `/home/jri/src/argtest/scripts/mutload_masks.py:115-185` —
  expected-load test, BED emission.
- `/home/jri/src/argtest/scripts/mutload_summary.py:173-179` —
  HTML summary expected-load.
- `/home/jri/src/argtest/scripts/argtest_common.py:107-158` —
  fork's `mutational_load` (windowed + drop sets).
- `/home/jri/src/argtest/scripts/argtest_common.py:372-399` —
  `resolve_mu_rate` rate fallback chain.
- `/home/jri/src/argtest/scripts/argtest_common.py:418-451` —
  `accessible_intervals_from_mu`, `tree_covered_accessible_bp`.
- `/home/jri/src/argtest/scripts/argtest_common.py:566-594` —
  `collapse_masked_intervals` (defined but unused by stat scripts;
  fork's coalescence script does NOT call it).

Upstream (raw URLs):
- `workflow/scripts/tree_statistics.py:31-100` — sim-on-trees stat block.
- `workflow/scripts/tree_statistics.py:103-142` — stratified stats.
- `workflow/scripts/diagnostics.py:30-187` — main plot path.
- `workflow/scripts/diagnostics.py:192-294` — stratified plots.
- `workflow/scripts/coalescence_rates.py:32-95` — log-spaced bins,
  collapse-masked-intervals, pair + cross.
- `workflow/scripts/plot_coalescence_rates.py:37-104` — plotting.
- `workflow/scripts/utils.py:126-143` — `mutational_load` (no windows,
  no drop logic).
- `workflow/scripts/utils.py:79-123` — `absorb_mutations_above_root`.
- `workflow/scripts/utils.py:158-191` — `find_genealogical_gaps`.
- `workflow/scripts/utils.py:194-217` — `multiply_ratemaps`,
  `extract_accessible_ratemap`.
- `workflow/scripts/chunk_chromosomes.py:441-475` — VCF site-mode stats.
- `workflow/scripts/chunk_chromosomes.py:642-708` — VCF stat output,
  stratified VCF stats.

---

## 7. Headline summary

The two pipelines are **methodologically aligned on the sim-on-trees
posterior predictive path** for pi, Tajima's D, SFS, and mutational
load: both run `msprime.sim_mutations(ts, rate=adjusted_mu, keep=False)`
and compute `mode="site"` stats on the result, with a per-base
denominator that combines mutation-rate accessibility and tree-edge
coverage. The fork dropped only the *unsimulated* branch-mode shortcut
(commit `c2192cf`).

They are **meaningfully diverged** in three places:

- Fork has no VCF-side stat path — "observed" comes from the
  SINGER-output TS, not the input data, which is fine for posterior
  predictive checks but breaks the upstream observed-vs-expected
  contract on the mutational-load and skyline plots.
- Fork's `coalescence_ne_plots_from_ts.py` does **not** call
  `collapse_masked_intervals` and explicitly warns its rates/Ne are
  biased on partial trees. Upstream's `coalescence_rates.py` always
  collapses first.
- Fork adds a new mutload-outlier filtering step
  (`scripts/mutload_masks.py`) and a Demes-graph posterior-predictive
  simulation pipeline (`coalescence_ne_plots_from_ts.py:275-351`),
  neither of which exist upstream.

Stratified analyses (cross-coalescence, per-stratum SFS, divergence
matrix) and the `repolarised`/`multimapped`/`omitted` plumbing are
upstream-only.
