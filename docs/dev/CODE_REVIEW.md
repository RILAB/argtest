# Code review: correctness, documentation, efficiency, and simplification

Reviewed 2026-08-09. Scope: `Snakefile`, `scripts/`, tests, configuration, and user-facing documentation. No implementation files were changed.

## Executive summary

The workflow is thoughtfully structured and has unusually good coverage of tree-sequence edge cases. The main opportunity is not a broad rewrite; it is to make the supported path narrower. There are two likely pipeline failures involving missing output directories, one sample-name matching bug, and a handful of silent fallbacks that can turn bad inputs into plausible-looking outputs. Separately, roughly 200–300 lines appear removable immediately because they implement unused or superseded paths. More substantial simplification is possible if the project explicitly assumes nonempty ARGs, positive chromosome lengths, leaf sample nodes, uniform ploidy, valid pipeline-produced BED files, and current pinned versions of `tskit`/`msprime`.

Recommended order:

1. Fix the two output-directory bugs and suffix stripping.
2. Replace silent input recovery with fail-fast behavior where a wrong result is worse than a failed job.
3. Delete unused tree-trimming helpers and `trim_samples_chunked.py`.
4. Declare the normal-input assumptions above and remove defensive branches that only support impossible or unsupported ARGs.
5. Consolidate mutation-load and mutation-map resolution code before doing further performance work.

## Correctness findings

### High: step 2 can fail when its chromosome output directory does not already exist

`scripts/find_low_access_regions.py:103-112` writes `out_path` without first creating `out_path.parent`. The Snakemake output is nested (`step2_low_access/<chrom>/<chrom>.low_access.bed`; `Snakefile:397-398`), and no earlier rule is guaranteed to create that chromosome directory. The log directory is only created after the BED write (`find_low_access_regions.py:119-123`), so it does not help.

Recommendation: create `out_path.parent` immediately before writing, consistent with the other pipeline scripts. Add a CLI test whose output parent does not exist.

### High: step 5 has the same missing-parent failure

`scripts/trim_samples.py:363-369` creates a directory only for its default output. With explicit `--out`—which is how `Snakefile:488-497` invokes it—the parent is not created before `dump_ts`. On a fresh run, `step5_trimmed_samples/<chrom>/` is not guaranteed to exist.

Recommendation: always run `out_path.parent.mkdir(parents=True, exist_ok=True)` after resolving `out_path`. This should be tested through the CLI, not only by calling `trim_samples_single_pass` directly.

### Medium: `suffix_to_strip` removes text anywhere in a name, not just a suffix

`scripts/argtest_common.py:61` uses `nm.replace(suffix_to_strip, "")`. For example, stripping `_A` changes `sample_A_old` even though `_A` is not a suffix, and repeated occurrences are all removed. The CLI and README consistently call this operation “suffix” stripping.

Recommendation: use `str.removesuffix`, or rename the option to make global replacement explicit. Prefer `removesuffix`; it matches user expectations and is simpler.

### Medium: missing or malformed mask inputs are silently accepted

`scripts/combine_remove_masks.py:33-49` treats a missing BED as empty and silently skips rows with fewer than three fields. In the pipeline all three inputs are declared Snakemake dependencies, so missing files should already be impossible; this recovery mainly hides standalone invocation errors and corrupted outputs. A partial combined mask is a successful-looking but scientifically wrong result.

Recommendation: remove both recovery branches and raise on missing/malformed inputs, as `load_remove_intervals` already does. This is less code and safer.

### Medium: broad exception handling can discard valid mutation-map errors

`scripts/trim_regions_single.py:47-54` catches every exception while locating and unpickling a mutation map and then proceeds with no rate metadata. This conflates “no sibling map exists” with ambiguity, corruption, incompatible pickle contents, and programming errors. `validation_plots_from_ts.py:171-196` has the same pattern twice. Results can therefore fall back to whole-sequence accessibility without making the reason visible.

Recommendation: catch `FileNotFoundError` only. Let ambiguity and corrupt data fail. Because the pipeline already embeds rate metadata, consider removing sibling-pickle fallback from validation entirely.

### Medium: duplicate individual IDs overwrite rather than combine nodes

`scripts/argtest_common.py:98-105` assigns `mapping[nm] = nodes`; duplicate metadata IDs cause earlier nodes to disappear. Elsewhere, `sample_names` and `aggregate_by_individual` explicitly group duplicate names, so behavior is inconsistent.

Recommendation: either enforce unique individual IDs once with a clear error (simplest) or extend the existing list. Given the requested tolerance for unlikely inputs, enforcing uniqueness is preferable.

### Low: production preconditions use `assert`

`scripts/argtest_common.py:116,624,633-634,682,688-690` uses assertions for input and API preconditions. Python removes these under `-O`, so invalid windows or incompatible rate maps can proceed and fail later in less obvious ways.

Recommendation: for truly supported user errors, use one explicit `ValueError`. For invariants guaranteed by the pipeline and pinned dependencies, remove the checks rather than retaining a large defensive layer.

## Efficiency findings

### High impact: expected mutation load is estimated with a full mutation simulation

`scripts/mutload_masks.py:118-123` calls `msprime.sim_mutations` once per chromosome/replicate and then traverses every simulated mutation and descendant set through `mutational_load`. This is likely step 3's dominant cost and introduces Monte Carlo noise into a quantity described as “expected” load. `mutload_summary.py:161-167` duplicates the same work.

The expected derived count can in principle be computed from branch lengths and the mutation rate map without materializing a second mutated tree sequence. That would avoid allocating all simulated sites/mutations, reduce runtime and memory, and make results deterministic. This is the one optimization worth benchmarking before micro-optimizing Python loops.

Recommendation: implement a branch-length expectation using tree statistics and rate-map overlap, then compare it numerically and by runtime against the current seeded simulation on representative ARGs. If the stochastic simulation is scientifically intentional, rename it from “expected” and document its variance.

### Medium impact: mutation counting builds Python sample lists for every mutation

`scripts/argtest_common.py:139-157` iterates trees, sites, and mutations, materializes `list(tree.samples(m.node))`, optionally remaps it with another Python list, then performs indexed addition. Runtime grows with mutation count times descendant count.

Recommendation: first check whether the pinned `tskit` can express this as a sample-count/general-statistic operation. If not, cache descendant arrays per `(tree.index, mutation.node)` when recurrent mutations share nodes, and avoid `list` conversion. Benchmark before changing: this is central scientific code, so clarity may be worth more than a modest speedup.

### Medium impact: validation loads each tree sequence repeatedly

`validation_plots_from_ts.py` first loads the first file during mutation-map discovery (`:171-180`), then loads every file again in `collect_stats` (`:332-334`). The Snakemake rule runs the entire program twice—cleaned and original (`Snakefile:762-778`). For large compressed ARGs, decompression and table allocation are material.

Recommendation: pass an already-loaded first tree sequence into collection, and consider one invocation that processes original/cleaned pairs. Do not merge these paths if it makes SLURM memory retention worse; measure peak RSS as well as elapsed time.

### Low/conditional: interval overlap assumes sorted, nonoverlapping accessibility intervals

`tree_covered_accessible_bp` and `overlap_lengths` use efficient two-pointer scans (`argtest_common.py:493-516,569-590`), but their contract is only in comments and callers may pass metadata directly. Sorting/merging defensively inside each call would cost time and obscure the contract.

Recommendation: explicitly document “sorted, disjoint intervals” as a required internal invariant and normalize once when ingesting external metadata. This is a good place to accept less protection for simpler, faster hot loops.

## Over-engineering and simplification opportunities

### Delete unused production helpers

The following functions have no production callers; they are referenced only by tests or by one another:

- `build_segments_with_drop_nodes`, `trim_ts_by_intervals`, and `assert_sample_ids_preserved` in `argtest_common.py:184-273`.
- `collapse_masked_and_low_access_windows` and `build_shared_mask` in `argtest_common.py:662-697`.

The first group appears to be the superseded implementation of sample/interval trimming. The second group appears to be a superseded combined masking path. Removing both groups would eliminate about 140 lines and several concepts from the shared module. Delete their tests at the same time; tests for dead code are maintenance, not coverage.

### Delete or formally support `trim_samples_chunked.py`

`scripts/trim_samples_chunked.py` is not called by the Snakefile, documented in the README, imported by other code, or tested directly. It duplicates the most complex part of `trim_samples.py`, but its surgery path has already diverged: unlike the primary implementation it does not call `_filter_trimmed_sample_mutations`. Keeping an unused alternative algorithm is high-risk over-engineering.

Recommendation: delete it. If real datasets demonstrably require chunking, make it the sole implementation and add equivalence tests before advertising it.

### Remove the internal-sample-node compatibility branch

`trim_samples.py:127-223` spends nearly 100 lines preserving mutations on target sample nodes that are also internal parents. Its own comment says pipeline sample nodes are leaves (`:135-136`). This is exactly the kind of highly unlikely scenario the review was asked to challenge.

Recommendation: declare “samples must be leaves” as an input contract and validate it once near loading if desired. Then simplify mutation removal to the leaf-only case. Do not keep a second algorithm embedded in the hot path for unsupported ARG structure.

### Narrow VCF export to the data model the pipeline actually uses

`export_vcf.py:53-80` supports missing individuals and mixed ploidy by silently falling back to one haploid VCF column per node. That changes sample identity and output semantics rather than merely improving robustness. The README advertises ploidy-aware individual output.

Recommendation: if pipeline inputs always assign samples to uniformly ploidy individuals, enforce that and fail otherwise. This removes the fallback and prevents surprising VCFs. Keep mixed-ploidy support only if it is a real input requirement and then implement it explicitly rather than silently changing representation.

### Stop supporting multiple mutation-rate discovery conventions

`infer_mu_base`/`infer_mu_path` (`argtest_common.py:519-566`) searches multiple filename stems across two directories, then performs a second glob-based fuzzy search. This is convenient but complex and can create ambiguity. The pipeline already has explicit config and embeds the rate map after trimming.

Recommendation: choose one source of truth: explicit config path at ingestion, then embedded metadata downstream. Removing heuristic discovery would simplify setup, error messages, staging in step 6, and several broad exception handlers.

### Consolidate duplicated window and expected-load helpers

`build_bp_windows`, `build_snp_windows`, and simulated expected-load logic are independently implemented in `mutload_masks.py` and `mutload_summary.py`. Validation has a third base-pair window builder. The versions already differ in where validation occurs.

Recommendation: put these small scientific primitives in one focused module. Keep HTML/report rendering out of `argtest_common.py`; that module is already a 697-line mixture of I/O, names, coordinates, rate maps, statistics, and legacy trimming.

### Reduce compatibility probing for the pinned dependency set

`validate_trimmed_ts` (`argtest_common.py:275-303`) probes several historical `tskit` APIs with nested `hasattr`, `TypeError`, and broad exception branches. `environment.yml` pins a specific Git commit, so multi-version compatibility is not a current requirement.

Recommendation: call the API supported by the pinned version directly. This removes about 25 lines and makes a changed dependency fail clearly. If broad library compatibility becomes a release goal, test a version matrix instead of detecting APIs at runtime.

### Treat report-only empty-data presentation as optional, not core complexity

`validation_plots_from_ts.py:216-247` and `pipeline_summary.py:237-277` contain multiple helpers for all-NaN, zero, empty, one-replicate, and zero-length presentation. These branches prevent plots from crashing on empty-mutation ARGs and summaries from dividing by zero, but they add noticeable surface area.

Recommendation: retain protection against missing optional plot files, which is operationally useful. Drop support for zero-length chromosomes, zero-sample ARGs, and entirely empty statistics if those inputs are outside the project contract. A single early validation error is simpler than making every plotting helper tolerate them.

## Documentation findings

### README gives the wrong default for `trim_samples.py`

`README.md:317-319` says the default is `results/<ts_stem>_trimmed.tsz`; the code writes `<input-parent>/trimmed/<ts_stem>_trimmed.tsz` (`trim_samples.py:339,366-368`). The argparse help correctly describes `results/...` too (`trim_samples.py:75`), so both user-facing descriptions disagree with behavior.

Recommendation: choose one behavior and make code, `--help`, and README agree.

### Workflow description names the wrong region-trimming script

`README.md:79` says the automated workflow trims regions with `trim_regions.py`, while the Snakefile calls `trim_regions_single.py` (`Snakefile:443-467`). The two scripts have materially different semantics: the former compacts coordinates, while the pipeline version preserves them.

Recommendation: update the workflow narrative to name `trim_regions_single.py` and emphasize coordinate preservation there. Consider moving `trim_regions.py` into an `extras/` directory or deleting it if coordinate compaction is no longer supported.

### The README's logging promise is too broad

`README.md:322` suggests detailed errors are available in corresponding logs. Several scripts only write a success summary after work completes, and broad exception handlers may suppress logging errors. Snakemake shell redirection captures some jobs, but standalone behavior is inconsistent.

Recommendation: state that Snakemake captures stdout/stderr under `logs/`, and describe per-script summary logs separately. Avoid implying structured error logs exist.

### Missing architectural boundary documentation

There is detailed operational documentation, but no short statement of supported invariants: positive sequence length, at least one sample, leaf sample nodes, unique individual IDs, sorted valid rate maps, and consistent individual/ploidy structure across replicates.

Recommendation: add a “Supported input assumptions” section. This enables removal of many scattered defensive checks while making the reduced robustness an explicit design decision rather than an accident.

## Testing and maintainability

- Tests are numerous and target the complex algorithms, which is a strength. However, they emphasize helper-level edge cases more than fresh-directory CLI execution; that is why both missing-parent bugs can survive.
- Add one small end-to-end CLI test per pipeline script using a new temporary output tree. These tests give more value than testing retired helpers.
- The configured environment could not run the suite during this review: `/opt/anaconda3/envs/argtest` activated with Python 3.14.3, but neither `pytest` nor `snakemake` was installed, despite both appearing in `environment.yml`. `compileall` completed successfully. Rebuilding or synchronizing the environment is necessary before relying on CI-like local validation.
- The current untracked files (`mutload.html`, `time_bins.txt`, and `time_bins.adjusted.txt`) predated this review and were not inspected as source or modified.

## Suggested reduced support contract

Adopting the following contract would enable the largest safe simplification:

- Every chromosome has positive length and every ARG has at least one sample.
- Sample nodes are leaves.
- Every sample belongs to exactly one individual; individual IDs are unique; ploidy is uniform within a tree sequence.
- BED files are pipeline-produced or otherwise syntactically valid; missing inputs are errors.
- Mutation maps enter through one explicit configured path and are embedded in tree-sequence metadata for downstream steps.
- Only the dependency versions in `environment.yml` are supported.
- Auxiliary visualization may fail clearly on an ARG with no mutations rather than manufacturing an empty plot.

This trades graceful behavior on exotic inputs for a smaller, faster, and more auditable scientific pipeline. It does not require weakening checks for conditions that commonly indicate wrong biological inputs, such as chromosome-label mismatches, inconsistent replicate sets, incompatible sequence lengths, or ambiguous mutation maps; those should remain hard failures.
