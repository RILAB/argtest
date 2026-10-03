# Simplification and reliability implementation plan

Derived from `CODE_REVIEW.md` and follow-up decisions on 2026-08-09.

## Decisions and scope

### Accepted design decisions

- Keep mutation simulation as the basis of expected mutational load.
- Keep the term **simulation-based expected load**. This is accurate enough as long as documentation states that it is estimated from one seeded mutation simulation, rather than an analytic expectation.
- Assume uniform ploidy within each tree sequence.
- Continue supporting multiple sample nodes per individual. Diploid individuals therefore contribute two sample nodes, with their loads pooled for per-individual analysis and their two alleles grouped in VCF output.
- Rename the current suffix-removal setting because it removes every occurrence of a string, not only a terminal suffix.
- Prefer explicit supported-input contracts and fail-fast errors over recovery for exotic or malformed inputs.
- Keep checks for biologically plausible setup mistakes such as chromosome mismatches and inconsistent replicate sets when those checks exist on live, supported pipeline paths. This does not commit to preserving checks that belong solely to deleted auxiliary tools. Mutation-map resolution will no longer raise a glob-derived ambiguity error: Phase 4 replaces fuzzy discovery with exact candidate lookups and a `FileNotFoundError` listing the paths tried.
- Preserve the current recurrent/back-mutation counting semantics. Do not replace `mutational_load` with a `ts.variants()` genotype-presence sweep merely for speed; that would be a scientific behavior change.
- Use `kept_intervals` metadata first and the embedded mutation map second as the common definition of accessible sequence in summaries and validation.

### Naming decision

Rename:

- Config: `suffix_to_strip` → `name_substring_to_remove`
- CLI: `--suffix-to-strip` → `--name-substring-to-remove`
- Python parameters/attributes: `suffix_to_strip` → `name_substring_to_remove`

The operation remains global substring removal with `name.replace(name_substring_to_remove, "")`. Do not change it to suffix-only behavior.

This should be a deliberate breaking rename rather than maintaining permanent aliases. During configuration parsing, detect the old `suffix_to_strip` key and raise a short migration error naming `name_substring_to_remove`; this protects users from silently getting the default while keeping the implementation simple. The old CLI spelling may likewise be rejected naturally by `argparse`; no compatibility alias is required.

## Phase 1: fix concrete workflow failures

### 1. Create output parents consistently

Files:

- `scripts/find_low_access_regions.py`
- `scripts/trim_samples.py`
- Corresponding tests

Work:

1. Resolve the final output path.
2. Create `out_path.parent` before the first write.
3. Keep log-directory creation separate and explicit.
4. Add CLI-level tests using output paths whose parent directories do not yet exist.

Acceptance criteria:

- Standalone step-2 CLI invocation succeeds when the requested output parent does not exist.
- Standalone `trim_samples.py` invocation succeeds when an explicit `--out` parent does not exist.
- Existing output files retain the same contents and naming.

Snakemake prepares declared output parent directories, so this is standalone CLI reliability rather than a pipeline-DAG failure. Keep the fix small and do not represent it as necessary for normal Snakemake execution.

### 2. Make mask combination fail fast

Files:

- `scripts/combine_remove_masks.py`
- `tests/test_combine_remove_masks.py`

Work:

- Raise `FileNotFoundError` for a missing input BED.
- Raise `ValueError` with file and line context for a nonblank, noncomment row with fewer than three fields.
- Continue accepting empty BED files as valid pipeline outputs.
- Continue ignoring blank lines and comments.
- Apply integral rounding where float windows are first emitted. In `mutload_masks.py`, write `floor(windows[w])` and `ceil(windows[w + 1])` for both masked and outlier BEDs so SNP-window boundaries neither lose covered sequence nor collapse to zero width.
- In `combine_remove_masks.py`, round every parsed interval outward before calling `merge_intervals`. Rounding after merging can create overlapping output rows from previously separated float intervals.
- Update `--inputs` help text: missing inputs are errors, not ignored files.
- Test a SNP window whose two float edges fall within the same integer and separated float intervals that overlap only after outward rounding.

Acceptance criteria:

- Missing or malformed mask inputs fail the job instead of producing a partial mask.
- Valid empty masks still combine successfully.
- Masked and outlier BED writers apply the same outward rounding rule because step 5 consumes the outlier BED directly.
- Combined output remains sorted and disjoint after rounding.

### 3. Write the correct chromosome into step-2 BED output

Files:

- `scripts/find_low_access_regions.py`
- `Snakefile`
- Relevant tests

Work:

- Add a required `--chrom` argument.
- Pass `{wildcards.chrom}` from the step-2 rule.
- Use that argument for BED column 1 instead of `ts_path.stem`, which is the replicate ID in the nested input layout.

Acceptance criteria:

- A step-2 BED generated from `<chrom>/<replicate>.trees` contains `<chrom>` in column 1.
- Standalone and pipeline invocation use the same explicit chromosome interface.

### 4. Preserve existing metadata during coordinate-preserving trimming

Files:

- `scripts/trim_regions_single.py`
- `tests/test_trim_regions_single.py`

Work:

- Merge `kept_intervals` and resolved mutation-map metadata into existing top-level metadata instead of replacing it wholesale.
- Add a test with unrelated input metadata and confirm it survives step 4.
- Continue allowing newly computed pipeline fields to replace stale fields of the same name.

Acceptance criteria:

- Existing provenance and application metadata survive trimming.
- `kept_intervals` and mutation-rate metadata reflect the trimmed output.

### 5. Remove the ineffective trimmed-tree validator

Files:

- `scripts/argtest_common.py`
- `scripts/trim_samples.py`
- `scripts/filter_min_samples.py`
- Related tests and log expectations

The current `validate_trimmed_ts` compatibility probing does not validate the pinned tskit tree sequence: its `check_index` fallback cannot obtain an index, returns early, and never reaches a useful `validate` call. It also produces a misleading warning in normal jobs.

Work:

- Delete `validate_trimmed_ts` and both call sites.
- Rely on construction of `tables.tree_sequence()` to reject malformed table collections.
- If a later invariant needs explicit checking, implement that invariant against the pinned API rather than restoring generic version probing.

Acceptance criteria:

- Normal step-5 and step-5b logs contain no “skipped check_index” warning.
- Malformed tables still fail while constructing a tree sequence.

## Phase 2: rename sample-name normalization everywhere

Files likely affected:

- `Snakefile`
- `config/snakemake.yaml`
- `scripts/argtest_common.py`
- `scripts/mutload_masks.py`
- `scripts/mutload_summary.py`
- `scripts/trim_samples.py`
- `scripts/export_vcf.py`
- Tests and fixtures referencing `suffix_to_strip`
- `README.md`, `CHANGELOG.md`, `MAKE_REALISTIC_EXAMPLE.md`, `NOTES.md`, `per_sample_trim_plan.md`, and other user-facing plans/notes where the option is current guidance

Work:

1. Rename the config, CLI flags, function parameters, local variables, and help text to `name_substring_to_remove` / `--name-substring-to-remove`.
2. Preserve the existing global `str.replace` behavior.
3. Add a Snakefile configuration check that reports the breaking rename when `suffix_to_strip` is present.
4. Add focused tests showing:
   - all occurrences are removed;
   - a middle occurrence is removed;
   - an empty configured substring leaves names unchanged;
   - matching remains consistent between mutation-load masks, sample trimming, and VCF naming.
5. Record the breaking rename in the changelog with a one-line migration example.

Acceptance criteria:

- `rg "suffix_to_strip|suffix-to-strip"` finds only changelog/migration text and any intentional rejection message.
- Step 3 and step 5 normalize IDs identically.
- Existing users receive an actionable error rather than silently using the default.

## Phase 3: audit the proposed individual/ploidy model without changing behavior

### Supported contract

For pipeline tree sequences:

- Every sample node belongs to one individual.
- Individual names after substring removal are unique.
- Every represented individual has the same number of sample nodes.
- One sample node per individual is haploid, two is diploid, and larger uniform values remain supported without special cases.
- Sample nodes are leaves in local trees.

The bundled datasets satisfy this proposed contract, but they are not sufficient evidence for the real amaranth/admix corpus. This phase is diagnostic only. It must not turn currently accepted inputs into hard failures or remove fallbacks.

### 1. Add a shared contract auditor

File:

- `scripts/argtest_common.py`

Work:

- Add a compact helper that reports individual-to-sample-node grouping and ploidy after name normalization.
- Report missing individual assignments, duplicate normalized names, mixed ploidy among represented individuals, and sample nodes used as edge parents.
- Ignore individuals with zero sample nodes when determining ploidy. They are unrepresented table rows, not zero-ploidy samples; `example_data/test_100trees.tsz` is a regression fixture for this distinction.
- Provide a warn-only audit command or mode that can scan tree sequences without changing pipeline output or exit status.
- Run it across every available bundled dataset and, before any enforcement, the real amaranth/admix ARG corpus.

Keep validation linear in the node/edge tables. Do not validate zero-length chromosomes, zero-sample ARGs, or other unsupported degenerate inputs repeatedly throughout the code.

### 2. Make enforcement a later, evidence-gated change

Files:

- `scripts/argtest_common.py`
- `scripts/export_vcf.py`
- `scripts/mutload_masks.py`
- `scripts/mutload_summary.py`
- `scripts/trim_samples.py`
- Relevant tests

Do not perform the following work in Phase 3. Schedule it only after the dead-code deletion and performance phases, and only if the audit finds no real-corpus violations or the user explicitly accepts the affected inputs becoming unsupported:

- Convert audit findings into clear entry-point errors.
- Replace `name_to_nodes_map` overwrite behavior with validated unique-name grouping.
- Keep pooling all sample-node loads for an individual.
- Remove the VCF fallback to one haploid column per node and always export by individual using validated uniform ploidy.
- Remove the internal-sample-node compatibility branch in mutation filtering after the leaf-sample contract is proven.

If the real corpus contains exceptions, retain the relevant small fallback or narrow the contract based on observed data rather than adding more generalized validation machinery.

Acceptance criteria:

- Haploid and diploid fixtures both pass.
- Individuals with zero sample nodes are ignored for uniform-ploidy calculation and do not trigger a warning or failure.
- Diploid loads are the sum of both sample nodes.
- Diploid VCF columns contain two phased alleles.
- The initial audit reports mixed-ploidy, missing-individual, duplicate-normalized-name, and internal-sample-node fixtures without failing.
- Hard failures and fallback removal land only after documented audit results for the real corpus.

## Phase 4: simplify mutation-map handling

Files:

- `scripts/argtest_common.py`
- `scripts/trim_regions_single.py`
- `scripts/validation_plots_from_ts.py`
- Step 6 staging logic in `Snakefile`
- Configuration and documentation

Work:

1. Preserve mutation-map resolution for raw inputs. Steps 2 and 3 read the original ARGs before any pipeline metadata is written, so exact sibling `*.mut_rate.p` lookup or an explicit scalar fallback remains required.
2. Simplify `infer_mu_path` to a documented set of exact candidate paths; remove suffix-derived fuzzy bases and the directory-wide `glob` fallback that can turn a missing map into a misleading ambiguity error.
3. Catch only `FileNotFoundError` where absence is genuinely supported. Allow corrupt pickle data, incompatible objects, and programming errors to fail.
4. Remove the production duck-typing shim that converts any pickled object with `position` and `rate` attributes into a `RateMap`; update fixtures to pickle the supported object.
5. Give `trim_regions_single.py` a `--mutation-rate` scalar fallback and replace its inline lookup/pickle block with `resolve_mu_rate(ts, args.ts, scalar_fallback=args.mutation_rate)`. Pass the same configured mutation-rate argument used by steps 2 and 3 from the step-4 Snakefile rule.
6. When step 4 resolves a scalar rate, convert it to a full-length `msprime.RateMap(position=[0, sequence_length], rate=[scalar])` before serialization. Embed the resolved map into the coordinate-preserving trimmed tree sequence while preserving unrelated metadata.
7. Confirm whether the bundled example intentionally uses its configured scalar rate: its maps at `example_data/sim_2chr_5rep/<chrom>.mut_rate.p` are not reachable from `trees/<chrom>/<rep>.trees` under the reduced exact candidate set. Either relocate the fixture maps or document scalar-rate behavior.
8. Require cleaned outputs after step 4 to use embedded metadata. Do not make this requirement apply retroactively to raw step-2/step-3 inputs.
9. Add `--ts FILE [FILE ...]` to `validation_plots_from_ts.py` and pass explicit original and cleaned files from Snakemake. Make `--ts` and `--ts-dir` mutually exclusive: `--ts` preserves caller order, while `--ts-dir` keeps natural sorting. Keep `--compare` directory-based or add a parallel `--compare-ts` and document the choice. The Snakefile must pass files in replicate order, with replicate IDs still derived by `extract_replicate_id`.
10. Cleaned validation uses embedded metadata; original files continue exact sibling-map resolution from their real paths. Delete step 6's temporary directories, symlink loops, inline Python map lookup, and `/tmp` dependency.

Do not add a step-0 rule merely to rewrite every raw ARG with embedded metadata. That would add substantial tree-sequence I/O and another DAG stage to eliminate a small, necessary raw-input lookup. Reconsider an explicit ingestion stage only if future input formats require broader normalization.

Acceptance criteria:

- Raw steps 2 and 3 resolve the intended exact sibling map or explicit scalar fallback.
- Scalar-only and map-backed configurations both embed a full rate map at step 4 for use by cleaned downstream outputs.
- Corrupt or unreadable maps produce hard failures; a missing map produces a `FileNotFoundError` naming the exact candidate paths.
- Step 6 no longer needs to infer, copy, or symlink a sibling mutation-rate pickle.
- A missing mutation map reports the exact paths tried rather than matching unrelated maps by prefix.
- Original validation resolves maps from original file paths; cleaned validation uses embedded metadata.

## Phase 5: delete superseded code

### Remove immediately after confirming references

- `scripts/trim_samples_chunked.py`
- `build_segments_with_drop_nodes`
- `trim_ts_by_intervals`
- `assert_sample_ids_preserved`
- `collapse_masked_and_low_access_windows`
- `build_shared_mask`
- `build_removal_segments` and the unused `remove_intervals` branch of `mutational_load`
- `and_ratemaps_binary`
- Tests that exist solely for these removed helpers

Before deletion, confirm with `rg` that no workflow, maintained script, or documentation imports them. Preserve tests of observable pipeline behavior; remove tests that only lock in dead implementations. This immediate list excludes `collapse_masked_intervals` and `ratemap_from_keep_intervals`, which remain live imports in `trim_regions.py` until that CLI and the `build_shared_mask` chain are removed in the ordered work below.

### Decide the status of coordinate-compacting `trim_regions.py`

The automated workflow uses `trim_regions_single.py`, which preserves chromosome coordinates. `trim_regions.py` is an auxiliary tool with different coordinate-compacting behavior. However, `trim_regions_single.py` currently imports `load_mask_intervals` and `complement_intervals` from it, so the module cannot be deleted directly.

Preferred action:

1. Move `load_mask_intervals` and `complement_intervals` into `argtest_common.py` (or a smaller shared interval module).
2. Update `trim_regions_single.py` and their tests to import the relocated helpers.
3. Verify coordinate-preserving step 4 before deleting `trim_regions.py`.
4. Delete `trim_regions.py` unless there is a known external coordinate-compaction use case. If retained, move it under a clearly labeled auxiliary location, rename it to communicate coordinate compaction, and keep it out of the suggested workflow.
5. Once the coordinate-compacting CLI and dead `build_shared_mask` chain are gone, delete their now-unreferenced helpers `collapse_masked_intervals` and `ratemap_from_keep_intervals` as well.
6. Delete `tests/test_trim_regions_split.py::test_trim_regions_applies_bed_mask` with the coordinate-compacting CLI. Retarget only tests of BED clipping/merging and interval complements to the relocated helpers.
7. Remove the obsolete `README.md` claim that `trim_regions.py` checks cross-replicate sequence lengths. Do not relocate that auxiliary-only check into the DAG: doing so would add new defensive complexity for an unobserved failure mode. Existing checks on live supported paths remain in scope.

Acceptance criteria:

- Shared code contains no functions with test-only callers.
- There is one supported sample-trimming implementation.
- The main workflow documents only coordinate-preserving trimming.
- `trim_regions_single.py` imports successfully after `trim_regions.py` is removed.
- Existing tests for BED clipping/merging and complement calculation follow the relocated helpers and continue to pass.

## Phase 6: reduce merge memory and redundant tree-sequence work

### 1. Batch chromosome concatenation

Files:

- `scripts/merge_treefiles_by_replicate.py`
- `tests/test_merge_treefiles.py`
- `measure_merge_mem.py` or a reproducible profiling script

Work:

- Confirm that the pinned tskit fork supports `first_ts.concatenate(*remaining)` with output identical to incremental concatenation.
- Load the chromosome tree sequences, collect their metadata inputs, and concatenate them in one batch.
- Preserve natural chromosome ordering and the existing `chrom_offsets`, rate-map, and kept-interval metadata behavior.
- Benchmark wall time and peak RSS on the pinned environment. The independent review measured about 34% lower peak RSS and 15% lower wall time on a representative 16-chromosome input; treat that as a target to reproduce, not an assumed guarantee.
- Adopt batched concatenation only if peak RSS on a representative real-corpus replicate is no worse than the incremental loop. If holding all inputs regresses memory, keep the incremental implementation and record the result; memory reduction, not wall-clock improvement, is the reason for this change.

### 2. Release the merged tree sequence before metadata reconstruction

- Release input tree sequences immediately after batched concatenation returns and before `merged.dump_tables()`, because `dump_tables()` is itself a whole-genome allocation. Drop `first_ts`, the remaining-input list, and loop-local references together.
- Structure metadata collection so it stores only the small extracted rate-map/interval/offset data needed after concatenation, not references to the input tree sequences.
- Read `existing = merged.metadata` before dumping tables. After `tables = merged.dump_tables()`, delete `merged` immediately so only the table collection remains live when `tables.tree_sequence()` constructs the final object.
- Measure peak RSS around this stage; retain the change only if it reduces or does not worsen memory under the pinned runtime.

### 3. Remove the redundant full-sequence copy in `mutational_load`

After deleting the unused interval-removal branch, iterate directly over the input `ts`. The current normal path calls `keep_intervals([(0, sequence_length)], simplify=False)`, which copies and reindexes the complete tree sequence without changing it.

Do not replace the mutation traversal with a genotype/`variants()` sweep in this phase. That alternative counts derived-state presence rather than every mutation below a sample and differs at recurrent/back-mutation sites.

### 4. Avoid duplicate coalescence input loading when memory permits

`coalescence_ne_plots_from_ts.py` currently loads post-burnin replicates once to derive quantile/log-spaced time windows and again to calculate statistics.

- Refactor grid construction to accept already-loaded tree sequences and reuse them for statistics when the resulting aggregate memory is acceptable.
- Benchmark peak RSS as well as I/O time on whole-genome merged ARGs. Holding every decompressed replicate concurrently may be worse than a second sequential read.
- If reuse materially increases peak memory, retain the two-pass behavior and document it as an intentional memory/I/O tradeoff rather than forcing an optimization that risks OOM.

### 5. Scan each HapMap once per invocation

`hapmap_low_rec_mask.py` currently scans the full HapMap once to collect chromosome names and again to load the selected chromosome; the per-chromosome Snakemake rule multiplies this across chromosomes.

- Parse the HapMap once into the information needed for name resolution and selected rows.
- Preserve current chromosome-alias resolution and output semantics.
- Consider a single all-chromosome job only as a separate DAG design change; it is not required to eliminate the duplicate scan within each invocation.

Acceptance criteria:

- Merged table contents and metadata match the current implementation.
- Peak merge RSS is measured and does not regress; the expected outcome is a substantial reduction.
- Fixed-seed mutation loads are bit-for-bit unchanged.
- `mutational_load` no longer calls `keep_intervals` in its supported path.
- Each HapMap rule invocation performs one input scan.
- Coalescence input reuse is adopted only if its peak-memory profile is acceptable.

## Phase 7: consolidate duplicated scientific helpers

Files:

- `scripts/argtest_common.py`, or a new small focused module such as `scripts/mutload_common.py`
- `scripts/mutload_masks.py`
- `scripts/mutload_summary.py`
- `scripts/validation_plots_from_ts.py`

Work:

- Consolidate base-pair window construction.
- Consolidate SNP-window construction.
- Consolidate simulation-based expected-load calculation.
- Consolidate the repeated “arange then clamp the final edge” window construction used by accessibility and coalescence code where the contracts are identical.
- Keep the seeded `msprime.sim_mutations` call.
- Keep the phrase “simulation-based expected load” in code help and documentation.
- Document that one seeded simulation draw estimates the expected load and therefore introduces simulation variance.

Prefer a focused mutation-load module over adding more unrelated responsibilities to the already broad `argtest_common.py`.

Acceptance criteria:

- One implementation exists for each window builder and expected-load calculation.
- The same seed, ARG, mutation map, and windows produce identical results before and after refactoring.
- Help text consistently says “simulation-based expected load.”

### Remove unreachable and dead compatibility branches during consolidation

- Remove the expected-name ordering check in `mutload_masks.py`; both orderings are produced from the same `names` object by the same function.
- Remove `pipeline_summary` fallbacks for absent `retained_by_rep`; `collect_retention` always populates it.
- Replace `getattr(args, "mutation_rate", None)` and `getattr(args, "log", None)` with direct argparse attributes.
- Remove the unused `seq_len` parameter from `_intersect_keep_with_drops`.
- Import the shared `load_ts` in the two HTML tree-view scripts rather than maintaining duplicate implementations.
- Reuse the shared interval merger in `score_realistic_example.py` unless keeping the example generator deliberately standalone is documented.

## Phase 8: narrow dependency compatibility and plotting edge cases

### Dependency API

- Add one cheap fail-fast check that the installed tskit provides the pinned fork capability required for partial-missing-data pair-coalescence normalization. Prefer an explicit version/feature marker; do not run a substantial synthetic calculation at every import.
- Test against the pinned environment rather than maintaining runtime compatibility branches for unsupported versions.
- After Phase 5 removes dead helpers, resolve the remaining production-assert policy. Convert `mutational_load`'s caller-supplied window-grid assertion to an explicit `ValueError`; delete assertions that only restate invariants guaranteed by the pinned dependency or live solely in removed rate-map helpers. Do not add replacements for deleted defensive checks.

### Empty and degenerate data

- Do not add special behavior for one-base chromosomes.
- Retain graceful handling for missing optional plot files.
- Decide whether mutation-free but otherwise valid ARGs are common enough to deserve empty plots. If not, fail once with a clear message and remove downstream all-zero/all-NaN plotting guards.

This decision should be based on real pipeline inputs, not hypothetical completeness. Mutation-free windows are normal and must remain supported; the question is only whether an entire ARG with no mutations is a supported reporting input.

### Plot lifecycle and duplication

- Replace `plt.clf()` with `plt.close(fig)` after each saved figure in validation and coalescence plotting scripts so pyplot does not retain cleared figure objects.
- Keep `safe_log_yscale`; it prevents a real Matplotlib failure when an SFS has no positive values.
- Table-drive the repeated diversity, segregating-sites, and Tajima's D scatter/skyline/trace blocks in `validation_plots_from_ts.py` using a small statistic specification table.
- Preserve each statistic's current NaN and accessibility behavior explicitly in the table or shared plotting helper; add output-level tests before deleting the duplicated blocks.

## Phase 9: documentation corrections

Files:

- `README.md`
- `CHANGELOG.md`
- CLI help strings
- Any current plan/notes referenced by README

Work:

- Correct the default output path for `trim_samples.py`, or change code to match the documented default; choose one behavior.
- Name `trim_regions_single.py` as the pipeline's coordinate-preserving trimming command.
- Clearly separate Snakemake stdout/stderr logs from per-script success-summary logs.
- Add a “Supported input assumptions” section containing the individual/ploidy model, leaf-sample requirement, positive chromosome lengths, nonempty ARGs, valid BED requirement, and pinned dependency policy.
- Explain the breaking `name_substring_to_remove` rename.
- State that expected load is estimated with one reproducible, seeded mutation simulation.
- Make `pipeline_summary.py` and `validation_plots_from_ts.py` use and describe the same accessibility precedence: `kept_intervals`, then embedded mutation-rate accessibility, then the explicitly documented fallback.
- Decide whether `mutload_summary.py --out` should honor the supplied path or be renamed `--out-name`; do not silently discard a caller-supplied directory.
- Move completed root-level design plans into `reports/` or remove them after confirming their durable decisions are represented in README/CHANGELOG.
- Fold the short sample-ID matching guidance from `NOTES.md` into the supported-input section if that file otherwise has no continuing purpose.

Recommended output-path decision: keep the implemented sibling directory (`<input-parent>/trimmed/<stem>_trimmed.tsz`) and update README/help text. It keeps standalone output near its input and avoids an implicit dependency on the repository root.

## Verification strategy

### Step 0: provision the declared environment

Before implementation or per-phase checks, rebuild or synchronize the `argtest` Conda environment from `environment.yml`. Confirm that the activated environment provides the pinned tskit fork, `pytest`, and `snakemake`. Stop if activation or dependency verification fails; do not silently run tests in base Python.

### Fast checks after each phase

- `python -m compileall` for scripts and tests.
- Targeted `pytest` modules for changed functionality.
- `snakemake -n -p --configfile config/snakemake.yaml` against the example data.
- `rg` checks for retired names and dead imports.
- `git diff --check`.

### End-to-end checks before completion

1. Run the full test suite.
2. Run the example pipeline from a completely absent output directory.
3. Run at least one haploid and one diploid fixture.
4. Compare pre/post outputs for a fixed random seed:
   - combined masks;
   - trimmed tree sequence statistics;
   - simulation-based expected loads;
   - outlier calls;
   - VCF sample columns and genotypes;
   - merged chromosome offsets and rate metadata.
5. Confirm step-2 BED chromosome labels and preservation of unrelated step-4 metadata.
6. Confirm pipeline summary and validation report the same accessible span for the same final tree sequence.
7. Record wall time and peak RSS for merge, step 3, step 6, and coalescence plotting to ensure simplification does not regress performance.

## Suggested commit sequence

Keep changes reviewable and behaviorally isolated:

1. Fix output parents, strict BED parsing, the step-2 chromosome label, and step-4 metadata preservation.
2. Batch merge concatenation, release the pre-metadata merged object, and record pinned-environment memory measurements.
3. Delete the ineffective trimmed-tree validator and its call sites.
4. Delete unused helpers, the unused mutation-removal branch, `trim_samples_chunked.py`, and the coordinate-compacting CLI in the dependency-safe order from Phase 5.
5. Rename `suffix_to_strip` to `name_substring_to_remove` with docs and migration error. Deleting dead code first avoids mechanically renaming a module and fixtures that disappear immediately afterward.
6. Add and run the warn-only individual/ploidy contract audit; do not remove fallbacks yet.
7. Make mutation-map propagation explicit, add direct step-6 file arguments, and remove broad/fuzzy fallbacks.
8. Remove the redundant whole-tree copy and consolidate mutation-load helpers without changing seeded results.
9. Review real-corpus audit results; only then enforce the individual/ploidy contract, simplify leaf-sample mutation filtering, and remove the VCF fallback if supported by evidence.
10. Consolidate plotting blocks, close figures correctly, and unify accessibility accounting.
11. Narrow pinned-version compatibility and supported empty-data behavior.
12. Finish documentation and run the full end-to-end verification matrix.

Each commit should leave tests passing and avoid mixing mechanical renames with scientific behavior changes, except where the rename's migration error and documentation must land atomically.
