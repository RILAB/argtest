# Configuration and input requirements

← [Back to README](../README.md)

For the input directory layout, see [Set up your data](../README.md#set-up-your-data).

## Required keys

The Snakemake config is in [config/snakemake.yaml](../config/snakemake.yaml). Edit it for your project and supply it with `--configfile`. The file has an inline comment for every option.

- `root_dir`: path to the chromosome-subdirectory root
- `hapmap`: single HapMap recombination map covering **all** chromosomes (one combined file, not one per chromosome — rows are grouped by the `Chromosome` column), used for step 1
- `fai`: FASTA index used for chromosome lengths
- `rec_fraction`: fraction of recombination-rate **intervals** (ranked by `Rate(cM/Mb)`, lowest first) to include in the low-recombination mask; e.g. `0.1` masks the bottom 10 % of intervals, while `0.0` writes empty low-recombination masks. Note this is a fraction of the *number of intervals between map markers*, not of base pairs — because the lowest-recombination intervals tend to be the longest, masking the bottom 5 % of intervals can remove well over 5 % of the genome by bp
- `low_access_window_size`: window size in bp for step 2
- `low_access_cutoff_bp`: minimum accessible bp per window for step 2
- exactly one of `mutload_window_size` or `mutload_snp_window` for step 3

## Optional keys

All optional keys have sensible defaults.

- `tree_pattern`: glob for treefiles within each chromosome directory (default: `"*"`), for example `"*.trees"` or `"*.tsz"`
- `tree_subdir`: optional subdirectory within each chromosome directory that holds the treefiles (default: unset → files live directly in the chromosome dir); e.g. `"trees"` for SINGER-style `chrN/trees/` layouts
- `mutload_cutoff`: outlier cutoff fraction for step 3 (default: `0.5`)
- `mutation_rate`: single scalar mutation rate (per bp per generation), the shared fallback for **both** step 3 and step 6 when no embedded or sibling ratemap is available
- `mutload_random_seed`: base seed for the per-replicate mutation simulation in step 3 (default: `1`)
- `mutload_fraction`: fraction threshold for writing mutation-masked BED rows in step 3
- `name_substring_to_remove`: substring removed from sample IDs before matching in step 3 and step 5 (default: `"_anchorwave"`). Removal uses global string replacement, so every occurrence is removed, not only a terminal suffix. The former `suffix_to_strip` key is no longer accepted; rename it explicitly in existing configs.
- `trim_individuals`: extra individual IDs removed genome-wide in step 5, **in addition** to the step-3 mutload outliers (e.g. introgressed samples). Comma-separated string (`"id1,id2"`) or a YAML list; normalized with `name_substring_to_remove` like step 3. Default: unset (trim only mutload outliers)
- `trim_remove_bed`: extra BED file(s) of per-individual intervals removed in step 5, in addition to the mutload outliers. Column 4 holds comma-separated sample IDs (or the filename stem if absent). Single path or a YAML list of paths; applied identically to every (chrom, rep). Default: unset
- `min_samples`: minimum number of non-isolated **retained sample nodes** (haploids, not individuals) a local tree must have; intervals below it are dropped in an optional **step 5b** that runs after step 5 and before the merge. Dropping uses `delete_intervals`, so sequence coordinates are **preserved** (the removed spans become empty gaps — there is no coordinate compaction) and the `kept_intervals` metadata is intersected with the surviving spans so downstream accessibility is not overestimated. Each `(chrom, rep)` also gets a diagnostic BED of the dropped intervals (`chrom start end retained_samples min_samples`) under `<out_dir>/step5b_min_samples/`. When set, the merge, VCF export, and validation steps automatically consume the step-5b filtered tree sequences instead of the raw step-5 output. Default: unset/null (skip step 5b entirely, so existing configs are unaffected). Sample pruning in step 5 is what creates the low-sample intervals this filter targets
- `allow_missing_replicates`: set to `true` to concatenate partial replicate sets (default: `false`)
- `burnin`: number of leading discovered replicates to discard before concatenation (default: `0`); must be smaller than the number of replicates discovered after applying `tree_pattern`
- `base_name`: prefix used for merged outputs (default: name of `root_dir`)
- `merged_out_suffix`: force a specific output suffix for merged files (`.ts`, `.trees`, `.tsz`); default is to inherit the suffix of the first input
- `out_dir`: output root for Snakemake products (default: `snakemake_out`; tilde is expanded)
- `run_validation`: master switch for step 6 (default: `true`); set `false` to skip the validation plots while keeping `mutation_rate` set for step 3. Step 6 also auto-skips when no rate source is available (neither `mutation_rate` nor `validation_sim_branch`)
- `validation_first_chrom_only`: run step 6 only on the first chromosome (default: `true`)
- `validation_window_size`: window size in bp for step-6 diversity, Tajima's D, and segregating-sites validation plots (default: `100000`); larger windows run faster and use less memory on large ARGs but give coarser QC curves
- `validation_sim_branch`: simulate site mutations on each ARG replicate with msprime for a posterior-predictive check instead of scaling branch statistics (default: `false`); can run without a scalar `mutation_rate` when every validated tree sequence has an embedded/sibling ratemap
- `emit_vcf`: if `true`, export one `.vcf.gz` per (chromosome, replicate) from the filtered per-chromosome tree sequences into `<out_dir>/vcf/` (default: `false`); see [VCF export](outputs.md#vcf-export)
- `vcf_reps`: restrict VCF output to specific replicate IDs (a subset of the post-`burnin` replicates); leave unset/null to emit every post-`burnin` replicate

## Where the mutation rate comes from

Steps 2 and 3 operate on raw ARGs, before step 4 embeds pipeline metadata, so they resolve a rate from exact sibling `*.mut_rate.p` candidates and then the scalar `mutation_rate` fallback. Step 4 embeds the resolved map (including a scalar rate represented across the sequence), and cleaned downstream ARGs use that embedded metadata. For reporting accessibility, `kept_intervals` takes precedence, followed by positive-rate intervals in the embedded mutation map, followed by the documented whole-sequence fallback when neither is available. A flat scalar gives step 3 a uniform-rate expectation rather than spatial correction; prefer a real map when local mutation-rate variation matters.

Step 3's **simulation-based expected load** is estimated from one mutation simulation on the input ARG. `mutload_random_seed` makes that draw reproducible, but a single draw still has simulation variance; it is not an analytic expectation or an average over many simulations.

## Supported input assumptions

The streamlined pipeline is intended for nonempty, positive-length tree sequences with valid half-open BED intervals and sample nodes assigned to named individuals. Per-individual mutation-load analysis pools all sample nodes belonging to an individual, so diploids may have two sample nodes. The proposed stricter contract is uniform ploidy among **represented** individuals and leaf sample nodes; individual-table rows with no sample nodes are ignored when assessing ploidy.

Those ploidy and leaf-sample conditions are currently audit targets, not a claim that every entry point already enforces them. Compatibility fallbacks should be removed only after the warn-only audit passes on the real amaranth/admix corpus (or affected inputs are explicitly declared unsupported). The project uses the pinned dependency versions in `environment.yml`; arbitrary tskit versions are not a supported compatibility target.

## File naming and what must match

The pipeline derives the chromosome label, the replicate ID, and (optionally) the mutation map from your directory layout and filenames. Two of these must line up with the *contents* of your input files; getting them wrong is a common setup failure.

**Startup chromosome validation.** During Snakefile evaluation, before the DAG is built or any local/SLURM job is submitted, the pipeline checks every chromosome label discovered from the ARG directory layout against both the HapMap `Chromosome` column and the first column of the `.fai`. If a label cannot be resolved, the workflow stops once with a `Chromosome naming mismatch detected before job execution` error that lists all unmatched pipeline chromosomes and the available reference names. For example, ARG directories named `chrom1` through `chrom10` do not match HapMap/FAI names `1` through `10` under the aliases described below, so the complete mismatch is reported at startup instead of failing later in ten separate step-1 jobs.

The workflow has no input VCF whose chromosome names need an independent check. When `emit_vcf: true`, each output VCF's `CHROM` value and contig header are generated from the already-validated ARG chromosome-directory label, so the ARG-derived output VCF label agrees by construction. The startup check therefore covers **ARG directory labels ↔ HapMap ↔ FAI**; it does not inspect chromosome metadata embedded inside a tree sequence.

- **Chromosome label** — the name of each chromosome subdirectory directly under `root_dir` (e.g. `chr1` in `<root>/chr1/...`; see the layout diagrams in [Set up your data](../README.md#set-up-your-data)). This one string is written as the BED chromosome column **and** is the key looked up in your hapmap and `.fai`. Keep it short and chromosome-like (`1`, `chr1`); **do not embed run descriptions in the directory name** — a directory called `chr.10.combined.snp.te.sorted` becomes the label verbatim and will not match a normal hapmap. (If your treefiles sit one level deeper, e.g. `chrN/trees/`, set `tree_subdir` rather than pushing `root_dir` down or lengthening the chromosome directory name.)

- **Replicate ID** — the treefile stem with a leading `<chromosome-label>.` prefix stripped, so `chr1/trees/chr1.2.tsz` is replicate `2`. One treefile per replicate per chromosome.

- **Hapmap `Chromosome` column** — the hapmap is one combined file given by the `hapmap:` path; **its filename is irrelevant**. For each chromosome the pipeline looks up the chromosome label in the `Chromosome` column, trying, in order: the label as-is, the substring after its first `.`, then `chr<that>` and `chr_<that>`. So a label `amaranth.1` matches a column value of `1`, `chr1`, or `chr_1`; a bare label `chr1` matches only `chr1` (there is no dot to strip, so `1` alone would *not* match). The **same matching is applied to the first column of the `.fai`.** A label like `chr.10.combined.snp.te.sorted` only ever reduces to `10.combined.snp.te.sorted`, which is why such names fail to match.

- **`*.mut_rate.p` files** (optional — see [Where the mutation rate comes from](#where-the-mutation-rate-comes-from)) — unlike the hapmap, these **are** discovered from the treefile path, independent of `root_dir`. For a treefile `<directory>/<stem>.<suffix>`, the exact candidates are `<stem>.mut_rate.p` and `<directory-name>.mut_rate.p`, searched first beside the treefile and then one directory above it. The first exact match wins. For example, `.../chr1/12.tsz` can use `.../chr1/12.mut_rate.p` (replicate-specific) or `.../chr1/chr1.mut_rate.p` (chromosome-wide). If none exists and `mutation_rate` is unset, the error lists every path tried. Broad glob and trailing-replicate heuristics are intentionally unsupported.
