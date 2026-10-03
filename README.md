# ARG Tree Sequence Utilities and Validation Plotting

Snakemake pipeline for post-processing, QC, and basic visualization of ARG tree sequences (`.ts`, `.trees`, `.tsz`). Written with the aid of [Codex](https://openai.com/codex/) and [Claude](https://claude.ai/).

If you use, please cite:

Ross-Ibarra, J. 2026. ARGtest: tools for QC and validation of ancestral recombination graphs. [doi: 10.5281/zenodo.19698118](https://doi.org/10.5281/zenodo.19698118)

![ARGtest overview: inferred ARGs are screened for low recombination, poor accessibility and aberrant mutation load; flagged regions are trimmed from all samples and outlier individuals are pruned only in the windows where they fail; outputs are cleaned ARGs, validation plots and an HTML summary report.](docs/graphical_abstract.svg)

## Documentation

This README covers installing and running the pipeline. The details live in [docs/](docs/):

- [The workflow](docs/workflow.md): what each step does and why, and caveats for statistics after sample pruning
- [Configuration](docs/configuration.md): every config key, where the mutation rate comes from, and which file names must match
- [Outputs](docs/outputs.md): the output layout, merged-file naming, finding a tree by position, and VCF export
- [Scripts](docs/scripts.md): every pipeline and auxiliary script, plus input formats and defaults
- [Realistic example dataset](docs/realistic_example.md) and [coalescence / Ne plots](docs/coalescence_ne_plots.md)
- [Changelog](CHANGELOG.md)

## Install

```bash
conda env create -f environment.yml
conda activate argtest
```

`tskit` is not installed from conda-forge: `environment.yml` pins it to
[nspope/tskit commit `73d8cd9`](https://github.com/nspope/tskit/commit/73d8cd922482475020ae01180cae95bf5abbf067)
and installs it with pip, so the build needs pip, git, and access to GitHub. The
pinned build is what provides the native partial-missing-data pair-coalescence
normalization (see [Auxiliary scripts](docs/scripts.md#auxiliary-scripts)) and it requires
Python ≥ 3.11 and NumPy ≥ 2. When upgrading from an environment built before
v1.9, recreate it rather than updating in place, so the conda-forge `tskit` is
replaced:

```bash
conda env remove -n argtest
conda env create -f environment.yml
```

**Which version to use:** work from the most recent tagged release rather than an
arbitrary commit on `main`. To check out the latest tag:

```bash
git fetch --tags
git checkout "$(git tag -l 'v*' | sort -V | tail -1)"
```

See [CHANGELOG.md](CHANGELOG.md) for a per-version breakdown of changes.

## Quick start

New here? After [installing](#install), run the whole pipeline end-to-end on the bundled example dataset:

```bash
# 1. Generate the example tree sequences. These are not committed; they
#    regenerate deterministically from the seed pinned in ground_truth.json
#    (3 chromosomes × 8 replicates × 16 diploids × 10 Mb; a few minutes).
python scripts/make_realistic_example.py --out-dir argtest-realistic-example \
    --n-chrom 3 --n-reps 8 --n-samples 16 --seq-length 10000000

# 2. Run the pipeline. The committed config/snakemake.yaml already points at
#    the example dataset above.
snakemake --cores 16 --configfile config/snakemake.yaml
```

Results land in `argtest-realistic-example-out/`. The example is a deliberately-flawed dataset (contaminated individuals, per-window sample pruning, an accessibility mask) with a ground-truth scorecard at [scoring_report.md](argtest-realistic-example-out/scoring_report.md) grading the pipeline's masks against the injected flaws; see [docs/realistic_example.md](docs/realistic_example.md) for the generator's CLI and ground-truth schema. For dry-run preview, run flags, and cluster execution, see [Run the pipeline](#run-the-pipeline) below.

## Set up your data

The pipeline expects one subdirectory per chromosome, holding one treefile per replicate:

```text
<root>/
  chr1/
    1.tsz
    2.tsz
    ...
  chr2/
    1.tsz
    2.tsz
    ...
```

The Snakemake workflow discovers treefiles with suffixes `.ts`, `.trees`, and `.tsz`. Replicate IDs are taken from the filename stem, so `chr1/1.tsz` is replicate `1` for chromosome `chr1`.

If the treefiles live in a subdirectory of each chromosome directory rather than directly inside it — for example SINGER output where each `chrN/trees/` holds the replicates — set `tree_subdir` to that subdirectory name and discovery looks there instead:

```text
<root>/
  chr1/
    trees/
      chr1.1.tsz
      chr1.2.tsz
      ...
```

The chromosome name still comes from the chromosome directory (`chr1`), and a leading `chrN.` prefix on the filename is stripped to get the replicate ID, so `chr1/trees/chr1.2.tsz` is replicate `2`.

Then copy [config/snakemake.yaml](config/snakemake.yaml) and edit it for your project. Every option has an inline comment. The keys you must set are:

- `root_dir`: the chromosome-subdirectory root shown above
- `hapmap`: one HapMap recombination map covering all chromosomes
- `fai`: FASTA index, for chromosome lengths
- `rec_fraction`: fraction of recombination-map intervals to mask as low recombination (step 1)
- `low_access_window_size` and `low_access_cutoff_bp`: window size and accessible-bp cutoff (step 2)
- exactly one of `mutload_window_size` or `mutload_snp_window` (step 3)

See [docs/configuration.md](docs/configuration.md) for the optional keys, where the mutation rate comes from, and how chromosome names must match between your directories, the HapMap and the `.fai`.

## Run the pipeline

From the repo root:

```bash
module load conda
conda activate argtest
```

With your config ready (see [Set up your data](#set-up-your-data)), preview the planned jobs with a dry run, then run for real:

```bash
snakemake -n -p --configfile config/snakemake.yaml
snakemake --cores 16 --rerun-incomplete --keep-going --configfile config/snakemake.yaml
```

`-n -p` prints the planned jobs and their commands without executing; `--rerun-incomplete` re-runs any jobs a previous interrupted run left half-finished; `--keep-going` lets independent jobs continue when one fails.

This is fine for the small example dataset. **For real datasets, prefer the [SLURM route](#running-on-a-slurm-cluster) below.** The `merge_replicates` step loads all chromosomes of a replicate into memory at once and can need ~50–128 GB each; with plain `--cores N`, Snakemake may launch several merges concurrently and OOM the machine. If you must run locally, cap memory with `--resources mem_mb=<node_RAM_in_mb>` so the merges serialize to fit.

### Running on a SLURM cluster

Add `--profile profiles/slurm` to submit each job to SLURM instead of running locally:

```bash
snakemake --profile profiles/slurm --configfile config/snakemake.yaml
```

This is the recommended way to run real datasets: besides parallelism, it fans the ~per-(chromosome, replicate) jobs out across the cluster, sending light steps to one partition and the memory-heavy `merge_replicates` / `step6_validation_plots` steps to a big-mem partition — so each merge gets its own right-sized node and avoids the local-run OOM noted above. **All cluster knobs live in your configfile** — `slurm_account`, `slurm_partition`, and the per-rule `resources:` block (mem/time/threads/partition). The profile file itself ([profiles/slurm/config.yaml](profiles/slurm/config.yaml)) is set-and-forget; you do not edit it.

Notes:
- **Partition names are cluster-specific.** The shipped config uses `"low"` and `"high"`, which are particular to our cluster; change `slurm_partition` (and any per-rule `partition` in the `resources:` block) in your configfile to match the partitions on your own cluster.
- **`out_dir` must be on a shared filesystem** (not node-local `/tmp`) — each job runs on a different node, so a `/tmp` output path silently loses results.
- The snakemake process here is just a controller (it submits jobs and polls), but it runs for the whole pipeline — launch it inside `tmux`/`screen` or a small `srun` so it survives disconnects.
- Per-job SLURM logs are written to `logs/slurm/`.
- Plain `snakemake --cores N …` (above) still runs everything locally; SLURM-only settings such as account, partition, memory, and walltime are ignored, while per-rule `threads` still affects local scheduling.

### Sandboxed environments

Some HPC or container setups mount `~/.cache` read-only. There, prefix the Snakemake command (dry-run or real run) with cache and temp-dir redirects to `/tmp`:

```bash
XDG_CACHE_HOME=/tmp/argtest-xdg-cache TMPDIR=/tmp/argtest-tmp \
    snakemake --cores 16 --configfile config/snakemake.yaml
```

On a normal machine where `~/.cache` is writable this is not needed.

## Outputs

By default, Snakemake writes outputs beneath `out_dir` with subdirectories for each stage:

```text
<out_dir>/
  step1_low_rec/
  step2_low_access/
  step3_mutload/
  step4_masks/
  step4_trimmed_regions/
  step5_trimmed_samples/
  step5b_min_samples/    # low-sample intervals dropped, if min_samples is set
  combined/
  step6_validation/    # step 6 validation plots + plot-data TSVs (original and cleaned), if configured
    <chrom>/{original,cleaned}/   # per-chromosome plots, windows.tsv, samples.tsv, sfs.tsv
    genomewide/{original,cleaned}/  # step 6b pooled expected-vs-observed plots + genomewide-windows.tsv
  vcf/                 # per-(chromosome, replicate) VCFs, if emit_vcf: true
  pipeline_summary.html
  logs/
```

Intermediate filenames include both chromosome and replicate information so they stay unique across the full workflow.

The final merged tree sequences are `combined/<base_name>.combined.<replicate>.<suffix>`, and `pipeline_summary.html` collects retention statistics and all plots in one page. See [docs/outputs.md](docs/outputs.md) for naming details, how to find the tree at a chromosome position in a merged file, and VCF export.

## Repository notes

- Generated `logs/` and `results/` are git-ignored.
- `.DS_Store` is git-ignored.
- Per-release changes are tracked in [CHANGELOG.md](CHANGELOG.md), keyed to the
  `v1.x` git tags.
- Developer notes, design plans and code reviews live in [docs/dev/](docs/dev/).
- The overview figure is generated by [docs/make_graphical_abstract.py](docs/make_graphical_abstract.py).

## Acknowledgements

None of this would be possible without the patient help and advice of Nate Pope. Any errors, bad code, or poor interpretations, however, are my responsibility alone. This repo also uses code from Nate Pope's [singer-snakemake](https://github.com/nspope/singer-snakemake).
