# Realistic synthetic example dataset

`scripts/make_realistic_example.py` generates a synthetic but realistic
tree-sequence dataset for end-to-end testing of the pipeline. It
simulates ARGs under a known demography, then injects three independent
flaws so the pipeline has something concrete to detect and clean up.
A `ground_truth.json` records exactly what was injected, enabling
precision/recall scoring of the pipeline's masks.

## Quick start

```bash
python scripts/make_realistic_example.py --out-dir /tmp/realistic-example
```

Defaults: 2 chromosomes × 5 replicates × 8 diploid individuals × 5 Mb,
with the two-bottleneck demography (35 ka and 9 ka), 2 contaminated
individuals (severity multipliers 2× and 5×), 50% of windows pruned,
33% of the genome masked as inaccessible, and a centromere-like low-rec
valley in the hapmap.

## The three injected flaws

| Flaw | Mechanism | Tests |
|---|---|---|
| Extraneous mutations on K individuals | Extra leaf-edge mutations sized by `(severity − 1) × mean per-individual load`; placed at random positions, assigned to a random sample node of the contaminated individual. With `--contam-hotspot-frac > 0`, a fraction of each individual's extras concentrates in a per-(chrom, individual) hotspot of `--contam-hotspot-size` bp (paralog-mapping mimic); rest stay uniform | `mutload_masks.py` outlier detection — both whole-genome (uniform mode) and per-window (hotspot mode) granularity |
| Per-window sample pruning | A fraction of windows is chosen; in each, a random subset of samples has its ancestry excised (edges removed; no `simplify`, so sample IDs stay stable). The same prune pattern is shared across MCMC-like replicates of one chromosome | `trim_samples.py` semantics, the partial-tree Ne correction in `compute_pair_coal`, and pipeline tolerance to isolated-sample regions |
| Genome-wide accessibility mask | `chrN.mut_rate.p` is an `msprime.RateMap` with rate = 0 over a fraction of the genome (default ~3 contiguous intervals totaling 33%). `sim_mutations` honors this — no mutations land in masked regions — so downstream accessibility detection has a clean signal | `find_low_access_regions.py`, the sim-based mutload expectation in `mutload_masks.py`, and `coalescence_ne_plots_from_ts.py`'s coordinate handling |

The severity baseline is the actual mean per-individual mutational load
computed via `argtest_common.mutational_load` (not the naive
`num_mutations / num_individuals`, which underestimates because most
mutations are carried by multiple individuals via shared ancestry). So
`--contam-severity 5` aims for the contaminated individual's load to
land roughly 5× the non-contaminated mean. Poisson noise and pruning
chopping injected leaf mutations introduce some variance — empirically
the actual ratio runs ~3–5× for severity 5 in small (~500 kb) tests
and tightens with longer sequence.

## MCMC-like replication

Real `singer` output is N tree-sequence files per chromosome,
representing posterior draws from MCMC. They share the same observed
mutation data but have different inferred ARG topologies. We can't
reproduce this exactly without running an inference step, so we
approximate:

| Property | Shared across replicates? | Why |
|---|---|---|
| Underlying coalescent realization | **No** (different seed per rep) | Mimics inference variance |
| Contaminated individuals + severity | **Yes** | Contamination is a property of the person, not the inference — same individual is dirty in every MCMC sample |
| Per-window prune mask | **Yes** | trim_samples decisions come from observed data, not inferred topology |
| `mut_rate.p` accessibility mask | **Yes** | Upstream of MCMC |

This means the pipeline sees: "across replicates, individual X looks
contaminated, window Y is masked, but the trees are slightly different"
— matching how real MCMC output behaves to the pipeline's averaging
code.

## Output layout

```
<out-dir>/
  chr1/
    rep_000.trees, rep_001.trees, ...   # MCMC-like replicates
  chr2/
    ...
  chr1.mut_rate.p      # msprime.RateMap, pickled; rate=0 over masked intervals
  chr2.mut_rate.p
  chr1.hapmap.tsv      # per-chrom hapmap (Chromosome, Position(bp), Rate(cM/Mb), Map(cM))
  chr2.hapmap.tsv
  all.hapmap.tsv       # combined hapmap across all chroms (what the Snakefile reads)
  sim.fai              # chrom name + length, for hapmap_low_rec_mask --fai
  ground_truth.json    # everything that was injected, for scoring
  README.md            # generated stub describing this specific run
```

Filenames and the per-chromosome subdirectory layout match what the
production pipeline's `infer_mu_path` and `find_tree_files` expect.
The combined `all.hapmap.tsv` and `sim.fai` are also written so the
dataset is ready to feed straight into the Snakemake workflow without
any post-processing.

## Running the pipeline on the generated dataset

```bash
snakemake --cores 4 --configfile config/snakemake.yaml \
  --config root_dir=<out-dir> \
           hapmap=<out-dir>/all.hapmap.tsv \
           fai=<out-dir>/sim.fai \
           tree_pattern="*.trees" \
           suffix_to_strip="" \
           burnin=0
```

Or copy `config/snakemake.yaml` to `<out-dir>/snakemake.yaml`, edit the
above fields once, and run `snakemake --configfile <out-dir>/snakemake.yaml`.

## CLI options

| Flag | Default | Effect |
|---|---|---|
| `--out-dir` | required | Where to write everything |
| `--n-chrom` | 2 | Chromosomes simulated |
| `--n-reps` | 5 | MCMC-like replicates per chromosome |
| `--n-samples` | 8 | Diploid individuals (= 2× sample nodes) |
| `--seq-length` | 5_000_000 | bp per chromosome |
| `--mutation-rate` | 3e-8 | Base mutation rate (per bp per generation) |
| `--recombination-rate` | 1e-8 | Recombination rate |
| `--n-contaminated` | 2 | Number of contaminated individuals |
| `--contam-severity` | `2,5,10` | Comma-separated severity multipliers, cycled across contaminated individuals |
| `--contam-hotspot-frac` | 0.0 | Fraction of each contaminated individual's extras placed in a per-(chrom, individual) hotspot (mimics paralog-mapping artifact). 0.0 = uniform; 0.8 puts 80% of extras in the hotspot |
| `--contam-hotspot-size` | 500000 | Width of each hotspot region in bp |
| `--prune-frac-windows` | 0.5 | Fraction of windows to prune |
| `--prune-frac-samples` | 0.25 | Per-window fraction of samples to prune |
| `--prune-window-size` | 200_000 | Window size for the prune mask in bp |
| `--mask-frac-genome` | 0.33 | Fraction of genome marked inaccessible |
| `--mask-n-intervals` | 3 | Number of contiguous mask intervals per chromosome |
| `--centromere` / `--no-centromere` | True | Include a low-rec valley in the hapmap |
| `--seed` | 42 | Base RNG seed; per-chrom and per-rep seeds derive deterministically from this |

Demography knobs (`--n-present`, `--n-ancestral`, `--n-bottleneck`,
`--expansion-end`, `--bot2-start`, `--bot2-end`, `--bot1-start`) match
`simulate_two_bottleneck_demography.py` and default to the same
two-bottleneck model documented there.

## Ground truth schema

```jsonc
{
  "seed": 42,
  "n_chrom": 2,
  "n_reps": 5,
  "n_samples_diploid": 8,
  "sequence_length": 5000000.0,
  "mutation_rate": 3e-08,
  "recombination_rate": 1e-08,
  "contaminated_individuals": [
    {"individual_id": 5, "severity": 2.0},
    {"individual_id": 6, "severity": 5.0}
  ],
  "contam_hotspot_frac": 0.0,
  "contam_hotspot_size": 500000.0,
  "chromosomes": [
    {
      "chrom": "chr1",
      "masked_intervals": [[117190.81, 199690.81], [367190.81, 449690.81]],
      "prune_intervals": [
        {"left": 0.0,     "right": 200000.0, "drop_sample_ids": [3, 9, 11, 15]},
        ...
      ],
      "contam_hotspots": [
        {"individual_id": 5, "left": 2350000.0, "right": 2850000.0},
        {"individual_id": 6, "left": 7100000.0, "right": 7600000.0}
      ]
    },
    ...
  ]
}
```

Contaminated individuals are global (same across all chromosomes and
replicates); masked intervals, prune intervals, and contamination
hotspots are per-chromosome but shared across replicates of that
chromosome. `contam_hotspots` is empty if `--contam-hotspot-frac` is 0.

## Scoring the pipeline

After running the Snakemake pipeline on a generated dataset, score it
against the ground truth:

```bash
python scripts/score_realistic_example.py \
  --ground-truth <out-dir>/ground_truth.json \
  --pipeline-out <pipeline_out_dir> \
  --report <pipeline_out_dir>/scoring_report.md
```

This emits both stdout text and a Markdown report covering:

- **Accessibility precision/recall** — `step2_low_access` BEDs vs
  `masked_intervals`.
- **Mutload outlier ranking** — per-chrom flag counts per individual
  from `step3_mutload/*.outliers.bed`; contaminated individuals should
  dominate by orders of magnitude.
- **Region-level mask vs prune intervals** — `mutation_masked.bed`
  overlap with the injected `prune_intervals`.
- **trim_samples recovery** — for each prune entry, count isolated
  samples at the window midpoint in the post-trim_samples ts.
- **Mutload per-(window, individual) FP rate** — split into total
  spurious, prune-explained spurious, and net FP (the actual
  statistical-noise floor of the per-window outlier test).

See the [example scoring report](#example-scoring-numbers-3-chrom-8-reps-16-dip-10-mb)
section below for representative numbers on a default-config run.

## Example scoring numbers (3 chrom × 8 reps × 16 dip × 10 Mb)

Default `make_realistic_example.py` config plus default
`config/snakemake.yaml` settings:

- **Accessibility**: precision 1.00, recall 0.96 — limited by the 50 kb
  step2 window grid not tiling the exact mask boundaries.
- **Mutload outlier individuals**: contaminated inds 1 and 12 are
  flagged ~9000-12000 times per chrom; the next most-flagged
  non-contaminated individual is flagged <50 times.
- **trim_samples recovery**: 100% — every pruned (rep, window, sample)
  is correctly isolated in the post-trim ts.
- **Mutload FP rate**: net FP rate (excluding prune-explained spurious)
  is **~0.17%** — about 1 in 600 flag-eligible cells is a false
  positive from natural Poisson noise alone.

Full numbers are in `argtest-realistic-example-out/scoring_report.md`
after running the scorer.

## Notes and limitations

- Injected mutations are placed at random positions across the whole
  genome (or within the hotspot when `--contam-hotspot-frac > 0`),
  **including** masked regions. This is a slight inconsistency with
  the "sim_mutations honors mut_rate.p" property of the natural
  mutations, but in practice mutload counts all carried mutations
  regardless of position, so it doesn't affect detection.
- Pruning happens after contamination injection. A small fraction of
  injected leaf mutations land on edges that subsequently get pruned —
  this reduces the effective severity slightly. For severity=5 the
  effect is small (~12% loss with default prune fractions); for
  severity=2 it can be enough to push some replicates below detection.
  If precise severity targeting matters, increase `--seq-length` or
  reduce `--prune-frac-windows`/`--prune-frac-samples`.
- The script writes `sim.fai` from the per-chromosome sequence lengths for use
  with `hapmap_low_rec_mask.py --fai`.
