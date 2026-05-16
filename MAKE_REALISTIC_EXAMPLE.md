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
| Extraneous mutations on K individuals | Extra leaf-edge mutations sized by `(severity − 1) × mean per-individual load`; placed at random positions, assigned to a random sample node of the contaminated individual | `mutload_masks.py` outlier detection |
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
  chr1.hapmap.tsv      # standard hapmap (Chromosome, Position(bp), Rate(cM/Mb), Map(cM))
  chr2.hapmap.tsv
  ground_truth.json    # everything that was injected, for scoring
  README.md            # generated stub describing this specific run
```

Filenames and the per-chromosome subdirectory layout match what the
production pipeline's `infer_mu_path` and `find_tree_files` expect, so
the dataset can be fed directly into the Snakemake workflow.

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
  "chromosomes": [
    {
      "chrom": "chr1",
      "masked_intervals": [[117190.81, 199690.81], [367190.81, 449690.81]],
      "prune_intervals": [
        {"left": 0.0,     "right": 200000.0, "drop_sample_ids": [3, 9, 11, 15]},
        ...
      ]
    },
    ...
  ]
}
```

Contaminated individuals are global (same across all chromosomes and
replicates); masked intervals and prune intervals are per-chromosome
but shared across replicates of that chromosome.

## Scoring the pipeline (sketch)

Given a pipeline run's mask BEDs and trimmed tree sequences, you can
score against ground truth:

- **Accessibility precision/recall** — intersect the pipeline's
  low-accessibility BED with `masked_intervals` for that chromosome.
- **Contamination precision/recall** — for each chromosome, take the
  union of per-replicate mutload-flagged individual IDs (from
  `mutload_masks.py` output) and compare to
  `contaminated_individuals`.
- **Pruning recovery** — for each `prune_intervals` window, check
  whether the post-trim_samples tree sequences have the expected samples
  isolated over that interval (`tree.parent(sample) == NULL`).

These would naturally live in a separate scoring script (not yet
written).

## Notes and limitations

- Injected mutations are placed at random positions across the whole
  genome, **including** masked regions. This is a slight inconsistency
  with the "sim_mutations honors mut_rate.p" property of the natural
  mutations, but in practice mutload counts all carried mutations
  regardless of position, so it doesn't affect detection.
- Pruning happens after contamination injection. A small fraction of
  injected leaf mutations land on edges that subsequently get pruned —
  this reduces the effective severity slightly. For severity=5 the
  effect is small (~12% loss with default prune fractions); for
  severity=2 it can be enough to push some replicates below detection.
  If precise severity targeting matters, increase `--seq-length` or
  reduce `--prune-frac-windows`/`--prune-frac-samples`.
- The script does **not** currently write a `.fai` file; if you need
  one for `hapmap_low_rec_mask.py --fai`, generate it from the per-chrom
  sequence lengths.
