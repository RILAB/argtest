#!/usr/bin/env python
"""Generate a realistic synthetic example dataset for end-to-end pipeline testing.

Base simulation: two-bottleneck demography (imported from
`simulate_two_bottleneck_demography.py`). Three flaws are injected
so the pipeline has something concrete to detect and clean up:

  1. **Extraneous mutations on K contaminated individuals.** Excess
     mutations are placed on the leaf edges of each contaminated sample
     node, scaled by a per-individual severity multiplier. Tests
     `mutload_masks.py`.
  2. **Per-window sample pruning.** A fraction of windows is selected;
     in each, a random subset of samples has its ancestry excised.
     Tests `trim_samples.py` semantics (and the partial-tree Ne fix).
  3. **Genome-wide accessibility mask.** A fraction of the genome is
     marked inaccessible in a `*.mut_rate.p` ratemap with rate=0 over
     those intervals. Tests `find_low_access_regions.py` and the
     mut-rate-aware sim mutload expectation.

A `ground_truth.json` records exactly what was injected so downstream
scoring can compute precision/recall on the pipeline's masks.

MCMC-like replication: across replicates within one chromosome, the
noise pattern (contaminated individuals, prune mask, accessibility
mask) is held constant while the underlying ARG is freshly simulated
with a different seed. Mimics how real singer MCMC output looks to
the pipeline.
"""
from __future__ import annotations

import argparse
import json
import pickle
import sys
from pathlib import Path

import numpy as np
import msprime
import tskit

# Reuse the demography builder + load helper from sibling project scripts.
sys.path.insert(0, str(Path(__file__).resolve().parent))
from simulate_two_bottleneck_demography import build_demography  # noqa: E402
from argtest_common import mutational_load  # noqa: E402


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out-dir", type=Path, required=True,
                   help="Directory to write the dataset into.")
    p.add_argument("--n-chrom", type=int, default=2,
                   help="Number of chromosomes (default: 2).")
    p.add_argument("--n-reps", type=int, default=5,
                   help="MCMC-like replicates per chromosome (default: 5).")
    p.add_argument("--n-samples", type=int, default=8,
                   help="Diploid individuals (default: 8 -> 16 sample nodes).")
    p.add_argument("--seq-length", type=float, default=5_000_000,
                   help="Sequence length per chromosome in bp (default: 5e6).")
    p.add_argument("--mutation-rate", type=float, default=3e-8,
                   help="Base mutation rate per bp per generation (default: 3e-8).")
    p.add_argument("--recombination-rate", type=float, default=1e-8,
                   help="Recombination rate per bp per generation (default: 1e-8).")
    p.add_argument("--n-contaminated", type=int, default=2,
                   help="Number of contaminated individuals (default: 2).")
    p.add_argument("--contam-severity", default="2,5,10",
                   help=("Comma-separated excess-mutation multipliers cycled across "
                         "contaminated individuals (default: '2,5,10')."))
    p.add_argument("--prune-frac-windows", type=float, default=0.5,
                   help="Fraction of windows to prune (default: 0.5).")
    p.add_argument("--prune-frac-samples", type=float, default=0.25,
                   help="Per-window fraction of samples to prune (default: 0.25).")
    p.add_argument("--prune-window-size", type=float, default=200_000.0,
                   help="Window size for the prune mask in bp (default: 2e5).")
    p.add_argument("--mask-frac-genome", type=float, default=0.33,
                   help="Fraction of genome to mark as inaccessible (default: 0.33).")
    p.add_argument("--mask-n-intervals", type=int, default=3,
                   help="Number of contiguous mask intervals per chromosome (default: 3).")
    p.add_argument("--centromere", action=argparse.BooleanOptionalAction, default=True,
                   help="Include a low-rec centromere valley in the hapmap (default: True).")
    p.add_argument("--seed", type=int, default=42,
                   help="Base random seed (default: 42).")
    # Demography knobs, default to the two-bottleneck model.
    p.add_argument("--n-present", type=int, default=200_000)
    p.add_argument("--n-ancestral", type=int, default=100_000)
    p.add_argument("--n-bottleneck", type=int, default=10_000)
    p.add_argument("--expansion-end", type=float, default=9000)
    p.add_argument("--bot2-start", type=float, default=15000)
    p.add_argument("--bot2-end", type=float, default=30000)
    p.add_argument("--bot1-start", type=float, default=35000)
    return p.parse_args()


def build_mu_ratemap(seq_length: float, masked_intervals, base_rate: float) -> msprime.RateMap:
    # Construct a stepwise ratemap that is base_rate over accessible
    # regions and 0 over masked intervals. `masked_intervals` is a
    # sorted list of (left, right) pairs that don't overlap.
    boundaries = [0.0]
    rates = []
    cur = 0.0
    for left, right in masked_intervals:
        if left > cur:
            boundaries.append(float(left))
            rates.append(base_rate)
        boundaries.append(float(right))
        rates.append(0.0)
        cur = float(right)
    if cur < seq_length:
        boundaries.append(float(seq_length))
        rates.append(base_rate)
    elif boundaries[-1] < seq_length:
        boundaries.append(float(seq_length))
        rates.append(0.0)
    return msprime.RateMap(position=boundaries, rate=rates)


def pick_mask_intervals(seq_length: float, frac: float, n_intervals: int,
                        rng: np.random.Generator):
    # Pick n_intervals non-overlapping intervals whose total length is
    # approximately `frac * seq_length`. Each interval gets ~ equal
    # length; placements are randomized but kept disjoint.
    total = frac * seq_length
    if total <= 0 or n_intervals <= 0:
        return []
    n_intervals = max(1, int(n_intervals))
    length_each = total / n_intervals
    # Slice the genome into n_intervals + 1 gaps; place an interval at the
    # start of each picked slot.
    gap_each = (seq_length - total) / n_intervals
    if gap_each < 0:
        # Mask is larger than the genome; clamp to one big interval.
        return [(0.0, float(seq_length))]
    starts = []
    cursor = rng.uniform(0, gap_each)
    for _ in range(n_intervals):
        starts.append(cursor)
        cursor += length_each + gap_each
    intervals = [(float(s), float(min(s + length_each, seq_length))) for s in starts]
    intervals.sort()
    return intervals


def pick_prune_windows(seq_length: float, window_size: float, frac_windows: float,
                       n_samples_total: int, frac_samples: float,
                       sample_ids: list[int], rng: np.random.Generator):
    # Build a list of (left, right, drop_sample_ids) for the per-window
    # pruning pattern. Same recipe as the bottleneck-test harness.
    n_windows = max(1, int(np.ceil(seq_length / window_size)))
    chosen = rng.choice(
        n_windows,
        size=max(1, int(round(frac_windows * n_windows))),
        replace=False,
    )
    n_drop = max(1, int(round(frac_samples * n_samples_total)))
    intervals = []
    for wi in sorted(chosen):
        left = float(wi * window_size)
        right = float(min(left + window_size, seq_length))
        drop = rng.choice(sample_ids, size=n_drop, replace=False).tolist()
        intervals.append((left, right, [int(x) for x in drop]))
    return intervals


def apply_prune(ts, prune_intervals):
    # Drop edges whose child is a target sample within a pruned window.
    # No simplify, so sample ids remain stable across chroms / reps.
    tables = ts.dump_tables()
    positions = sorted(
        set(float(l) for l, _, _ in prune_intervals)
        | set(float(r) for _, r, _ in prune_intervals)
    )
    for pos in positions:
        n_initial = tables.edges.num_rows
        for i in range(n_initial):
            e = tables.edges[i]
            if e.left < pos < e.right:
                tables.edges[i] = e.replace(right=pos)
                tables.edges.append(e.replace(left=pos))
    keep = np.ones(tables.edges.num_rows, dtype=bool)
    left_arr = tables.edges.left
    right_arr = tables.edges.right
    child_arr = tables.edges.child
    for l, r, drop in prune_intervals:
        if not drop:
            continue
        within = (left_arr >= l) & (right_arr <= r)
        in_drop = np.isin(child_arr, np.asarray(drop, dtype=int))
        keep &= ~(within & in_drop)
    tables.edges.keep_rows(keep)
    tables.sort()
    tables.build_index()
    tables.compute_mutation_parents()
    return tables.tree_sequence()


def inject_contamination(ts, contam_indiv_ids, severities, rng: np.random.Generator):
    # For each contaminated individual, drop extra leaf mutations onto a
    # random sample node of that individual. Baseline is the actual mean
    # per-individual load (computed from mutational_load and aggregated to
    # individual level), so severity=k aims for final load ~ k × mean load
    # on that individual.
    if not contam_indiv_ids:
        return ts
    n_indiv = ts.num_individuals
    if n_indiv == 0:
        return ts
    sample_load = mutational_load(ts)
    per_ind_load = np.zeros(n_indiv, dtype=float)
    for u, val in enumerate(sample_load):
        ind_id = ts.node(u).individual
        if ind_id != tskit.NULL:
            per_ind_load[ind_id] += float(val)
    baseline = float(per_ind_load.mean()) if n_indiv > 0 else 0.0
    tables = ts.dump_tables()
    site_positions = set(float(p) for p in np.asarray(tables.sites.position))
    indiv_nodes = {int(i): [int(u) for u in ts.individual(i).nodes]
                   for i in contam_indiv_ids}
    for indiv_id, sev in zip(contam_indiv_ids, severities):
        excess = max(0, sev - 1.0)
        n_extra = int(rng.poisson(excess * baseline))
        nodes = indiv_nodes[int(indiv_id)]
        if n_extra == 0 or not nodes:
            continue
        for _ in range(n_extra):
            pos = float(rng.uniform(0.0, ts.sequence_length))
            # Jitter off any exact collisions with existing sites.
            while pos in site_positions:
                pos = min(pos + 1e-9, ts.sequence_length - 1e-9)
            site_positions.add(pos)
            u = int(nodes[rng.integers(0, len(nodes))])
            site_id = tables.sites.add_row(position=pos, ancestral_state="0")
            tables.mutations.add_row(site=site_id, node=u, derived_state="1")
    tables.sort()
    tables.build_index()
    tables.compute_mutation_parents()
    return tables.tree_sequence()


def write_hapmap(out_path: Path, chrom_name: str, seq_length: float,
                 base_rate_cm_per_mb: float, centromere: bool,
                 rng: np.random.Generator):
    # Tab-separated, header Chromosome / Position(bp) / Rate(cM/Mb) / Map(cM).
    # If centromere is True, drop the rate to ~5% of baseline over the
    # middle 10% of the chromosome — gives hapmap_low_rec_mask.py
    # something to detect.
    step = max(10_000.0, seq_length / 200.0)  # ~200 anchor rows
    positions = np.arange(0, seq_length + step, step)
    if positions[-1] > seq_length:
        positions[-1] = seq_length
    if centromere:
        cen_left = 0.45 * seq_length
        cen_right = 0.55 * seq_length
    else:
        cen_left = cen_right = -1.0
    rates_cm_per_mb = np.full(len(positions), base_rate_cm_per_mb, dtype=float)
    for i, p in enumerate(positions):
        if cen_left <= p <= cen_right:
            rates_cm_per_mb[i] = base_rate_cm_per_mb * 0.05
    # Cumulative map distance in cM. cM = (Rate cM/Mb) * (segment Mb).
    map_cm = np.zeros(len(positions), dtype=float)
    for i in range(1, len(positions)):
        seg_mb = (positions[i] - positions[i - 1]) / 1e6
        # Each segment uses the rate at the LEFT boundary (segment starts there).
        map_cm[i] = map_cm[i - 1] + rates_cm_per_mb[i - 1] * seg_mb
    lines = ["Chromosome\tPosition(bp)\tRate(cM/Mb)\tMap(cM)"]
    for p, r, m in zip(positions, rates_cm_per_mb, map_cm):
        lines.append(f"{chrom_name}\t{int(p)}\t{r:.4f}\t{m:.6f}")
    out_path.write_text("\n".join(lines) + "\n")


def main():
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    severities = [float(s) for s in args.contam_severity.split(",") if s.strip()]
    if not severities:
        raise ValueError("--contam-severity must list at least one multiplier.")

    demography = build_demography(
        n_present=args.n_present,
        n_ancestral=args.n_ancestral,
        n_bottleneck=args.n_bottleneck,
        expansion_end=args.expansion_end,
        bot2_start=args.bot2_start,
        bot2_end=args.bot2_end,
        bot1_start=args.bot1_start,
    )

    # Global noise: which individuals are contaminated and at what severity.
    # Deterministic so the same individuals are dirty in every chrom / rep.
    global_rng = np.random.default_rng(args.seed)
    n_diploid = int(args.n_samples)
    contam_ids = global_rng.choice(
        n_diploid, size=min(args.n_contaminated, n_diploid), replace=False,
    ).tolist()
    contam_ids.sort()
    contam_severity = [severities[i % len(severities)] for i in range(len(contam_ids))]

    truth = {
        "seed": int(args.seed),
        "n_chrom": int(args.n_chrom),
        "n_reps": int(args.n_reps),
        "n_samples_diploid": n_diploid,
        "sequence_length": float(args.seq_length),
        "mutation_rate": float(args.mutation_rate),
        "recombination_rate": float(args.recombination_rate),
        "contaminated_individuals": [
            {"individual_id": int(i), "severity": float(s)}
            for i, s in zip(contam_ids, contam_severity)
        ],
        "chromosomes": [],
    }

    for chrom_i in range(args.n_chrom):
        chrom_name = f"chr{chrom_i + 1}"
        chrom_dir = args.out_dir / chrom_name
        chrom_dir.mkdir(parents=True, exist_ok=True)

        # Chrom-specific RNG so changing --n-chrom doesn't reshuffle prior chroms.
        chrom_rng = np.random.default_rng(args.seed * 7919 + chrom_i)

        # 1) Accessibility mask -> mut_rate.p
        masked_intervals = pick_mask_intervals(
            args.seq_length, args.mask_frac_genome, args.mask_n_intervals, chrom_rng,
        )
        mu_ratemap = build_mu_ratemap(args.seq_length, masked_intervals, args.mutation_rate)
        with open(args.out_dir / f"{chrom_name}.mut_rate.p", "wb") as fh:
            pickle.dump(mu_ratemap, fh)

        # 2) Hapmap with optional low-rec centromere.
        write_hapmap(
            args.out_dir / f"{chrom_name}.hapmap.tsv",
            chrom_name=chrom_name,
            seq_length=args.seq_length,
            base_rate_cm_per_mb=1.0,
            centromere=args.centromere,
            rng=chrom_rng,
        )

        # 3) Per-window prune pattern. Sample ids are 0..2*n_diploid-1.
        sample_ids = list(range(2 * n_diploid))
        prune_intervals = pick_prune_windows(
            args.seq_length, args.prune_window_size,
            args.prune_frac_windows, len(sample_ids), args.prune_frac_samples,
            sample_ids, chrom_rng,
        )

        # 4) MCMC-like replicates: fresh ARG topology per rep, shared noise pattern.
        for rep_i in range(args.n_reps):
            rep_seed = int(args.seed * 1_000_003 + chrom_i * 1009 + rep_i)
            ts = msprime.sim_ancestry(
                samples={"pop": n_diploid},
                demography=demography,
                sequence_length=args.seq_length,
                recombination_rate=args.recombination_rate,
                random_seed=rep_seed,
            )
            # Use the ratemap (with zeros over masked regions) so observed
            # mutations honor the accessibility mask.
            ts = msprime.sim_mutations(
                ts, rate=mu_ratemap, random_seed=rep_seed + 1, keep=False,
            )
            # Inject contamination (rep-specific noise so the extra count varies
            # across MCMC iterates, like inference noise).
            inject_rng = np.random.default_rng(rep_seed + 2)
            ts = inject_contamination(ts, contam_ids, contam_severity, inject_rng)
            # Apply per-window pruning.
            ts = apply_prune(ts, prune_intervals)
            out_ts = chrom_dir / f"rep_{rep_i:03d}.trees"
            ts.dump(out_ts)

        truth["chromosomes"].append({
            "chrom": chrom_name,
            "masked_intervals": [[float(l), float(r)] for l, r in masked_intervals],
            "prune_intervals": [
                {"left": float(l), "right": float(r), "drop_sample_ids": d}
                for l, r, d in prune_intervals
            ],
        })

    truth_path = args.out_dir / "ground_truth.json"
    truth_path.write_text(json.dumps(truth, indent=2) + "\n")

    readme = args.out_dir / "README.md"
    readme.write_text(_README_TEMPLATE.format(
        out_dir=args.out_dir.name,
        n_chrom=args.n_chrom,
        n_reps=args.n_reps,
        n_samples=n_diploid,
        seq_length=int(args.seq_length),
        contam_ids=", ".join(str(i) for i in contam_ids),
        contam_severity=", ".join(f"{s:g}" for s in contam_severity),
        mask_frac=args.mask_frac_genome,
        prune_windows=args.prune_frac_windows,
        prune_samples=args.prune_frac_samples,
        seed=args.seed,
    ))

    print(f"Wrote {args.n_chrom} chroms x {args.n_reps} reps to {args.out_dir}")
    print(f"Contaminated individuals: {contam_ids}  severities: {contam_severity}")
    print(f"Ground truth: {truth_path}")


_README_TEMPLATE = """# {out_dir}

Synthetic example dataset generated by `scripts/make_realistic_example.py`.

- {n_chrom} chromosomes x {n_reps} replicates
- {n_samples} diploid individuals ({n_samples} x 2 = sample nodes per chrom)
- {seq_length} bp per chromosome
- Two-bottleneck demography (see scripts/simulate_two_bottleneck_demography.py)
- seed = {seed}

## Injected flaws

1. Contaminated individuals: {contam_ids} with severity multipliers {contam_severity}.
2. Per-window sample pruning: {prune_samples:.0%} of samples pruned in {prune_windows:.0%}
   of windows (per-chrom RNG).
3. Accessibility mask: ~{mask_frac:.0%} of genome zeroed in chrN.mut_rate.p.

See ground_truth.json for the exact masked intervals and prune sets per chrom.

## Layout

- `chrN/rep_NNN.trees`  -- MCMC-like replicates (fresh ARG topology per rep)
- `chrN.mut_rate.p`      -- msprime.RateMap (rate=0 over inaccessible intervals)
- `chrN.hapmap.tsv`      -- recombination map (low-rec valley around midpoint)
- `ground_truth.json`    -- what was injected, for precision/recall scoring
"""


if __name__ == "__main__":
    main()
