#!/usr/bin/env python
"""Simulate replicate ARGs under a two-bottleneck + expansion demography.

Forward-time history (assuming 1 year per generation):
    (-inf, 35000 ya] : N = 100_000   (ancestral)
    (35000, 30000] ya : N = 10_000   (bottleneck #1)
    (30000, 15000] ya : N = 100_000  (recovery)
    (15000,  9000] ya : N = 10_000   (bottleneck #2)
    (9000,    0] ya   : exponential growth from 100_000 to 200_000 (present)

In msprime backward-time terms this is one population whose size and
growth rate change at t in {9000, 15000, 30000, 35000} (generations
before present).
"""
from __future__ import annotations

import argparse
import math
from pathlib import Path

import msprime


def build_demography(
    n_present: int,
    n_ancestral: int,
    n_bottleneck: int,
    expansion_end: float,
    bot2_start: float,
    bot2_end: float,
    bot1_start: float,
) -> msprime.Demography:
    """Build the demography. Time arguments are forward-time labels in
    generations before present (e.g. `bot2_start` is the forward-time start of
    the more-recent bottleneck), but msprime events are scheduled in backward
    time. add_population_parameters_change(time=T, initial_size=N) sets the
    size for the epoch [T, next_T) — i.e., for times deeper into the past
    starting at T.
    """
    if not (0 < expansion_end < bot2_start < bot2_end < bot1_start):
        raise ValueError(
            "Time points must satisfy 0 < expansion_end < bot2_start < "
            "bot2_end < bot1_start."
        )
    # Want N(t=expansion_end) = n_ancestral going backward, where
    # N(t) = n_present * exp(-g * t). Solve for g.
    growth_rate = math.log(n_present / n_ancestral) / expansion_end

    d = msprime.Demography()
    d.add_population(name="pop", initial_size=n_present, growth_rate=growth_rate)
    # Backward t=expansion_end (forward: 9 ka): enter bottleneck #2.
    d.add_population_parameters_change(
        time=expansion_end,
        initial_size=n_bottleneck,
        growth_rate=0,
        population="pop",
    )
    # Backward t=bot2_start (forward: 15 ka): exit bottleneck #2 going back.
    d.add_population_parameters_change(
        time=bot2_start,
        initial_size=n_ancestral,
        growth_rate=0,
        population="pop",
    )
    # Backward t=bot2_end (forward: 30 ka): enter bottleneck #1.
    d.add_population_parameters_change(
        time=bot2_end,
        initial_size=n_bottleneck,
        growth_rate=0,
        population="pop",
    )
    # Backward t=bot1_start (forward: 35 ka): exit bottleneck #1; deep
    # ancestral N continues to infinity.
    d.add_population_parameters_change(
        time=bot1_start,
        initial_size=n_ancestral,
        growth_rate=0,
        population="pop",
    )
    return d


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out-dir", type=Path, required=True,
                   help="Directory to write replicate .trees files into.")
    p.add_argument("--n-replicates", type=int, default=10,
                   help="Number of independent ARG replicates (default: 10).")
    p.add_argument("--samples", type=int, default=20,
                   help="Diploid samples (default: 20).")
    p.add_argument("--sequence-length", type=float, default=1e7,
                   help="Sequence length in bp (default: 10 Mb).")
    p.add_argument("--recombination-rate", type=float, default=1e-8,
                   help="Recombination rate per bp per gen (default: 1e-8).")
    p.add_argument("--mutation-rate", type=float, default=3e-8,
                   help="Mutation rate per bp per gen for sim_mutations "
                        "(default: 3e-8). Set to 0 to skip adding mutations.")
    p.add_argument("--n-present", type=int, default=200_000,
                   help="Present-day population size (default: 200000).")
    p.add_argument("--n-ancestral", type=int, default=100_000,
                   help="Ancestral / inter-bottleneck size (default: 100000).")
    p.add_argument("--n-bottleneck", type=int, default=10_000,
                   help="Size during each bottleneck (default: 10000).")
    p.add_argument("--expansion-end", type=float, default=9000,
                   help="Generations ago when expansion phase begins going "
                        "forward (default: 9000).")
    p.add_argument("--bot2-start", type=float, default=15000,
                   help="Generations ago when bottleneck #2 began going "
                        "forward (default: 15000).")
    p.add_argument("--bot2-end", type=float, default=30000,
                   help="Generations ago when bottleneck #1 ended going "
                        "forward (default: 30000).")
    p.add_argument("--bot1-start", type=float, default=35000,
                   help="Generations ago when bottleneck #1 began going "
                        "forward (default: 35000).")
    p.add_argument("--seed", type=int, default=42,
                   help="Base random seed (default: 42).")
    p.add_argument("--prefix", default="rep",
                   help="Output filename prefix (default: 'rep').")
    p.add_argument("--demes-out", type=Path, default=None,
                   help="Optional path to dump the resolved demography as a "
                        "Demes YAML file.")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    # Note: in build_demography these are interpreted as msprime backward-time
    # event times. The script's argparse names match the forward-time labels,
    # but the order of times we pass to msprime must be increasing.
    demography = build_demography(
        n_present=args.n_present,
        n_ancestral=args.n_ancestral,
        n_bottleneck=args.n_bottleneck,
        expansion_end=args.expansion_end,
        bot2_start=args.bot2_start,
        bot2_end=args.bot2_end,
        bot1_start=args.bot1_start,
    )

    if args.demes_out is not None:
        import demes
        graph = demography.to_demes()
        args.demes_out.parent.mkdir(parents=True, exist_ok=True)
        demes.dump(graph, str(args.demes_out))
        print(f"Wrote demography to {args.demes_out}")

    for i in range(args.n_replicates):
        seed = args.seed + 2 * i
        ts = msprime.sim_ancestry(
            samples={"pop": args.samples},
            demography=demography,
            sequence_length=args.sequence_length,
            recombination_rate=args.recombination_rate,
            random_seed=seed,
        )
        if args.mutation_rate > 0:
            ts = msprime.sim_mutations(
                ts, rate=args.mutation_rate, random_seed=seed + 1, keep=False
            )
        out_path = args.out_dir / f"{args.prefix}_{i:03d}.trees"
        ts.dump(out_path)
        print(f"Wrote {out_path} (trees={ts.num_trees}, sites={ts.num_sites})")


if __name__ == "__main__":
    main()
