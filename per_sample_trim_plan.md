# Per-sample (per-node) trim granularity — design notes (deferred)

Planning notes for project todo **#6** ("Make `trim_samples.py` removal
granularity explicitly per-sample (per-node), robust to ploidy"). Deferred —
this is a bigger change than a quick edit. Verified against the code on
2026-06-26.

## Goal

Make sample removal explicitly **per sample node (per haplotype)**, robust to
ploidy. Today removal is keyed by **individual name** and expands to *all* of
that individual's nodes — fine on the haploid `admix` data (1 node/individual),
but on a truly diploid/polyploid tree sequence it removes **every** haplotype of
the individual (= the whole diploid), silently.

## Verified current behavior (2026-06-26)

- **Detection — step 3, `scripts/mutload_masks.py`:** observed load is aggregated
  to individuals via `aggregate_by_individual(load, names)` (line 143), and the
  sim-based expected matrix likewise (line 122). The outlier test is per
  `(window, individual)`. BED **column 4** is the comma-joined list of
  individual names `unique_names` (line 179); columns 5/6 are the matching
  observed / expected loads.
- **Trim — step 5, `scripts/trim_samples.py`:** `name_to_nodes_map`
  (`argtest_common.py:97`) maps a name → `ts.individual(id).nodes` (**all** nodes
  of the individual). `trim_samples_single_pass` (line 126+) does
  `nodes = name_to_nodes.get(name, [])` and removes ancestry for **all** of them
  over the interval.
- **Net:** the chain is per-**individual** end to end. "Per-sample" behavior only
  coincides on haploid data where each individual has exactly one node.

## Design forks (need a decision before implementing)

1. **BED col-4 label scheme** (the originally-flagged decision):
   - *node ids (ints)* — rejected: node ids are not stable across tree sequences
     (they differ per chrom/replicate), so a label can't be carried between files.
   - *`name#haplotype`* (e.g. `ind5#0`, `ind5#1`) — **recommended**: stable,
     human-readable, maps to the k-th sample node of the individual. A **bare
     `name`** (no `#`) keeps meaning "all nodes of that individual," so haploid
     `admix` output is byte-for-byte unchanged (backward compatible).

2. **Detection granularity on diploid** (the subtler issue): the mutload test
   currently **pools** an individual's haplotypes *before* the outlier test. To
   flag/trim a single haplotype, step 3 must be able to identify a single node —
   i.e. **not pool** (test per sample node), which is a statistical change that
   also affects the expected-load sim aggregation. Decide between:
   - keep per-individual detection and trim all nodes (status quo — not actually
     per-sample), or
   - test per sample node (no pooling) so an individual haplotype can be flagged
     and trimmed. Likely wants a config knob to select pooled vs per-node mode.

3. **`suffix_to_strip` / `get_individual_name` interaction:** ensure `#hap`
   parsing happens *after* suffix stripping so the two don't collide.

## Backward-compatibility requirement

On haploid `admix`-style data (1 node/individual) bare-name labels must continue
to mean "the single node" and produce identical output. Add a regression test on
a haploid fixture asserting no change.

## Implementation sketch (once decisions are made)

- `argtest_common.py`: add a `name_to_node` resolver that understands `name#k`
  → the k-th node of the individual; keep `name_to_nodes_map` for bare names.
- `mutload_masks.py`: emit col-4 labels at the chosen granularity (node-level
  `name#k` when not pooling).
- `trim_samples.py`: parse col-4 labels and map each to specific node(s).
- `Snakefile`: no interface change unless a config knob is added for pooling mode.
- **Tests:** extend the trim tests to a synthetic **diploid** ts (2 nodes/
  individual) and assert only the targeted haplotype is trimmed; keep a haploid
  regression test.

## Related

- `project_trim_samples_bias` (memory) — coalescence-rate / Ne bias from
  partial-sample trees; the same low-sample-tree problem motivates todo #5.
- `project_kept_intervals` (memory) — any tree-mutating filter must consider
  whether to rewrite `kept_intervals` metadata.
- Todo **#5** (`--min-samples` filter) overlaps in motivation (both guard against
  low-sample trees distorting stats), but operates at a different layer.
