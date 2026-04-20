# Mutational load trace: singer-snakemake vs argtest

Side-by-side comparison of how each codebase calculates and plots the per-sample
mutational load trace across MCMC replicates.

---

## 1. The `mutational_load` function

Both codebases use **identical logic**.

| | singer-snakemake (`utils.py`) | argtest (`scripts/argtest_common.py`) |
|---|---|---|
| Iterate sites | `ts.first()` + `tree.seek(pos)` | `sub.trees()` + `tree.sites()` |
| Filter | `m.edge != tskit.NULL` | `m.edge != tskit.NULL` |
| Count | `tree.samples(m.node)` | `tree.samples(m.node)` |
| Extra capability | — | optional `remove_intervals` / `name_to_nodes` for sample trimming |

**singer-snakemake `utils.py`:**
```python
def mutational_load(ts, windows=None):
    genome_windows = np.array([0, ts.sequence_length]) if windows is None else windows
    mutations_window = np.digitize(ts.sites_position[ts.mutations_site], genome_windows) - 1
    load = np.zeros((genome_windows.size - 1, ts.num_samples))
    tree = ts.first(sample_lists=True)
    for s in ts.sites():
        tree.seek(s.position)
        for m in s.mutations:
            if m.edge != tskit.NULL:
                window = mutations_window[m.id]
                samples = list(tree.samples(m.node))
                load[window, samples] += 1.0
    return load.squeeze(0) if windows is None else load
```

**argtest `argtest_common.py`:**
```python
def mutational_load(ts, windows=None, remove_intervals=None, name_to_nodes=None):
    genome_windows = np.array([0, ts.sequence_length]) if windows is None else windows
    load = np.zeros((genome_windows.size - 1, ts.num_samples))
    segments = [(0.0, ts.sequence_length, frozenset())]
    if remove_intervals and name_to_nodes:
        segments = build_removal_segments(remove_intervals, ts.sequence_length)
    for left, right, drop_names in segments:
        sub = ts.keep_intervals([(left, right)], simplify=False)
        # ... optional per-sample drop logic ...
        for tree in sub.trees(sample_lists=True):
            for s in tree.sites():
                for m in s.mutations:
                    if m.edge != tskit.NULL:
                        samples = list(tree.samples(m.node))
                        load[window, samples] += 1.0
    return load.squeeze(0) if windows is None else load
```

**Verified numerically identical** on 150 replicates of `amaranth/amaranth.1`:
global max |argtest − singer| = 0.0

---

## 2. Normalization

| | singer-snakemake (`tree_statistics.py`) | argtest (`validation_plots_from_ts.py`) |
|---|---|---|
| Denominator | `np.sum(accessible_bp)` — **per replicate** | **per replicate** for cleaned trees; **fixed** for original trees |
| Source | `inaccessible.p` ratemap × `extract_accessible_ratemap(trees)` | `kept_intervals` metadata (per-replicate, set by step 4) → `mu_intervals` from `.mut_rate.p` → `ts.sequence_length` |
| Tree-coverage adjustment | Yes: empty tree intervals (no edges) excluded each replicate | No |

**singer-snakemake `tree_statistics.py`:**
```python
accessible = msprime.RateMap(position=inaccessible.position, rate=1 - inaccessible.rate)
accessible = multiply_ratemaps(accessible, extract_accessible_ratemap(trees))
accessible_bp = np.diff(accessible.get_cumulative_mass(windows.position))

observed_load = mutational_load(trees)
observed_load /= np.sum(accessible_bp)   # denominator varies per replicate
```

**argtest `validation_plots_from_ts.py`:**
```python
# Priority: kept_intervals > mu_intervals from ratemap > sequence_length
kept_intervals = ts.metadata.get("kept_intervals") if ts.metadata else None
if kept_intervals is not None:
    total_accessible = float(sum(r - l for l, r in kept_intervals))
elif mu_intervals is not None:
    total_accessible = float(np.sum(mu_intervals[:, 1] - mu_intervals[:, 0]))
else:
    total_accessible = float(ts.sequence_length)

load = mutational_load(ts) / total_accessible   # same denominator every replicate
```

---

## 3. Plotting the trace

| | singer-snakemake (`diagnostics.py`) | argtest (`validation_plots_from_ts.py`) |
|---|---|---|
| X-axis | MCMC iteration (`replicate_index × mcmc_thin`) | Replicate file index (0-based) |
| X-axis scale | 0 – 1490 (with `mcmc_thin=10`) | 0 – 149 |
| Per-sample lines | `for i, x in enumerate(observed_load)` | `for i in range(load_vals.shape[0])` |
| Burnin marker | none on trace | vertical dashed line at `burnin` |
| Comparison overlay | none | optional second dataset (original vs cleaned) |

**singer-snakemake `diagnostics.py`:**
```python
mcmc_iterates = np.arange(num_mcmc) * mcmc_thin   # e.g. 0, 10, 20, ..., 1490

for i, x in enumerate(observed_load):             # observed_load: [sample, rep]
    plt.plot(mcmc_iterates, x, "-", color="black", linewidth=0.5, alpha=0.25)
plt.xlabel("MCMC iteration")
plt.ylabel("Derived mutations / base in each sample")
```

**argtest `validation_plots_from_ts.py`:**
```python
for i in range(pri["load_vals"].shape[0]):         # load_vals: [sample, rep]
    ax.plot(pri["keep_idx"], pri["load_vals"][i], "-", color=PRI_SITE,
            linewidth=0.5, alpha=0.2)
ax.axvline(pri["burnin"], color=PRI_BRANCH, linestyle="--", linewidth=1)
ax.set_xlabel("Replicate index")
ax.set_ylabel("Derived mutations / base in each sample")
```

---

## Summary

The two plots are **visually identical** when run on the same data. Confirmed on
`amaranth/amaranth.1` (150 replicates, 14 samples): the mutation counts are
bit-for-bit equal. Apparent differences between the plots are explained by:

1. **X-axis label only** — replicate index vs MCMC iteration (factor of 10)
2. **Normalization source** — argtest uses a fixed accessible-bp denominator;
   singer-snakemake recomputes it per replicate including tree-coverage masking.
   In practice this produces the same y-axis range on this dataset.
