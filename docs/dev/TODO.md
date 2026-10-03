# TODO

Open work items. See `CHANGELOG.md` for released changes and
`REVIEW_FOLLOWUP.md` for the review findings these came from.

---

## 1. Benchmark merge peak RSS — may require reverting a v1.10 change

**Status:** open · **Blocked on:** HPC + real ARGs · **Could trigger:** v1.11 revert

v1.10 changed `merge_group` from incremental concatenation
(`merged = merged.concatenate(ts)` in a loop) to a single batched
`first.concatenate(*rest)` call, on the theory that it avoids materialising a
progressively larger intermediate. **This was never measured on real data.**

Batching can plausibly be *worse*: it holds every input tree sequence alive
simultaneously, whereas the incremental loop releases each one as it goes. The
theoretical peak is similar either way (incremental holds
`merged_{N-1} + merged_N + ts_N`; batched holds all inputs plus the output), so
any win comes from tskit's internal allocation behaviour — which is exactly why
it needs measuring rather than reasoning about. Codex independently flagged the
same risk during review.

`measure_merge_mem.py` now supports both strategies. Run each as its **own
process** — `ru_maxrss` is a high-water mark and never decreases, so measuring
both in one process reports only the larger:

```
HPC_MEM=192G ~/.claude/bin/hpc_run 'python measure_merge_mem.py --mode incremental --rep 0 --glob "<one replicate, all chroms>"'
HPC_MEM=192G ~/.claude/bin/hpc_run 'python measure_merge_mem.py --mode batched     --rep 0 --glob "<same glob>"'
```

Use post-step-5/5b tree sequences — that is what `merge_replicates` actually
consumes. Budget generously; per `CLAUDE.md` the merge peak is far higher than
any other step.

**Then:**
- Batched wins or ties → record the numbers in the v1.10 changelog entry, which
  currently makes *no* performance claim precisely because this is unmeasured.
- Incremental wins → revert `merge_group` to the loop and ship v1.11.

---

## 2. Run the ARG contract audit on the real corpus — a live bug may be hiding

**Status:** open · **Blocked on:** real amaranth/admix ARGs · **Cheap:** read-only

`scripts/audit_arg_contract.py` has only ever run against fixtures. It is the
evidence gate for every deferred enforcement decision (uniform ploidy, leaf
sample nodes, removing the VCF one-column-per-node fallback) — but more
urgently, it is the only way to detect a **live correctness bug**:

`name_to_nodes_map` is last-wins on duplicate normalized individual names, while
`aggregate_by_individual` *pools* them. If two individuals collide after
`name_substring_to_remove` is applied, step 3 flags a pooled individual and step
5 trims only one of them — silently, with no error. With the default removal
string `_anchorwave`, a collision is plausible if any two sample IDs differ only
by that substring.

v1.10 made the audit *report* duplicate names but deliberately left the
underlying behaviour in place, pending this evidence.

```
python scripts/audit_arg_contract.py --name-substring-to-remove _anchorwave <real ARGs...>
```

**Then:**
- Clean → record it, and the deferred contract enforcement becomes safe to do.
- Duplicate names reported → this is a real bug; fix `name_to_nodes_map` to
  merge rather than overwrite, or make duplicates a hard error.

---

## Also pending

- **`v1.10` is tagged locally but not pushed.** `git push origin v1.10` when
  ready.
- **Gate 1 (mutation-map discovery) is closeable any time** — not blocking, but
  unrun on real data. `python scripts/check_mu_paths.py --root-dir <root_dir>`,
  **without** `--mutation-rate`: with the scalar supplied a discovery miss
  reports as `scalar` rather than `MISSING`, which is the case you want to see.
  Since v1.10 also stamps `mu_source` into every ARG and reports it in the
  pipeline summary, a completed run can be audited after the fact instead.
