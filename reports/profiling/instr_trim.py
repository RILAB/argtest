import time, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, "scripts")
from argtest_common import load_ts, load_remove_intervals, name_to_nodes_map, merge_intervals
import trim_samples as T

def stamp(msg, t0):
    print(f"[{time.time()-t0:8.1f}s] {msg}", flush=True)

t0 = time.time()
ts = load_ts(Path("admix_tested/step4_trimmed_regions/combined.4/2.tsz"))
stamp(f"loaded ts edges={ts.num_edges} nodes={ts.num_nodes} trees={ts.num_trees}", t0)

ri = load_remove_intervals(["admix_tested/step3_mutload/combined.4/2.outliers.bed"])
stamp(f"loaded remove intervals names={len(ri)}", t0)

name_to_nodes = name_to_nodes_map(ts, suffix_to_strip="_anchorwave")
seq_len = float(ts.sequence_length)
node_intervals = {}
for name, spans in ri.items():
    nodes = name_to_nodes.get(name, [])
    if not nodes: continue
    pairs = list(zip(spans["starts"], spans["ends"]))
    for node in nodes:
        node_intervals.setdefault(int(node), []).extend(pairs)
stamp(f"built node_intervals nodes={len(node_intervals)}", t0)

tables = ts.dump_tables()
stamp("dump_tables", t0)

edges = tables.edges
e_left, e_right, e_parent, e_child = edges.left, edges.right, edges.parent, edges.child
target_children = np.fromiter(node_intervals.keys(), dtype=e_child.dtype)
is_target = np.isin(e_child, target_children)
out_left=[e_left[~is_target]]; out_right=[e_right[~is_target]]; out_parent=[e_parent[~is_target]]; out_child=[e_child[~is_target]]
stamp(f"partitioned target edges (target={int(is_target.sum())})", t0)

for node, pairs in node_intervals.items():
    sel = np.flatnonzero(e_child == node)
    if sel.size == 0: continue
    rem = np.asarray(merge_intervals(pairs), dtype=np.float64)
    keep_left = np.concatenate(([0.0], rem[:,1])); keep_right = np.concatenate((rem[:,0], [seq_len]))
    ne = keep_left < keep_right; keep_left, keep_right = keep_left[ne], keep_right[ne]
    fL,fR,fP = T._intersect_edges_with_keep(e_left[sel], e_right[sel], e_parent[sel], keep_left, keep_right)
    out_left.append(fL); out_right.append(fR); out_parent.append(fP); out_child.append(np.full(fL.size, node, dtype=e_child.dtype))
stamp("vectorized interval-subtraction over all target children", t0)

edges.set_columns(left=np.concatenate(out_left), right=np.concatenate(out_right),
                  parent=np.concatenate(out_parent).astype(e_parent.dtype), child=np.concatenate(out_child).astype(e_child.dtype))
stamp(f"set_columns new_edges={edges.num_rows}", t0)
tables.sort(); stamp("sort", t0)
tables.simplify(); stamp(f"simplify edges_now={tables.edges.num_rows}", t0)
tables.edges.squash(); stamp("squash", t0)
tables.build_index(); stamp("build_index", t0)
tables.compute_mutation_parents(); stamp("compute_mutation_parents", t0)
out = tables.tree_sequence()
stamp("tree_sequence()", t0)
import tszip
tszip.compress(out, "/tmp/argtest-prof/rep2_instr.tsz"); stamp("tszip.compress (dump)", t0)
print("DONE total", round(time.time()-t0,1), "s")
