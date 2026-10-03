import time, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, "scripts")
from argtest_common import load_ts, load_remove_intervals, name_to_nodes_map

def t():
    return time.time()

ts = load_ts(Path("admix_tested/step4_trimmed_regions/combined.4/2.tsz"))
print(f"loaded: edges={ts.num_edges} nodes={ts.num_nodes} trees={ts.num_trees} muts={ts.num_mutations}", flush=True)

# (A) one plain simplify on the CLEAN table
tabs = ts.dump_tables()
t0 = t(); tabs.sort(); t1 = t(); tabs.simplify(); t2 = t()
print(f"CLEAN: sort={t1-t0:.1f}s  simplify={t2-t1:.1f}s", flush=True)

# (B) one OLD-style single-interval removal (split at 2 positions, drop a few leaf edges, simplify)
def old_remove_ancestry(ts, samples, left, right):
    def split_edges_at(tables, position):
        for i, edge in enumerate(tables.edges):
            if edge.left < position < edge.right:
                tables.edges[i] = edge.replace(right=position)
                tables.edges.append(edge.replace(left=position))
    tables = ts.dump_tables()
    ta = t(); split_edges_at(tables, left); split_edges_at(tables, right); tb = t()
    drop = np.logical_and.reduce([np.isin(tables.edges.child, samples),
                                  tables.edges.left >= left, tables.edges.right <= right])
    tables.edges.keep_rows(~drop)
    tc = t(); tables.sort(); td = t(); tables.edges.drop_metadata(); tables.simplify(); te = t()
    print(f"OLD one interval: 2x split_edges_at(python loop)={tb-ta:.1f}s  sort={td-tc:.1f}s  simplify={te-td:.1f}s", flush=True)
    return tables.tree_sequence()

n2n = name_to_nodes_map(ts, suffix_to_strip="_anchorwave")
samp = n2n.get("ind0", [])
old_remove_ancestry(ts, samp, 1_000_000.0, 1_050_000.0)
print("DONE", flush=True)
