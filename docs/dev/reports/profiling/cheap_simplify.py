import time, sys
from pathlib import Path
sys.path.insert(0, "scripts")
from argtest_common import load_ts

CHEAP = dict(filter_nodes=False, filter_sites=False, filter_individuals=False,
             filter_populations=False, record_provenance=False)

ts = load_ts(Path("admix_tested/step4_trimmed_regions/combined.4/2.tsz"))
L = ts.sequence_length
print(f"full: edges={ts.num_edges} nodes={ts.num_nodes} trees={ts.num_trees} muts={ts.num_mutations}", flush=True)

def timed_simplify(tabs, label, **kw):
    t0 = time.time(); tabs.sort(); t1 = time.time(); tabs.simplify(**kw); t2 = time.time()
    print(f"  {label}: sort={t1-t0:.1f}s simplify={t2-t1:.1f}s -> edges={tabs.edges.num_rows} "
          f"nodes={tabs.nodes.num_rows} sites={tabs.sites.num_rows}", flush=True)

# 40 Mb chunk: default vs cheap
sub = ts.delete_intervals([[40e6, L]], simplify=False).trim()
print(f"\n40Mb chunk: edges={sub.num_edges} trees={sub.num_trees} muts={sub.num_mutations}", flush=True)
timed_simplify(sub.dump_tables(), "DEFAULT")
timed_simplify(sub.dump_tables(), "CHEAP  ", **CHEAP)

# Full chromosome: cheap only (default would be ~2h)
print(f"\nFULL chromosome:", flush=True)
timed_simplify(ts.dump_tables(), "CHEAP  ", **CHEAP)
print("DONE", flush=True)
