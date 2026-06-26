import time, sys
from pathlib import Path
import numpy as np
import tskit
sys.path.insert(0, "scripts")
from argtest_common import load_ts, load_remove_intervals, name_to_nodes_map, merge_intervals
import trim_samples as T

ts = load_ts(Path("admix_tested/step4_trimmed_regions/combined.4/2.tsz"))
L = ts.sequence_length
ri_full = load_remove_intervals(["admix_tested/step3_mutload/combined.4/2.outliers.bed"])
n2n = name_to_nodes_map(ts, suffix_to_strip="_anchorwave")
print(f"full: edges={ts.num_edges} trees={ts.num_trees} muts={ts.num_mutations} L={L:.0f}", flush=True)

def messy_tables(sub_ts, ri):
    """edge-surgery only (no simplify); returns trimmed tables ready to simplify."""
    seq_len = float(sub_ts.sequence_length)
    node_int = {}
    for name, spans in ri.items():
        for nd in n2n.get(name, []):
            pairs = [(s, e) for s, e in zip(spans["starts"], spans["ends"]) if s < seq_len]
            node_int.setdefault(int(nd), []).extend(pairs)
    tb = sub_ts.dump_tables()
    e = tb.edges
    eL, eR, eP, eC = e.left, e.right, e.parent, e.child
    tgt = np.fromiter(node_int.keys(), dtype=eC.dtype)
    m = np.isin(eC, tgt)
    oL=[eL[~m]]; oR=[eR[~m]]; oP=[eP[~m]]; oC=[eC[~m]]
    for nd, pairs in node_int.items():
        sel = np.flatnonzero(eC == nd)
        if sel.size == 0: continue
        merged = merge_intervals(pairs)
        if not merged:
            oL.append(eL[sel]); oR.append(eR[sel]); oP.append(eP[sel]); oC.append(eC[sel]); continue
        rem = np.asarray(merged, float)
        kl = np.concatenate(([0.0], rem[:,1])); kr = np.concatenate((rem[:,0], [seq_len]))
        ne = kl < kr; kl, kr = kl[ne], kr[ne]
        fL,fR,fP = T._intersect_edges_with_keep(eL[sel], eR[sel], eP[sel], kl, kr)
        oL.append(fL); oR.append(fR); oP.append(fP); oC.append(np.full(fL.size, nd, dtype=eC.dtype))
    e.set_columns(left=np.concatenate(oL), right=np.concatenate(oR),
                  parent=np.concatenate(oP).astype(eP.dtype), child=np.concatenate(oC).astype(eC.dtype))
    return tb

for chunk in (5e6, 10e6, 20e6, 40e6):
    sub = ts.delete_intervals([[chunk, L]], simplify=False).trim()
    print(f"\n--- chunk={chunk/1e6:.0f}Mb: edges={sub.num_edges} trees={sub.num_trees} muts={sub.num_mutations} ---", flush=True)

    # clean simplify WITH mutations
    t0=time.time(); tb=sub.dump_tables(); tb.sort(); tb.simplify(); print(f"  clean simplify (with {sub.num_mutations} muts): {time.time()-t0:.1f}s", flush=True)

    # clean simplify WITHOUT mutations
    t0=time.time(); tb=sub.dump_tables(); tb.mutations.clear(); tb.sites.clear(); tb.sort(); tb.simplify(); print(f"  clean simplify (no muts): {time.time()-t0:.1f}s", flush=True)

    # messy (trimmed) simplify WITH mutations
    tb=messy_tables(sub, ri_full); t0=time.time(); tb.sort(); tb.simplify(); print(f"  MESSY simplify (with muts): {time.time()-t0:.1f}s", flush=True)

    # messy simplify WITHOUT mutations
    tb=messy_tables(sub, ri_full); tb.mutations.clear(); tb.sites.clear(); t0=time.time(); tb.sort(); tb.simplify(); print(f"  MESSY simplify (no muts): {time.time()-t0:.1f}s", flush=True)

print("\nDONE", flush=True)
