"""Measure peak RSS of merging one replicate's TRIMMED (step5) chroms.

Mirrors merge_treefiles_by_replicate.merge_group: load all chroms, sequential
concatenate, then dump_tables (the metadata-rebuild path = another full copy).
Reports ru_maxrss peak. Writes nothing. Run on `low` so it doesn't consume the
high-QOS budget the live pipeline is using.
"""
import glob, re, resource, sys
from pathlib import Path
sys.path.insert(0, "scripts")
from argtest_common import load_ts

rep = sys.argv[1] if len(sys.argv) > 1 else "0"
pattern = f"admix_tested/step5_trimmed_samples/combined.*/{rep}.tsz"

def chrom_num(p):
    m = re.search(r"combined\.(\d+)", p)
    return int(m.group(1)) if m else 0

paths = sorted(glob.glob(pattern), key=chrom_num)
def peak():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6  # kB -> GB

print(f"replicate {rep}: {len(paths)} chrom files matched", flush=True)
if not paths:
    sys.exit("no files matched")

tseqs = [load_ts(Path(p)) for p in paths]
tot = sum(ts.tables.nbytes for ts in tseqs) / 1e9
print(f"sum per-chrom table nbytes: {tot:.2f} GB | peak so far: {peak():.2f} GB", flush=True)

merged = tseqs[0]
for ts in tseqs[1:]:
    merged = merged.concatenate(ts)
print(f"after sequential concatenate: edges={merged.num_edges} | peak so far: {peak():.2f} GB", flush=True)

_ = merged.dump_tables()
print(f"after dump_tables rebuild | peak so far: {peak():.2f} GB", flush=True)
print(f"\n>>> PEAK RSS for trimmed merge of replicate {rep}: {peak():.2f} GB <<<", flush=True)
print(f"merged genome: {merged.sequence_length/1e6:.0f} Mbp, {merged.num_edges} edges", flush=True)
