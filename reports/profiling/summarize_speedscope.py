#!/usr/bin/env python
"""Summarize a py-spy speedscope profile: top frames by self and total time.

self  = samples where the frame is the leaf (top of stack) -> where CPU actually is
total = samples where the frame appears anywhere -> inclusive cost of that call
"""
import json
import sys
from collections import defaultdict

path = sys.argv[1] if len(sys.argv) > 1 else "reports/profiling/trim_rep2.speedscope.json"
with open(path) as fh:
    doc = json.load(fh)

frames = doc["shared"]["frames"]


def label(i):
    f = frames[i]
    name = f.get("name", "?")
    file = f.get("file", "")
    line = f.get("line", "")
    tail = file.split("/")[-1] if file else ""
    return f"{name}  ({tail}:{line})" if tail else name


self_w = defaultdict(float)
total_w = defaultdict(float)
grand = 0.0

for prof in doc["profiles"]:
    samples = prof["samples"]
    weights = prof.get("weights", [1] * len(samples))
    for stack, w in zip(samples, weights):
        if not stack:
            continue
        grand += w
        self_w[stack[-1]] += w          # leaf = top of stack in speedscope
        for fi in set(stack):
            total_w[fi] += w

def show(title, table):
    print(f"\n=== {title} (grand total = {grand:.0f}) ===")
    for fi, w in sorted(table.items(), key=lambda kv: -kv[1])[:20]:
        print(f"{100*w/grand:6.1f}%  {w:10.0f}  {label(fi)}")

show("TOP by SELF time (leaf)", self_w)
show("TOP by TOTAL time (inclusive)", total_w)
