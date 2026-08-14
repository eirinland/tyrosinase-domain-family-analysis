#!/usr/bin/env python3
"""Regenerate nc_patterns_*.tsv from helix_and_gap_filtered_structures.tsv."""
import csv
from collections import Counter, defaultdict
from pathlib import Path

HERE = Path(__file__).parent
SRC = HERE / "helix_and_gap_filtered_structures.tsv"

HIS_COLS = ["CuA_His1", "CuA_His2", "CuA_His3", "CuB_His1", "CuB_His2", "CuB_His3"]
CLS_MAP = {"binuclear": "binuclear", "mononuclear": "mononuclear", "no_cu": "degenerate"}
MIN_N = 2

with open(SRC) as f:
    rows = list(csv.DictReader(f, delimiter="\t"))

groups = defaultdict(list)
for r in rows:
    cls = CLS_MAP[r["classification"]]
    pat = tuple(r[c] for c in HIS_COLS)
    groups[(cls, pat)].append(r)

def top_genera(members, top_n=3):
    gc = Counter(r.get("genus", "") for r in members)
    parts = [f"{g}({n})" for g, n in gc.most_common(top_n) if g]
    return "; ".join(parts)

def n_species(members):
    return len({r.get("species", "") for r in members} - {""})

def n_phyla(members):
    return len({r.get("phylum", "") for r in members} - {""})

header_tier = HIS_COLS + ["n", "n_species", "n_phyla", "Top genera"]
header_compact_full = ["Classification"] + header_tier
header_compact_ms = ["Classification", "H1", "H2", "H3", "H4", "H5", "H6", "n", "n_s", "Top genera"]

def write_tier(filename, cls_name, include_gaps=False):
    tier_groups = [(pat, members) for (cls, pat), members in groups.items()
                   if cls == cls_name and len(members) >= MIN_N]
    if not include_gaps:
        tier_groups = [(p, m) for p, m in tier_groups if "---" not in p]
    tier_groups.sort(key=lambda x: -len(x[1]))
    with open(HERE / filename, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(header_tier)
        for pat, members in tier_groups:
            w.writerow(list(pat) + [len(members), n_species(members),
                                     n_phyla(members), top_genera(members)])
    print(f"  {filename}: {len(tier_groups)} patterns")

for cls_name, fname in [("binuclear", "nc_patterns_binuclear.tsv"),
                         ("mononuclear", "nc_patterns_mononuclear.tsv"),
                         ("degenerate", "nc_patterns_degenerate.tsv")]:
    write_tier(fname, cls_name, include_gaps=True)

all_groups = [(cls, pat, members) for (cls, pat), members in groups.items()
              if len(members) >= MIN_N]
all_groups.sort(key=lambda x: ({"binuclear": 0, "mononuclear": 1, "degenerate": 2}[x[0]], -len(x[2])))

with open(HERE / "nc_patterns_compact.tsv", "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(header_compact_ms)
    for cls, pat, members in all_groups:
        w.writerow([cls] + list(pat) + [len(members), n_species(members),
                                         top_genera(members)])
print(f"  nc_patterns_compact.tsv: {len(all_groups)} patterns")

no_gap = [(c, p, m) for c, p, m in all_groups if "---" not in p]
with open(HERE / "nc_patterns_compact_nogap.tsv", "w", newline="") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(header_compact_ms)
    for cls, pat, members in no_gap:
        w.writerow([cls] + list(pat) + [len(members), n_species(members),
                                         top_genera(members)])
print(f"  nc_patterns_compact_nogap.tsv: {len(no_gap)} patterns (for Table S7)")

print(f"\nTotal structures: {len(rows)}")
print(f"Classification: {Counter(CLS_MAP[r['classification']] for r in rows)}")
