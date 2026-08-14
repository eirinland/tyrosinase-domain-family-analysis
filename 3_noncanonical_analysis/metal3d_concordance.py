#!/usr/bin/env python3
"""Concordance: AllMetal3D metal-site call vs the chemical-plausibility classification,
on the REVISED non-canonical pool (noncanonical_analysis.tsv, 1,687 structures).

Purpose: validate whether AllMetal3D agrees with the residue-based (chemical-plausibility)
classification before adopting AllMetal3D as the primary sorter. Reports a concordance
matrix (rows = chemical-plausibility class, cols = AllMetal3D Cu-site count) + agreement.

AllMetal3D sources (per AF3 Cu atom: closest predicted Cu within DIST_THRESHOLD):
  metal3d/results/combined_all.tsv          original run (main format)
  metal3d/results/metal3d_*.tsv      revised-pipeline additions
                                            (run run_metal3d_new443.sh, then merge_results.py)
  metal3d/results/metal3d_neither_*.tsv,
  metal3d/results/combined_mono.tsv         PmTYR-anchored fallback for the neither/mono sub-runs

Run after the new443 GPU job finishes. Missing-data structures are reported, not silently dropped.
"""
import csv, glob, os
from collections import Counter, defaultdict

BASEDIR = os.path.dirname(os.path.abspath(__file__))
DIST_THRESHOLD = 3.0
PLAUSIBLE_LIGANDS = {"GLU", "ASP", "CYS", "TYR", "MET"}


def tier(h1, h2, h3):
    res = [h1, h2, h3]
    n = res.count("HIS")
    if n == 3:
        return "canonical"
    if n == 2 and [r for r in res if r != "HIS"][0] in PLAUSIBLE_LIGANDS:
        return "plausible"
    return "divergent"


def overall_class(a, b):
    if a != "divergent" and b != "divergent":
        return "binuclear"
    if a != "divergent" or b != "divergent":
        return "mononuclear"
    return "neither"


def has_gap(r):
    return any(r[p] == "---" for p in
              ["CuA_His1", "CuA_His2", "CuA_His3", "CuB_His1", "CuB_His2", "CuB_His3"])


# ---- chemical-plausibility classification (the structures being validated) ----
struct = {}
with open(os.path.join(BASEDIR, "noncanonical_analysis.tsv")) as f:
    for r in csv.DictReader(f, delimiter="\t"):
        if has_gap(r):
            continue
        acc = r["accession"].split("_taxID_")[0]
        cua = tier(r["CuA_His1"], r["CuA_His2"], r["CuA_His3"])
        cub = tier(r["CuB_His1"], r["CuB_His2"], r["CuB_His3"])
        struct[acc] = {"cls": overall_class(cua, cub), "cua_tier": cua, "cub_tier": cub,
                       "cu1_assign": r.get("cu1_assignment", ""),
                       "cu2_assign": r.get("cu2_assignment", "")}

# ---- AllMetal3D, main format (combined_all + new443): keyed (acc -> {af3_cu_index: row}) ----
main = defaultdict(dict)
main_files = [os.path.join(BASEDIR, "metal3d/results/combined_all.tsv")] + \
             sorted(glob.glob(os.path.join(BASEDIR, "metal3d/results/metal3d_*.tsv")))
for fn in main_files:
    if not os.path.exists(fn):
        continue
    with open(fn) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r.get("status") == "ok":
                main[r["accession"].split("_taxID_")[0]][r.get("af3_cu_index", "")] = r

# ---- PmTYR-anchored fallback (neither/mono sub-runs, different column layout) ----
fallback = {}
for fn in sorted(glob.glob(os.path.join(BASEDIR, "metal3d/results/metal3d_neither_*.tsv"))) + \
          [os.path.join(BASEDIR, "metal3d/results/combined_mono.tsv")]:
    if not os.path.exists(fn):
        continue
    with open(fn) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r.get("status") == "ok":
                fallback[r["accession"].split("_taxID_")[0]] = r


def cu_at_site(acc, site):
    """True/False if AllMetal3D predicts Cu within threshold at canonical site CuA/CuB; None if no data."""
    s = struct[acc]
    if acc in main:
        for idx, key in (("1", "cu1_assign"), ("2", "cu2_assign")):
            if s[key] == site and idx in main[acc]:
                d = main[acc][idx].get("closest_cu_dist", "")
                try:
                    return d != "" and float(d) < DIST_THRESHOLD
                except ValueError:
                    return False
    if acc in fallback:
        d = fallback[acc].get(f"{site}_closest_cu_dist", "")
        try:
            return d != "" and float(d) < DIST_THRESHOLD
        except ValueError:
            return False
    return None


# ---- concordance matrix ----
matrix = defaultdict(Counter)
have_data = no_data = 0
M3D = {2: "binuclear", 1: "mononuclear", 0: "neither"}
for acc, s in struct.items():
    a, b = cu_at_site(acc, "CuA"), cu_at_site(acc, "CuB")
    if a is None and b is None:
        no_data += 1
        continue
    have_data += 1
    matrix[s["cls"]][M3D[sum(1 for x in (a, b) if x)]] += 1

print(f"Non-canonical pool (no-gap classifiable): {len(struct)}")
print(f"  with AllMetal3D data: {have_data}   missing: {no_data}"
      + ("   <- run run_metal3d_new443.sh + merge to fill" if no_data else ""))
print(f"  distance threshold: {DIST_THRESHOLD} A\n")
print("CONCORDANCE   rows = chemical-plausibility   cols = AllMetal3D Cu-site count")
print(f"  {'chem \\ m3d':14}{'binuclear':>11}{'mononuclear':>13}{'neither':>9}{'agree':>8}")
order = ["binuclear", "mononuclear", "neither"]
tot = agree = 0
for c in order:
    row = matrix[c]
    n = sum(row.values())
    tot += n
    agree += row[c]
    pct = 100 * row[c] / n if n else 0
    print(f"  {c:14}{row['binuclear']:>11}{row['mononuclear']:>13}{row['neither']:>9}{pct:>7.0f}%")
print(f"\nOverall agreement (diagonal): {agree}/{tot} = {100*agree/tot if tot else 0:.0f}%")
