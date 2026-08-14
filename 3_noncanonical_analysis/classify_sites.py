#!/usr/bin/env python3
"""Classify noncanonical structures by chemical plausibility of Cu coordination.

v2 (2026-06-25): for 2-His sites with a coordinatable substitution
(GLU/ASP/CYS/TYR/MET), the substituting residue's nearest coordinating atom
must be within 4.0 A of the canonical Cu position (helix-anchored Kabsch+ICP,
B2ZB02 frame). If *_coord_dist columns are absent (pre-v2 TSV), the check is
skipped and behaviour is identical to v1."""

import csv
from collections import Counter, defaultdict

PLAUSIBLE_LIGANDS = {"GLU", "ASP", "CYS", "TYR", "MET"}
COORD_DIST_CUTOFF = 4.0

def get_coord_dists(row, site):
    dists = []
    for i in range(1, 4):
        v = row.get(f"{site}_His{i}_coord_dist", "")
        dists.append(float(v) if v else None)
    return dists

def site_tier(his1, his2, his3, coord_dists=None):
    residues = [his1, his2, his3]
    n_his = residues.count("HIS")
    if n_his == 3:
        return "canonical"
    if n_his == 2:
        non_his_idx = next(i for i, r in enumerate(residues) if r != "HIS")
        non_his = residues[non_his_idx]
        if non_his in PLAUSIBLE_LIGANDS:
            if coord_dists and coord_dists[non_his_idx] is not None:
                if coord_dists[non_his_idx] > COORD_DIST_CUTOFF:
                    return "divergent"
            return "plausible"
    return "divergent"

def overall_class(cua_tier, cub_tier):
    if cua_tier == "canonical" and cub_tier == "canonical":
        return "Canonical binuclear"
    if cua_tier != "divergent" and cub_tier != "divergent":
        return "Plausible binuclear"
    if cua_tier != "divergent" or cub_tier != "divergent":
        return "Plausible mononuclear"
    return "Divergent"

GAP_POSITIONS = ["CuA_His1", "CuA_His2", "CuA_His3", "CuB_His1", "CuB_His2", "CuB_His3"]

def gap_signature(row):
    return "+".join(p for p in GAP_POSITIONS if row[p] == "---")

def category(row):
    a1, a2 = row["cu1_assignment"], row["cu2_assignment"]
    assignments = {a1, a2}
    if "CuA" in assignments and "CuB" in assignments:
        return "binuclear"
    if ("CuA" in assignments or "CuB" in assignments) and "neither" in assignments:
        return "mononuclear"
    if a1 == "neither" and a2 == "neither":
        return "neither"
    return "other"

with open("noncanonical_analysis.tsv") as f:
    reader = csv.DictReader(f, delimiter="\t")
    rows = list(reader)

classifiable = rows
n_gap = sum(1 for r in rows if gap_signature(r))
n_orient_fail = 0
oc = Counter()
for r in classifiable:
    cua_cd = get_coord_dists(r, "CuA")
    cub_cd = get_coord_dists(r, "CuB")
    cua_t = site_tier(r["CuA_His1"], r["CuA_His2"], r["CuA_His3"], cua_cd)
    cub_t = site_tier(r["CuB_His1"], r["CuB_His2"], r["CuB_His3"], cub_cd)
    cua_noorient = site_tier(r["CuA_His1"], r["CuA_His2"], r["CuA_His3"])
    cub_noorient = site_tier(r["CuB_His1"], r["CuB_His2"], r["CuB_His3"])
    if cua_t != cua_noorient or cub_t != cub_noorient:
        n_orient_fail += 1
    oc[overall_class(cua_t, cub_t)] += 1

binu = oc["Canonical binuclear"] + oc["Plausible binuclear"]
print("=" * 60)
print(f" CHEMICAL-PLAUSIBILITY CLASSIFICATION (primary)")
print(f" {len(classifiable)} classified ({n_gap} carry an unmappable His position, scored divergent)")
print(f" {n_orient_fail} downgraded by coord-atom orientation check (>{COORD_DIST_CUTOFF} A from Cu)")
print("=" * 60)
print(f"  BINUCLEAR  (both sites canonical or plausible): {binu}")
print(f"       Canonical binuclear (both 3-His):  {oc['Canonical binuclear']}")
print(f"       Plausible binuclear:               {oc['Plausible binuclear']}")
print(f"  MONONUCLEAR (one site canonical/plausible):     {oc['Plausible mononuclear']}")
print(f"  NEITHER     (both sites divergent):             {oc['Divergent']}")

CLASS_MAP = {"Canonical binuclear": "binuclear", "Plausible binuclear": "binuclear",
             "Plausible mononuclear": "mononuclear", "Divergent": "no_cu"}
out_fields = list(rows[0].keys()) + ["CuA_tier", "CuB_tier", "gap_signature", "classification"]
out_rows = []
for r in classifiable:
    cua_cd = get_coord_dists(r, "CuA")
    cub_cd = get_coord_dists(r, "CuB")
    cua = site_tier(r["CuA_His1"], r["CuA_His2"], r["CuA_His3"], cua_cd)
    cub = site_tier(r["CuB_His1"], r["CuB_His2"], r["CuB_His3"], cub_cd)
    oc_label = overall_class(cua, cub)
    out = dict(r)
    out["CuA_tier"] = cua
    out["CuB_tier"] = cub
    out["gap_signature"] = gap_signature(r)
    out["classification"] = CLASS_MAP[oc_label]
    out_rows.append(out)

with open("helix_and_gap_filtered_structures.tsv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=out_fields, delimiter="\t")
    w.writeheader()
    w.writerows(out_rows)

cls_counts = Counter(r["classification"] for r in out_rows)
print(f"\nWritten helix_and_gap_filtered_structures.tsv: {len(out_rows)} rows")
for k, v in cls_counts.most_common():
    print(f"  {k}: {v}")

sig_struct = Counter(r["gap_signature"] for r in out_rows if r["gap_signature"])
sig_species = defaultdict(set)
for r in out_rows:
    if r["gap_signature"]:
        sig_species[r["gap_signature"]].add(r.get("species", ""))
print(f"\nGap-signature groups ({sum(sig_struct.values())} structures with >=1 unmapped His):")
print(f"  {'signature':28} {'n':>4} {'species':>8}")
for sig, n in sig_struct.most_common():
    print(f"  {sig:28} {n:>4} {len(sig_species[sig]):>8}")

print("\n--- secondary: breakdown by AF3 Cu-atom assignment (not used for the final call) ---")
for cat_name in ["binuclear", "mononuclear", "neither"]:
    subset = [r for r in rows if category(r) == cat_name]
    print(f"\n{'='*60}")
    print(f" {cat_name.upper()}: {len(subset)} structures")
    print(f"{'='*60}")

    counts = Counter()
    detail = Counter()
    for r in subset:
        cua_cd = get_coord_dists(r, "CuA")
        cub_cd = get_coord_dists(r, "CuB")
        cua = site_tier(r["CuA_His1"], r["CuA_His2"], r["CuA_His3"], cua_cd)
        cub = site_tier(r["CuB_His1"], r["CuB_His2"], r["CuB_His3"], cub_cd)
        oc = overall_class(cua, cub)
        counts[oc] += 1
        detail[(oc, cua, cub)] += 1

    for oc in ["Canonical binuclear", "Plausible binuclear", "Plausible mononuclear", "Divergent"]:
        n = counts.get(oc, 0)
        pct = 100 * n / len(subset) if subset else 0
        print(f"\n  {oc}: {n} ({pct:.1f}%)")
        for (o, ca, cb), c in sorted(detail.items(), key=lambda x: -x[1]):
            if o == oc:
                print(f"    CuA={ca}, CuB={cb}: {c}")
