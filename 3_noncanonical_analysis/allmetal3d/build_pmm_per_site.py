#!/usr/bin/env python3
"""Rebuild pmm_per_site.tsv using mapped residue numbers at all 6 canonical positions.

For each PMM Cu prediction, matches its ligand residue numbers against the
6 canonical-position resnums from noncanonical_analysis.tsv (not just His).
Assigns to CuA or CuB based on ligand overlap (±1 residue tolerance).
"""
import csv
from pathlib import Path

HERE = Path(__file__).parent
NC_DIR = HERE.parent
NC_TSV = NC_DIR / "noncanonical_analysis.tsv"
PMM_RAW = NC_DIR / "pinmymetal" / "pmm_nc_results.tsv"
OUT = HERE / "pmm_per_site.tsv"
OUT_FECONI = HERE / "pmm_feconi_per_site.tsv"

CUA_COLS = ["CuA_His1_resnum", "CuA_His2_resnum", "CuA_His3_resnum"]
CUB_COLS = ["CuB_His1_resnum", "CuB_His2_resnum", "CuB_His3_resnum"]

# Load canonical-position resnums per structure
canon_pos = {}
with open(NC_TSV) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        acc = row["accession"]
        cua = {int(row[c]) for c in CUA_COLS if row.get(c, "").strip()}
        cub = {int(row[c]) for c in CUB_COLS if row.get(c, "").strip()}
        canon_pos[acc] = (cua, cub)

def overlap(ligand_resnums, canon_resnums, tol=1):
    return sum(1 for lr in ligand_resnums
               if any(abs(lr - cr) <= tol for cr in canon_resnums))

stats = {"offsite": 0, "no_canon": 0}

def assign_site(lig_resnums, cua_pos, cub_pos):
    oa = overlap(lig_resnums, cua_pos)
    ob = overlap(lig_resnums, cub_pos)
    if oa == 0 and ob == 0:
        return None
    if oa > ob:
        return "CuA"
    if ob > oa:
        return "CuB"
    return "CuA" if min(lig_resnums) < (min(cua_pos | cub_pos) + max(cua_pos | cub_pos)) / 2 else "CuB"

# Collect all on-site predictions keyed by (acc, site)
# Each key gets the best Cu and best FECONI; then we pick one winner per key
all_hits = {}  # (acc, site) -> {"copper": prob, "feconi": prob}

with open(PMM_RAW) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        acc = row["accession"]
        if acc not in canon_pos:
            continue
        cua_pos, cub_pos = canon_pos[acc]
        if not cua_pos and not cub_pos:
            stats["no_canon"] += 1
            continue

        probs_raw = row.get("all_probs", "")
        ligands_raw = row.get("ligands", "")
        if not probs_raw or not ligands_raw:
            continue

        site_probs = probs_raw.split(";")
        site_ligands = ligands_raw.split(" | ")

        for i, prob_entry in enumerate(site_probs):
            if i >= len(site_ligands):
                break
            metal, prob = prob_entry.rsplit(":", 1)
            metal = metal.strip()
            prob = float(prob)

            lig_resnums = []
            for lig in site_ligands[i].split(";"):
                parts = lig.strip().split(",")
                if len(parts) >= 1 and parts[0].strip().isdigit():
                    lig_resnums.append(int(parts[0].strip()))
            if not lig_resnums:
                continue

            site = assign_site(lig_resnums, cua_pos, cub_pos)
            if site is None:
                stats["offsite"] += 1
                continue

            key = (acc, site)
            if key not in all_hits:
                all_hits[key] = {}

            if metal.lower() == "copper":
                all_hits[key]["copper"] = max(all_hits[key].get("copper", 0), prob)
            elif metal.upper() == "FECONI":
                all_hits[key]["feconi"] = max(all_hits[key].get("feconi", 0), prob)

# Deduplicate: each (acc, site) is either Cu or FECONI, not both.
# Cu wins if present (PMM did find copper); FECONI only where Cu is absent.
cu_results = []
feconi_results = []
for (acc, site), metals in all_hits.items():
    if "copper" in metals:
        cu_results.append({"accession": acc, "site": site, "pmm_cu_prob": f"{metals['copper']:.4f}"})
    elif "feconi" in metals:
        feconi_results.append({"accession": acc, "site": site, "pmm_feconi_prob": f"{metals['feconi']:.4f}"})

with open(OUT, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["accession", "site", "pmm_cu_prob"], delimiter="\t")
    w.writeheader()
    w.writerows(cu_results)

with open(OUT_FECONI, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["accession", "site", "pmm_feconi_prob"], delimiter="\t")
    w.writeheader()
    w.writerows(feconi_results)

cua_cu = sum(1 for r in cu_results if r["site"] == "CuA")
cub_cu = sum(1 for r in cu_results if r["site"] == "CuB")
cua_fe = sum(1 for r in feconi_results if r["site"] == "CuA")
cub_fe = sum(1 for r in feconi_results if r["site"] == "CuB")
n_both = sum(1 for m in all_hits.values() if "copper" in m and "feconi" in m)
print(f"Cu:     {len(cu_results)} entries (CuA={cua_cu}, CuB={cub_cu})")
print(f"FECONI: {len(feconi_results)} entries (CuA={cua_fe}, CuB={cub_fe})")
print(f"Both Cu+FECONI at same site (Cu wins): {n_both}")
print(f"Off-site: {stats['offsite']}, no canon resnums: {stats['no_canon']}")
