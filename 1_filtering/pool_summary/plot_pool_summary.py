#!/usr/bin/env python3
"""Summary statistics figure for canonical and non-canonical PPO pools."""

import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["font.sans-serif"] = ["Helvetica Neue", "Helvetica", "Arial"]
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import pandas as pd
import os

from pathlib import Path as _P
BASE = str(_P(__file__).resolve().parents[2])   # repo root
OUTDIR = f"{BASE}/1_filtering/pool_summary"

C_CAN = "#0570B0"
C_NCAN = "#41AB5D"
C_DISC = "#969696"

canon = pd.read_csv(f"{BASE}/canonical_criteria_all_ca.csv")
helix_cov = pd.read_csv(f"{BASE}/1_filtering/helix_coverage.csv")
helix_plddt = pd.read_csv(f"{BASE}/1_filtering/helix_plddt.csv")
overall_plddt = pd.read_csv(f"{BASE}/1_filtering/filtering_edges/overall_plddt.csv")
taxonomy = pd.read_csv(f"{BASE}/taxonomy_lookup.csv")

foldseek = pd.read_csv(
    f"{BASE}/1_filtering/foldseek/foldseek_multiref_best.tsv", sep="\t"
)
foldseek["accession"] = foldseek["query"].str.replace(r"_taxID_.*", "", regex=True)

helix_avg = helix_plddt.copy()
plddt_cols = ["a1_plddt", "a2_plddt", "a3_plddt", "a4_plddt"]
helix_avg["avg_helix_plddt"] = helix_avg[plddt_cols].replace(0, np.nan).mean(axis=1)

df = canon.merge(helix_cov[["accession", "n_helices", "qtmscore"]], on="accession", how="left")
df = df.merge(helix_avg[["accession", "avg_helix_plddt"]], on="accession", how="left")
df = df.merge(overall_plddt[["accession", "overall_plddt"]], on="accession", how="left")
df = df.merge(foldseek[["accession", "alnlen", "rmsd", "lddt"]], on="accession", how="left")
df = df.merge(taxonomy[["accession", "kingdom"]], on="accession", how="left")


final = pd.read_csv(f"{BASE}/1_filtering/final_pools/three_pool_assignment_final.csv")
df = df.merge(final[["accession", "pool"]], on="accession", how="left")
df["pool"] = df["pool"].map({"canonical": "Canonical", "noncanonical": "Non-canonical", "discarded": "Discarded"}).fillna("Discarded")

can = df[df["pool"] == "Canonical"]
ncan = df[df["pool"] == "Non-canonical"]
disc = df[df["pool"] == "Discarded"]

print(f"Canonical: {len(can)}, Non-canonical: {len(ncan)}, Discarded: {len(disc)}")

fig, axes = plt.subplots(3, 2, figsize=(10, 12))
fig.patch.set_facecolor("white")

hist_kw = dict(alpha=0.75, edgecolor="white", linewidth=0.3, density=True)

# --- A: Cu-Cu distance (clipped to 0-20 A for readability) ---
ax = axes[0, 0]
cu_can = can.loc[can["has_2cu"], "cu_dist"]
cu_ncan = ncan.loc[ncan["has_2cu"], "cu_dist"]
bins_cu = np.linspace(0, 20, 60)
ax.hist(cu_can, bins=bins_cu, color=C_CAN, label=f"Canonical (n={len(cu_can):,})", **hist_kw)
ax.hist(cu_ncan, bins=bins_cu, color=C_NCAN, label=f"Non-canonical (n={len(cu_ncan):,})", **hist_kw)
ax.axvline(2.8, color="#333333", ls="--", lw=0.8, alpha=0.6)
ax.axvline(5.5, color="#333333", ls="--", lw=0.8, alpha=0.6)
ax.set_xlim(0, 20)
ax.set_xlabel("Cu–Cu distance (Å)")
ax.set_ylabel("Density")
ax.set_title("A  Cu–Cu distance", loc="left", fontweight="bold")
ax.legend(fontsize=7.5, frameon=False, loc="upper right")

# --- B: Number of coordinating His (percentage) ---
ax = axes[0, 1]

his_vals = list(range(0, 9))
width = 0.35
x_pos = np.arange(len(his_vals))
can_pct = [100 * len(can[can["n_coord_his"] == v]) / len(can) for v in his_vals]
ncan_pct = [100 * len(ncan[ncan["n_coord_his"] == v]) / len(ncan) for v in his_vals]
ax.bar(x_pos - width/2, can_pct, width, color=C_CAN, edgecolor="white", linewidth=0.3, label="Canonical")
ax.bar(x_pos + width/2, ncan_pct, width, color=C_NCAN, edgecolor="white", linewidth=0.3, label="Non-canonical")
ax.set_xticks(x_pos)
ax.set_xticklabels([str(v) for v in his_vals])
ax.set_xlabel("Coordinating His (NE2 ≤ 3.5 Å from Cu)")
ax.set_ylabel("Percentage (%)")
ax.set_title("B  Coordinating histidines", loc="left", fontweight="bold")
ax.legend(fontsize=7.5, frameon=False)

# --- C: Per-His min CA pLDDT ---
ax = axes[1, 0]
bins_plddt = np.linspace(0, 100, 60)
ax.hist(can["min_plddt"].dropna(), bins=bins_plddt, color=C_CAN, label="Canonical", **hist_kw)
ax.hist(ncan["min_plddt"].dropna(), bins=bins_plddt, color=C_NCAN, label="Non-canonical", **hist_kw)
ax.axvline(70, color="#333333", ls="--", lw=0.8, alpha=0.6)
ax.set_xlabel("Min per-His CA pLDDT")
ax.set_ylabel("Density")
ax.set_title("C  Per-His CA pLDDT (minimum)", loc="left", fontweight="bold")
ax.legend(fontsize=7.5, frameon=False)

# --- D: Average helix pLDDT ---
ax = axes[1, 1]
bins_hplddt = np.linspace(60, 100, 50)
ax.hist(can["avg_helix_plddt"].dropna(), bins=bins_hplddt, color=C_CAN, label="Canonical", **hist_kw)
ax.hist(ncan["avg_helix_plddt"].dropna(), bins=bins_hplddt, color=C_NCAN, label="Non-canonical", **hist_kw)
ax.axvline(70, color="#333333", ls="--", lw=0.8, alpha=0.6)
ax.set_xlabel("Average core helix pLDDT")
ax.set_ylabel("Density")
ax.set_title("D  Core helix pLDDT", loc="left", fontweight="bold")
ax.legend(fontsize=7.5, frameon=False)

# --- E: Foldseek qtmscore ---
ax = axes[2, 0]
bins_qtm = np.linspace(0.2, 1.0, 50)
ax.hist(can["qtmscore"].dropna(), bins=bins_qtm, color=C_CAN, label="Canonical", **hist_kw)
ax.hist(ncan["qtmscore"].dropna(), bins=bins_qtm, color=C_NCAN, label="Non-canonical", **hist_kw)
ax.set_xlabel("Foldseek qTMscore (best reference)")
ax.set_ylabel("Density")
ax.set_title("E  Structural similarity to references", loc="left", fontweight="bold")
ax.legend(fontsize=7.5, frameon=False)

# --- F: Taxonomy breakdown (percentage within each pool) ---
ax = axes[2, 1]
kingdoms = ["Fungi", "Animals", "Plants", "Bacteria", "Oomycota", "Other Eukaryota", "Archaea"]
can_with_tax = can["kingdom"].dropna()
ncan_with_tax = ncan["kingdom"].dropna()
can_tax = can_with_tax.value_counts()
ncan_tax = ncan_with_tax.value_counts()
x_pos = np.arange(len(kingdoms))
can_pct_k = [100 * can_tax.get(k, 0) / len(can_with_tax) if len(can_with_tax) > 0 else 0 for k in kingdoms]
ncan_pct_k = [100 * ncan_tax.get(k, 0) / len(ncan_with_tax) if len(ncan_with_tax) > 0 else 0 for k in kingdoms]
ax.bar(x_pos - width/2, can_pct_k, width, color=C_CAN, edgecolor="white", linewidth=0.3,
       label=f"Canonical (n={len(can_with_tax):,})")
ax.bar(x_pos + width/2, ncan_pct_k, width, color=C_NCAN, edgecolor="white", linewidth=0.3,
       label=f"Non-canonical (n={len(ncan_with_tax):,})")
ax.set_xticks(x_pos)
ax.set_xticklabels(kingdoms, rotation=30, ha="right", fontsize=9)
ax.set_xlabel("Kingdom")
ax.set_ylabel("Percentage (%)")
ax.set_title("F  Taxonomic distribution", loc="left", fontweight="bold")
ax.legend(fontsize=7.5, frameon=False)

for a in axes.flat:
    a.spines["top"].set_visible(False)
    a.spines["right"].set_visible(False)
    a.tick_params(labelsize=9)

fig.suptitle(
    f"Canonical (n={len(can):,}) vs Non-canonical (n={len(ncan):,}) PPO pools",
    fontsize=14, fontweight="bold", y=0.99,
)

plt.tight_layout(rect=[0, 0, 1, 0.95])

for fmt in ["pdf", "png"]:
    out = os.path.join(OUTDIR, f"pool_summary.{fmt}")
    fig.savefig(out, dpi=300, facecolor="white", bbox_inches="tight")
    print(f"Saved: {out}")

plt.close()

print("\n=== Pool summary statistics ===")
for label, sub in [("Canonical", can), ("Non-canonical", ncan)]:
    print(f"\n{label} (n={len(sub):,}):")
    cu_sub = sub[sub["has_2cu"]]
    print(f"  Cu-Cu dist: {cu_sub['cu_dist'].median():.2f} A (median), {cu_sub['cu_dist'].mean():.2f} +/- {cu_sub['cu_dist'].std():.2f}")
    print(f"  Coord His:  {sub['n_coord_his'].median():.0f} (median), range {sub['n_coord_his'].min():.0f}-{sub['n_coord_his'].max():.0f}")
    print(f"  Min His CA pLDDT: {sub['min_plddt'].median():.1f} (median)")
    print(f"  Avg helix pLDDT: {sub['avg_helix_plddt'].median():.1f} (median)")
    print(f"  qTMscore: {sub['qtmscore'].median():.3f} (median)")
    tax = sub["kingdom"].value_counts()
    top = ", ".join(f"{k}: {v}" for k, v in tax.head(6).items())
    print(f"  Taxonomy: {top}")
    print(f"  No Cu: {len(sub[~sub['has_2cu']])}, 1 Cu: {len(sub[sub['n_cu'] == 1])}")
