#!/usr/bin/env python3
"""AF3 vs Metal3D vs PinMyMetal (Cu + FECONI): binary prediction rate, grouped bar chart."""
import csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from pathlib import Path

HERE = Path(__file__).parent
NC_DIR = HERE.parent
TSV = HERE / "af3_vs_m3d_exclusive.tsv"
PMM_SITE_TSV = HERE / "pmm_per_site.tsv"
PMM_FE_TSV = HERE / "pmm_feconi_per_site.tsv"
OUT_PNG = HERE / "af3_vs_m3d_bar.png"
OUT_PDF = HERE / "af3_vs_m3d_bar.pdf"

for fname in ["Helvetica Neue", "Helvetica", "Arial"]:
    matches = fm.findSystemFonts()
    if any(fname.lower().replace(" ", "") in f.lower().replace(" ", "") for f in matches):
        plt.rcParams["font.family"] = "sans-serif"
        plt.rcParams["font.sans-serif"] = [fname]
        break

GROUPS = ["3His", "2His+coord", "2His+inert", "1His", "0His"]
GROUP_LABELS = ["3 His\n(canonical)", "2 His +\ncoord.", "2 His +\nnoncoord.", "1 His", "0 His"]

def lighten(c, f=0.45):
    return tuple(ci + (1 - ci) * f for ci in c)

AF3_FILL = lighten((147/255, 196/255, 214/255))
M3D_FILL = lighten((131/255, 180/255, 131/255))
PMM_FILL = lighten((218/255, 183/255, 225/255))
FE_FILL  = lighten((237/255, 230/255, 160/255))
AF3_EDGE = (40/255, 110/255, 170/255)
M3D_EDGE = (40/255, 105/255, 40/255)
PMM_EDGE = (100/255, 40/255, 120/255)
FE_EDGE  = (130/255, 125/255, 55/255)

data = {g: {"af3_yes": 0, "m3d_yes": 0, "pmm_yes": 0, "fe_yes": 0, "total": 0} for g in GROUPS}
acc_site_group = {}

with open(TSV) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        g = row["chem_group"]
        if g not in data:
            continue
        data[g]["total"] += 1
        acc_site_group[(row["accession"], row["site"])] = g
        if row["af3_plddt"]:
            data[g]["af3_yes"] += 1
        if row["m3d_cu_prob"]:
            data[g]["m3d_yes"] += 1

with open(PMM_SITE_TSV) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        acc, site = row["accession"], row["site"]
        g = acc_site_group.get((acc, site))
        if g is None:
            continue
        data[g]["pmm_yes"] += 1

with open(PMM_FE_TSV) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        acc, site = row["accession"], row["site"]
        g = acc_site_group.get((acc, site))
        if g is None:
            continue
        data[g]["fe_yes"] += 1

x = np.arange(len(GROUPS))
w = 0.19

af3_pct = [100 * data[g]["af3_yes"] / data[g]["total"] for g in GROUPS]
m3d_pct = [100 * data[g]["m3d_yes"] / data[g]["total"] for g in GROUPS]
pmm_pct = [100 * data[g]["pmm_yes"] / data[g]["total"] for g in GROUPS]
fe_pct  = [100 * data[g]["fe_yes"]  / data[g]["total"] for g in GROUPS]
totals = [data[g]["total"] for g in GROUPS]

fig, ax = plt.subplots(figsize=(9, 4.5))
fig.subplots_adjust(left=0.09, right=0.95, bottom=0.18, top=0.88)

ax.bar(x - 1.5*w, af3_pct, w, color=AF3_FILL, edgecolor=AF3_EDGE, lw=1.0, label="AF3", zorder=3)
ax.bar(x - 0.5*w, m3d_pct, w, color=M3D_FILL, edgecolor=M3D_EDGE, lw=1.0, label="AllMetal3D", zorder=3)
ax.bar(x + 0.5*w, pmm_pct, w, color=PMM_FILL, edgecolor=PMM_EDGE, lw=1.0, label="PMM Cu", zorder=3)
ax.bar(x + 1.5*w, fe_pct,  w, color=FE_FILL,  edgecolor=FE_EDGE,  lw=1.0, label="PMM Fe/Co/Ni", zorder=3)

for i, (a, m, p, fe, n) in enumerate(zip(af3_pct, m3d_pct, pmm_pct, fe_pct, totals)):
    ax.text(i - 1.5*w, a + 1.5, f"{a:.0f}%", ha="center", va="bottom", fontsize=6.5, fontweight="bold", color="black")
    ax.text(i - 0.5*w, m + 1.5, f"{m:.0f}%", ha="center", va="bottom", fontsize=6.5, fontweight="bold", color="black")
    ax.text(i + 0.5*w, p + 1.5, f"{p:.0f}%", ha="center", va="bottom", fontsize=6.5, fontweight="bold", color="black")
    ax.text(i + 1.5*w, fe + 1.5, f"{fe:.0f}%", ha="center", va="bottom", fontsize=6.5, fontweight="bold", color="black")
    ax.text(i, -6, f"n={n}", ha="center", va="top", fontsize=8, color="gray")

ax.set_ylim(-10, 115)
ax.set_xlim(-0.6, len(GROUPS) - 0.4)
ax.set_xticks(x)
ax.set_xticklabels(GROUP_LABELS, fontsize=9)
ax.set_ylabel("Sites with metal predicted (%)", fontsize=10)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.22), ncol=4, fontsize=8.5, frameon=False)
ax.axhline(50, color="gray", ls="--", lw=0.6, alpha=0.4)

fig.patch.set_edgecolor("black")
fig.patch.set_linewidth(1.0)
plt.savefig(OUT_PNG, dpi=300, bbox_inches="tight", edgecolor=fig.get_edgecolor())
plt.savefig(OUT_PDF, bbox_inches="tight", edgecolor=fig.get_edgecolor())
print(f"Saved {OUT_PNG} and {OUT_PDF}")
for g in GROUPS:
    print(f"  {g}: AF3={data[g]['af3_yes']}, M3D={data[g]['m3d_yes']}, PMM={data[g]['pmm_yes']}, FE={data[g]['fe_yes']}, total={data[g]['total']}")
