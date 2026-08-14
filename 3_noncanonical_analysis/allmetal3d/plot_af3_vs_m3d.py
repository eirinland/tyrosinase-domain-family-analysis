#!/usr/bin/env python3
"""AF3 vs AllMetal3D vs PinMyMetal (Cu + FECONI): merged violin figure.
Four side-by-side violins per chem_group. AF3 pLDDT normalised to 0–1."""
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
OUT_PNG = HERE / "af3_vs_m3d_figure.png"
OUT_PDF = HERE / "af3_vs_m3d_figure.pdf"

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

data = {g: {"af3": [], "m3d": [], "pmm": [], "feconi": [], "n_total": 0} for g in GROUPS}
acc_site_group = {}

with open(TSV) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        g = row["chem_group"]
        if g not in data:
            continue
        data[g]["n_total"] += 1
        acc_site_group[(row["accession"], row["site"])] = g
        af3 = float(row["af3_plddt"]) / 100 if row["af3_plddt"] else None
        m3d = float(row["m3d_cu_prob"]) if row["m3d_cu_prob"] else None
        if af3 is not None:
            data[g]["af3"].append(af3)
        if m3d is not None:
            data[g]["m3d"].append(m3d)

with open(PMM_SITE_TSV) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        acc, site = row["accession"], row["site"]
        g = acc_site_group.get((acc, site))
        if g is None:
            continue
        prob = float(row["pmm_cu_prob"])
        if prob > 0:
            data[g]["pmm"].append(prob)

with open(PMM_FE_TSV) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        acc, site = row["accession"], row["site"]
        g = acc_site_group.get((acc, site))
        if g is None:
            continue
        prob = float(row["pmm_feconi_prob"])
        if prob > 0:
            data[g]["feconi"].append(prob)

fig, ax = plt.subplots(figsize=(11, 5))
fig.subplots_adjust(left=0.08, right=0.95, bottom=0.18, top=0.90)

vw = 0.20
jitter_w = 0.06
rng = np.random.RandomState(42)

METHODS = [
    ("af3",    -0.36, AF3_FILL, AF3_EDGE, "AF3 (pLDDT/100)"),
    ("m3d",    -0.12, M3D_FILL, M3D_EDGE, "AllMetal3D"),
    ("pmm",     0.12, PMM_FILL, PMM_EDGE, "PMM Cu"),
    ("feconi",  0.36, FE_FILL,  FE_EDGE,  "PMM Fe/Co/Ni"),
]

for i, g in enumerate(GROUPS):
    ntot = data[g]["n_total"]
    for key, dx, fill, edge, label in METHODS:
        vals = data[g][key]
        pos = i + dx

        if vals:
            parts = ax.violinplot([vals], positions=[pos], showmedians=False,
                                  showextrema=False, widths=vw)
            for pc in parts["bodies"]:
                pc.set_facecolor(fill)
                pc.set_edgecolor(edge)
                pc.set_linewidth(1.0)
                pc.set_alpha(0.85)
                pc.set_zorder(2)

            jx = rng.uniform(-jitter_w, jitter_w, len(vals))
            ax.scatter(pos + jx, vals, s=4, alpha=0.4, color=edge,
                       edgecolors="none", zorder=4, rasterized=True)

            med = np.median(vals)
            ax.plot([pos - 0.08, pos + 0.08], [med, med], color="black", lw=1.5, zorder=5)

    ax.text(i, -0.08, f"n={ntot}", ha="center", va="top", fontsize=7.5, color="#555555")

ax.set_ylim(-0.14, 1.05)
ax.set_xlim(-0.6, len(GROUPS) - 0.4)
ax.set_xticks(range(len(GROUPS)))
ax.set_xticklabels(GROUP_LABELS, fontsize=9)
ax.set_ylabel("Confidence (pLDDT/100 or probability)", fontsize=10)
ax.axhline(0.5, color="gray", ls="--", lw=0.7, alpha=0.4)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

from matplotlib.patches import Patch
legend_handles = [Patch(facecolor=fill, edgecolor=edge, linewidth=1.0, alpha=0.85, label=l)
                  for _, _, fill, edge, l in METHODS]
ax.legend(handles=legend_handles, loc="upper right", fontsize=8, framealpha=0.9,
          edgecolor="#cccccc", handlelength=1.2, handleheight=0.9)

fig.patch.set_edgecolor("black")
fig.patch.set_linewidth(1.0)
plt.savefig(OUT_PNG, dpi=300, bbox_inches="tight", edgecolor=fig.get_edgecolor())
plt.savefig(OUT_PDF, bbox_inches="tight", edgecolor=fig.get_edgecolor())
print(f"Saved {OUT_PNG} and {OUT_PDF}")
for g in GROUPS:
    print(f"  {g}: AF3={len(data[g]['af3'])}, M3D={len(data[g]['m3d'])}, PMM={len(data[g]['pmm'])}, FECONI={len(data[g]['feconi'])}, total={data[g]['n_total']}")
