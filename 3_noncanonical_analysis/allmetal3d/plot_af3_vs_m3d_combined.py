#!/usr/bin/env python3
"""Combined bar + violin figure: Cu/FECONI prediction rate (top) and confidence (bottom)."""
import csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from matplotlib.patches import Patch
from pathlib import Path

HERE = Path(__file__).parent
NC_DIR = HERE.parent
TSV = HERE / "af3_vs_m3d_exclusive.tsv"
PMM_SITE_TSV = HERE / "pmm_per_site.tsv"
PMM_FE_TSV = HERE / "pmm_feconi_per_site.tsv"
OUT_PNG = HERE / "af3_vs_m3d_combined.png"
OUT_PDF = HERE / "af3_vs_m3d_combined.pdf"

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

# ---------- Load data ----------
bar_data = {g: {"af3_yes": 0, "m3d_yes": 0, "pmm_yes": 0, "fe_yes": 0, "total": 0} for g in GROUPS}
vio_data = {g: {"af3": [], "m3d": [], "pmm": [], "feconi": [], "n_total": 0} for g in GROUPS}
acc_site_group = {}

with open(TSV) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        g = row["chem_group"]
        if g not in bar_data:
            continue
        bar_data[g]["total"] += 1
        vio_data[g]["n_total"] += 1
        acc_site_group[(row["accession"], row["site"])] = g
        if row["af3_plddt"]:
            bar_data[g]["af3_yes"] += 1
            vio_data[g]["af3"].append(float(row["af3_plddt"]) / 100)
        if row["m3d_cu_prob"]:
            bar_data[g]["m3d_yes"] += 1
            vio_data[g]["m3d"].append(float(row["m3d_cu_prob"]))

with open(PMM_SITE_TSV) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        acc, site = row["accession"], row["site"]
        g = acc_site_group.get((acc, site))
        if g is None:
            continue
        prob = float(row["pmm_cu_prob"])
        bar_data[g]["pmm_yes"] += 1
        if prob > 0:
            vio_data[g]["pmm"].append(prob)

with open(PMM_FE_TSV) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        acc, site = row["accession"], row["site"]
        g = acc_site_group.get((acc, site))
        if g is None:
            continue
        prob = float(row["pmm_feconi_prob"])
        bar_data[g]["fe_yes"] += 1
        if prob > 0:
            vio_data[g]["feconi"].append(prob)

# ---------- Figure ----------
fig, (ax_bar, ax_vio) = plt.subplots(2, 1, figsize=(9, 7.5),
    gridspec_kw={"height_ratios": [1, 1.1], "hspace": 0.08})

x = np.arange(len(GROUPS))
bw = 0.19
vw = 0.20
jitter_w = 0.06
rng = np.random.RandomState(42)

METHODS = [
    ("af3",    -0.36, AF3_FILL, AF3_EDGE, "AF3 (pLDDT/100)"),
    ("m3d",    -0.12, M3D_FILL, M3D_EDGE, "AllMetal3D"),
    ("pmm",     0.12, PMM_FILL, PMM_EDGE, "PMM Cu"),
    ("feconi",  0.36, FE_FILL,  FE_EDGE,  "PMM Fe/Co/Ni"),
]
# ----- Top: bar chart -----
af3_pct = [100 * bar_data[g]["af3_yes"] / bar_data[g]["total"] for g in GROUPS]
m3d_pct = [100 * bar_data[g]["m3d_yes"] / bar_data[g]["total"] for g in GROUPS]
pmm_pct = [100 * bar_data[g]["pmm_yes"] / bar_data[g]["total"] for g in GROUPS]
fe_pct  = [100 * bar_data[g]["fe_yes"]  / bar_data[g]["total"] for g in GROUPS]

ax_bar.bar(x - 1.5*bw, af3_pct, bw, color=AF3_FILL, edgecolor=AF3_EDGE, lw=1.0, zorder=3)
ax_bar.bar(x - 0.5*bw, m3d_pct, bw, color=M3D_FILL, edgecolor=M3D_EDGE, lw=1.0, zorder=3)
ax_bar.bar(x + 0.5*bw, pmm_pct, bw, color=PMM_FILL, edgecolor=PMM_EDGE, lw=1.0, zorder=3)
ax_bar.bar(x + 1.5*bw, fe_pct,  bw, color=FE_FILL,  edgecolor=FE_EDGE,  lw=1.0, zorder=3)

for i, (a, m, p, fe) in enumerate(zip(af3_pct, m3d_pct, pmm_pct, fe_pct)):
    ax_bar.text(i - 1.5*bw, a + 1.5, f"{a:.0f}%", ha="center", va="bottom", fontsize=6, fontweight="bold", color="black")
    ax_bar.text(i - 0.5*bw, m + 1.5, f"{m:.0f}%", ha="center", va="bottom", fontsize=6, fontweight="bold", color="black")
    ax_bar.text(i + 0.5*bw, p + 1.5, f"{p:.0f}%", ha="center", va="bottom", fontsize=6, fontweight="bold", color="black")
    ax_bar.text(i + 1.5*bw, fe + 1.5, f"{fe:.0f}%", ha="center", va="bottom", fontsize=6, fontweight="bold", color="black")

ax_bar.set_ylim(0, 115)
ax_bar.set_xlim(-0.6, len(GROUPS) - 0.4)
ax_bar.set_xticks(x)
ax_bar.set_xticklabels([])
ax_bar.set_ylabel("Sites with metal predicted (%)", fontsize=10)
ax_bar.axhline(50, color="gray", ls="--", lw=0.6, alpha=0.4)
ax_bar.spines["top"].set_visible(False)
ax_bar.spines["right"].set_visible(False)
ax_bar.spines["bottom"].set_visible(False)
ax_bar.tick_params(bottom=False)

# ----- Bottom: violin plot -----
for i, g in enumerate(GROUPS):
    ntot = vio_data[g]["n_total"]
    for key, dx, fill, edge, label in METHODS:
        vals = vio_data[g][key]
        pos = i + dx
        if vals:
            parts = ax_vio.violinplot([vals], positions=[pos], showmedians=False,
                                      showextrema=False, widths=vw)
            for pc in parts["bodies"]:
                pc.set_facecolor(fill)
                pc.set_edgecolor(edge)
                pc.set_linewidth(1.0)
                pc.set_alpha(0.85)
                pc.set_zorder(2)
            jx = rng.uniform(-jitter_w, jitter_w, len(vals))
            ax_vio.scatter(pos + jx, vals, s=3, alpha=0.4, color=edge,
                           edgecolors="none", zorder=4, rasterized=True)
            med = np.median(vals)
            ax_vio.plot([pos - 0.06, pos + 0.06], [med, med], color="black", lw=1.5, zorder=5)

    ax_vio.text(i, -0.08, f"n={ntot}", ha="center", va="top", fontsize=7.5, color="#555555")

ax_vio.set_ylim(-0.14, 1.05)
ax_vio.set_xlim(-0.6, len(GROUPS) - 0.4)
ax_vio.set_xticks(x)
ax_vio.set_xticklabels(GROUP_LABELS, fontsize=9)
ax_vio.set_ylabel("Confidence (pLDDT/100 or probability)", fontsize=10)
ax_vio.axhline(0.5, color="gray", ls="--", lw=0.7, alpha=0.4)
ax_vio.spines["top"].set_visible(False)
ax_vio.spines["right"].set_visible(False)

# ----- Shared legend -----
legend_handles = [Patch(facecolor=fill, edgecolor=edge, linewidth=1.0, alpha=0.85, label=l)
                  for _, _, fill, edge, l in METHODS]
ax_bar.legend(handles=legend_handles, fontsize=8.5, loc="upper right", frameon=False, ncol=2)

ax_bar.text(-0.08, 1.02, "A", transform=ax_bar.transAxes, fontsize=14, fontweight="bold", va="top")
ax_vio.text(-0.08, 1.02, "B", transform=ax_vio.transAxes, fontsize=14, fontweight="bold", va="top")

fig.patch.set_edgecolor("black")
fig.patch.set_linewidth(1.0)

plt.savefig(OUT_PNG, dpi=300, bbox_inches="tight", edgecolor=fig.get_edgecolor())
plt.savefig(OUT_PDF, bbox_inches="tight", edgecolor=fig.get_edgecolor())
print(f"Saved {OUT_PNG} and {OUT_PDF}")
for g in GROUPS:
    d = bar_data[g]
    print(f"  {g}: AF3={d['af3_yes']}, M3D={d['m3d_yes']}, PMM={d['pmm_yes']}, FE={d['fe_yes']}, total={d['total']}")
