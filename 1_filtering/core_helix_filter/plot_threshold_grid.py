#!/usr/bin/env python3
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["font.sans-serif"] = ["Helvetica Neue", "Helvetica", "Arial"]
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os, platform

# BASE is resolved by repo_paths so this script runs from a fresh clone.
# Set PPO_BASE to override. Original hardcoded block kept below for
# provenance.
# if platform.system() == "Darwin":
#     BASE = "/Users/eirinlandsem/Library/Mobile Documents/com~apple~CloudDocs/Proteinkjemi_PhD/Skriving/Thesis/manuscripts:articles/Manuscript_TyrY/New_bioinf/bioinf_redo/Super_reference_pipeline"
# else:
#     BASE = "/cluster/work/projects/nn1003k/eirin/bioinf/Super_reference_pipeline"
import os as _os
import sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))))
from repo_paths import BASE

OUTDIR = f"{BASE}/1_filtering/core_helix_filter"
SWEEP = f"{OUTDIR}/sweep_grid.tsv"

cols = ["d3_range", "d4_range", "TP", "FP", "FN", "TN", "Prec", "Rec", "F1"]
df = pd.read_csv(SWEEP, sep="\t", header=None, names=cols)

d3_vals = sorted(df["d3_range"].unique(), key=lambda s: (float(s.split("-")[0]), float(s.split("-")[1])))
d4_vals = sorted(df["d4_range"].unique(), key=lambda s: (float(s.split("-")[0]), float(s.split("-")[1])))

grid = np.full((len(d4_vals), len(d3_vals)), np.nan)
for _, row in df.iterrows():
    i = d4_vals.index(row["d4_range"])
    j = d3_vals.index(row["d3_range"])
    grid[i, j] = row["F1"]

fig, ax = plt.subplots(figsize=(9, 7.5))
fig.patch.set_facecolor("white")

im = ax.imshow(grid, cmap="YlOrRd_r", aspect="auto", vmin=0.82, vmax=0.93,
               origin="lower")

ax.set_xticks(range(len(d3_vals)))
ax.set_xticklabels(d3_vals, rotation=45, ha="right", fontsize=7.5)
ax.set_yticks(range(len(d4_vals)))
ax.set_yticklabels(d4_vals, fontsize=7.5)

ax.set_xlabel("d$_3$ range (C$\\alpha_{i{\\to}i+3}$, Å)", fontsize=11)
ax.set_ylabel("d$_4$ range (C$\\alpha_{i{\\to}i+4}$, Å)", fontsize=11)

for i in range(len(d4_vals)):
    for j in range(len(d3_vals)):
        val = grid[i, j]
        if not np.isnan(val):
            tp = int(df[(df["d3_range"] == d3_vals[j]) & (df["d4_range"] == d4_vals[i])]["TP"].iloc[0])
            color = "white" if val < 0.87 else "#333333"
            ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=5.5, color=color)

prod_d3 = "4.8-6.4"
prod_d4 = "5.4-7.4"
if prod_d3 in d3_vals and prod_d4 in d4_vals:
    pj, pi = d3_vals.index(prod_d3), d4_vals.index(prod_d4)
    ax.add_patch(plt.Rectangle((pj - 0.5, pi - 0.5), 1, 1, linewidth=2.5,
                                edgecolor="#333333", facecolor="none", linestyle="--"))
    ax.annotate("Original\n(F1=0.826)", (pj, pi), xytext=(pj - 2.5, pi + 2),
                fontsize=8, fontweight="bold", color="#333333",
                arrowprops=dict(arrowstyle="->", color="#333333", lw=1.2),
                ha="center", va="center")

final_d3 = "4.0-6.4"
final_d4 = "4.8-8.2"
if final_d3 in d3_vals and final_d4 in d4_vals:
    fj, fi = d3_vals.index(final_d3), d4_vals.index(final_d4)
    ax.add_patch(plt.Rectangle((fj - 0.5, fi - 0.5), 1, 1, linewidth=2.5,
                                edgecolor="#0570B0", facecolor="none"))
    ax.annotate("Final\n(F1=0.924)", (fj, fi), xytext=(fj + 3, fi - 2),
                fontsize=8, fontweight="bold", color="#0570B0",
                arrowprops=dict(arrowstyle="->", color="#0570B0", lw=1.2),
                ha="center", va="center")

cbar = fig.colorbar(im, ax=ax, shrink=0.75, pad=0.02)
cbar.set_label("F1 score", fontsize=10)
cbar.ax.tick_params(labelsize=8)

ax.set_title("M7 threshold sensitivity (241-structure benchmark)\nFP = 0 across all 400 parameter combinations",
             fontsize=12, fontweight="bold", pad=12)

ax.text(0.02, 0.02, "Benchmark: 162 positive, 79 negative\nPrecision = 1.000 everywhere",
        transform=ax.transAxes, fontsize=8, color="#555555", va="bottom")

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

plt.tight_layout()
for fmt in ["pdf", "png"]:
    out = os.path.join(OUTDIR, f"threshold_sensitivity.{fmt}")
    fig.savefig(out, dpi=300, facecolor="white", bbox_inches="tight")
    print(f"Saved: {out}")
plt.close()
