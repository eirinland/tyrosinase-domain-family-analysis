#!/usr/bin/env python3
"""Continuous pLDDT colour key matching the PyMOL coloring used for the structure
figures: a smooth gradient anchored to the AlphaFold confidence bins so each band
reads as its own colour -- <50 orange | 50-70 yellow | 70-90 light blue |
90-100 dark blue (#0570B0 = network-figure outer ring). Pastel, soft blends at the
50/70/90 boundaries, no green. (general pymol scripts/change_color_gradient_all_structures_pLDDT.py).
Writes plddt_legend.pdf (+ .png)."""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

plt.rcParams["font.family"] = "Helvetica Neue"
plt.rcParams["pdf.fonttype"] = 42  # keep text editable in Illustrator

O = (0.957, 0.694, 0.514)    # #F4B183 pastel orange
Y = (0.984, 0.890, 0.604)    # #FBE39A pastel yellow
L = (0.651, 0.808, 0.890)    # #A6CEE3 pastel light blue
D = (0.0196, 0.4392, 0.6902) # #0570B0 dark blue (network ring)
# (position on 0-1 = pLDDT/100, colour): solid bins, soft blends at 50/70/90
ANCHORS = [(0.00, O), (0.45, O), (0.55, Y), (0.65, Y),
           (0.75, L), (0.85, L), (0.95, D), (1.00, D)]
cmap = LinearSegmentedColormap.from_list("plddt", ANCHORS, N=256)
grad = np.linspace(0, 1, 256).reshape(-1, 1)

fig = plt.figure(figsize=(1.05, 2.6))
ax = fig.add_axes([0.08, 0.06, 0.36, 0.82])
ax.imshow(grad, origin="lower", extent=[0, 1, 0, 100], aspect="auto", cmap=cmap)
ax.set_xticks([])
ax.set_xlim(0, 1)
ax.set_ylim(0, 100)
ax.yaxis.set_ticks_position("right")
ax.set_yticks([0, 50, 70, 90, 100])  # the AlphaFold bin boundaries
ax.tick_params(axis="y", length=3, width=0.6, labelsize=8, pad=2)
for s in ax.spines.values():
    s.set_linewidth(0.6)
ax.set_title("pLDDT", fontsize=9, pad=6)

here = os.path.dirname(os.path.abspath(__file__))
out = os.path.join(here, "plddt_legend.pdf")
fig.savefig(out, transparent=True)
fig.savefig(out.replace(".pdf", ".png"), dpi=300, transparent=True)
print("wrote", out)
