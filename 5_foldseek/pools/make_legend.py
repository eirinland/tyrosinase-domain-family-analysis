import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
from matplotlib.patches import FancyBboxPatch, Circle, Wedge

CLADE_COLORS = {
    'a': '#D9A2EB', 'b': '#C9A2EB', 'c': '#B9A2EB', 'd': '#AAA2EB',
    'e': '#A2AAEB', 'f': '#A2B9EB', 'g': '#A2C9EB', 'h': '#A2D9EB',
    'i': '#A5E8B6', 'j': '#B1EAA3', 'k': '#D1EBA2', 'l': '#EDE6A0',
}
TAX_COLORS = {
    'Fungi': '#EAD2B9',
    'Metazoa': '#83B3F4',
    'Bacteria': '#A9D8E8',
    'Viridiplantae': '#5A9163',
    'Other Eukaryota': '#BAA3C5',
    'Unknown': '#D9D9D9',
}
BLUE = '#0570B0'
GREEN = '#41AB5D'

fig = plt.figure(figsize=(7.0, 4.8))
col_centers = [0.155, 0.48, 0.82]

for cx, title in zip(col_centers, ['Thioether fraction', 'Taxonomy', 'Phylogenetic clade']):
    fig.text(cx, 0.97, title, ha='center', va='top', fontsize=10, fontweight='bold')

# Use same node radius in all three; donut extends beyond
r_node = 1.0
lim = 1.4  # accommodate donut overshoot
circle_y = 0.76

for i, cx in enumerate(col_centers):
    ax = fig.add_axes([cx - 0.09, circle_y - 0.09, 0.18, 0.18])
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect('equal')
    ax.axis('off')

    if i == 0:
        ax.add_patch(Wedge((0, 0), r_node, 90, 270, fc=BLUE, ec='none'))
        ax.add_patch(Wedge((0, 0), r_node, 270, 90, fc=GREEN, ec='none'))
        ax.add_patch(Circle((0, 0), r_node, fc='none', ec='#333333', lw=1.2))

    elif i == 1:
        ax.add_patch(Circle((0, 0), r_node, fc='#DDDDDD', ec='#333333', lw=1.2))
        r_inner = r_node * 0.88
        r_outer = r_node * 1.18
        tax_fracs = [0.30, 0.25, 0.12, 0.15, 0.10, 0.08]
        tax_cols = list(TAX_COLORS.values())
        start = 90
        for frac, col in zip(tax_fracs, tax_cols):
            angle = frac * 360
            ax.add_patch(Wedge((0, 0), r_outer, start, start + angle,
                               width=r_outer - r_inner, fc=col, ec='white', lw=0.8))
            start += angle

    elif i == 2:
        ax.add_patch(Circle((0, 0), r_node, fc='#DDDDDD', ec='#333333', lw=1.2))
        r_pie = r_node * 0.92
        pie_fracs = [0.12, 0.09, 0.07, 0.08, 0.10, 0.07, 0.08, 0.07, 0.09, 0.08, 0.07, 0.08]
        pie_cols = list(CLADE_COLORS.values())
        start = 90
        for frac, col in zip(pie_fracs, pie_cols):
            angle = frac * 360
            ax.add_patch(Wedge((0, 0), r_pie, start, start + angle, fc=col, ec='white', lw=0.5))
            start += angle

# --- Gradient bar ---
ax1g = fig.add_axes([col_centers[0] - 0.02, 0.08, 0.04, 0.52])
blue_rgb = np.array(mcolors.hex2color(BLUE))
green_rgb = np.array(mcolors.hex2color(GREEN))
gradient = np.array([blue_rgb * (1 - t) + green_rgb * t for t in np.linspace(0, 1, 256)])
ax1g.imshow(gradient.reshape(256, 1, 3), aspect='auto', origin='lower', extent=[0, 1, 0, 1])
ax1g.set_xticks([])
ax1g.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
ax1g.set_yticklabels(['0', '0.25', '0.5', '0.75', '1.0'], fontsize=7.5)
ax1g.yaxis.tick_right()
ax1g.tick_params(axis='y', length=2, pad=2)
for spine in ax1g.spines.values():
    spine.set_linewidth(0.5)

# --- Taxonomy key ---
key_top = 0.58
key_step = 0.082
box_w = 0.014
box_h = 0.018

for i, (name, col) in enumerate(TAX_COLORS.items()):
    y = key_top - i * key_step
    ax_box = fig.add_axes([col_centers[1] - 0.07, y - box_h/2, box_w, box_h])
    ax_box.set_xlim(0, 1); ax_box.set_ylim(0, 1); ax_box.axis('off')
    ax_box.add_patch(FancyBboxPatch((0.05, 0.05), 0.9, 0.9, boxstyle="round,pad=0.02",
                                     fc=col, ec='#888888', lw=0.5))
    fig.text(col_centers[1] - 0.07 + box_w + 0.015, y, name,
             va='center', ha='left', fontsize=8)

# --- Clade key ---
for i in range(12):
    col_idx = i // 6
    row_idx = i % 6
    x_base = col_centers[2] - 0.1 + col_idx * 0.11
    y = key_top - row_idx * key_step
    ax_box = fig.add_axes([x_base, y - box_h/2, box_w, box_h])
    ax_box.set_xlim(0, 1); ax_box.set_ylim(0, 1); ax_box.axis('off')
    ax_box.add_patch(FancyBboxPatch((0.05, 0.05), 0.9, 0.9, boxstyle="round,pad=0.02",
                                     fc=list(CLADE_COLORS.values())[i], ec='#888888', lw=0.5))
    fig.text(x_base + box_w + 0.01, y, list(CLADE_COLORS.keys())[i],
             va='center', ha='left', fontsize=8)

for ext in ['pdf', 'svg']:
    fig.savefig(f'/tmp/network_legend.{ext}', bbox_inches='tight', dpi=300, transparent=True)
print('Saved')
