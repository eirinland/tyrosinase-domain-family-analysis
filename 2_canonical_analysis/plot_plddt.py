#!/usr/bin/env python3
"""Violin plot of per-position pLDDT at the 10 active-site variable positions."""

import csv
from collections import defaultdict
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

POSITION_ORDER = [
    "Gly46", "Phe65", "Trp68", "Glu195", "Asn205",
    "Arg209", "Val218", "Ala221", "Phe227", "His230",
]

data = defaultdict(list)
with open("plddt_per_position.csv") as f:
    for row in csv.DictReader(f):
        data[row["position"]].append(float(row["plddt"]))

positions = [p for p in POSITION_ORDER if p in data]
values = [data[p] for p in positions]

fig, ax = plt.subplots(figsize=(8, 4))
vp = ax.violinplot(values, showmedians=True, showextrema=False)
for body in vp['bodies']:
    body.set_facecolor('#4878CF')
    body.set_alpha(0.7)
vp['cmedians'].set_color('black')

ax.set_xticks(range(1, len(positions) + 1))
ax.set_xticklabels(positions, rotation=45, ha='right', fontsize=9)
ax.set_ylabel('pLDDT')
ax.set_ylim(0, 100)
ax.axhline(70, color='grey', linestyle='--', linewidth=0.8, alpha=0.5)
ax.text(len(positions) + 0.3, 71, 'pLDDT 70', fontsize=7, color='grey', va='bottom')

n = len(values[0])
ax.set_title(f'Per-position pLDDT at active-site variable positions (n = {n:,})')

plt.tight_layout()
plt.savefig('plddt_violin.pdf', dpi=300)
plt.savefig('plddt_violin.png', dpi=300)
print(f"Saved plddt_violin.pdf/png ({len(positions)} positions, {n} structures)")
