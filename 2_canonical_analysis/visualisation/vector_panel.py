"""Side panel: top vectors as a residue table with activity and taxonomy."""

import csv
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches

NODE_CSV = Path(__file__).parent / 'network_nodes.csv'

POSITIONS = ['Gly46','Phe65','Trp68','Glu195','Asn205',
             'Arg209','Val218','Ala221','Phe227','His230','thioether']
DISPLAY =   ['Gly\n46','Phe\n65','Trp\n68','Glu\n195','Asn\n205',
             'Arg\n209','Val\n218','Ala\n221','Phe\n227','His\n230','Cys']

CONSENSUS = {'Gly46':'~','Phe65':'F','Trp68':'W','Glu195':'E','Asn205':'N',
             'Arg209':None,'Val218':None,'Ala221':None,'Phe227':'F','His230':'H',
             'thioether':None}

ACTIVITY_COLORS = {
    'TYR':'#58A0DF','CaOx':'#9BD0F0','AUS':'#B5EAE4','oAPO':'#BEF0BE',
    'oMP':'#7BC486','DHICA ox':'#FFDBBB','DCT':'#F9BE9E','hemocyanin':'#D0B5E4',
}

KINGDOM_COLORS = {
    'Fungi':'#EAD2B9','Animals':'#83B3F4','Plants':'#5A9163',
    'Bacteria':'#A9D8E8','Oomycota':'#BAA3C5',
}

N_TOP = 20


def main():
    with open(NODE_CSV) as f:
        rows = list(csv.DictReader(f))

    top = sorted(rows, key=lambda r: -int(r['count']))[:N_TOP]

    n_rows = len(top)
    n_cols = len(POSITIONS)

    fig, ax = plt.subplots(figsize=(5.5, 0.38 * n_rows + 1.8))

    for i, r in enumerate(top):
        y = n_rows - 1 - i
        for j, pos in enumerate(POSITIONS):
            val = r[pos]
            cons = CONSENSUS[pos]
            if cons is not None and val != cons:
                bg = '#FFE0E0'
            else:
                bg = '#F5F5F5'
            rect = mpatches.FancyBboxPatch((j - 0.42, y - 0.38), 0.84, 0.76,
                                            boxstyle="round,pad=0.04",
                                            facecolor=bg, edgecolor='#CCCCCC', linewidth=0.5)
            ax.add_patch(rect)
            ax.text(j, y, val, ha='center', va='center', fontsize=9,
                    fontfamily='monospace', fontweight='bold' if bg == '#FFE0E0' else 'normal')

    for i, r in enumerate(top):
        y = n_rows - 1 - i
        count = int(r['count'])
        ax.text(-0.8, y, f'{count:,}', ha='right', va='center', fontsize=8, color='#444444')

        act = r['activity'] if r['activity'] != 'none' else ''
        if act:
            acts = act.split('|')
            label = acts[0]
            color = ACTIVITY_COLORS.get(label, '#888888')
            ax.text(n_cols - 0.3, y, label, ha='left', va='center', fontsize=7.5,
                    fontweight='bold', color=color,
                    bbox=dict(boxstyle='round,pad=0.2', facecolor=color, alpha=0.2, edgecolor='none'))

        kingdom = r['top_kingdom']
        kc = KINGDOM_COLORS.get(kingdom, '#CCCCCC')
        ax.plot(-1.1, y, 'o', color=kc, markersize=6, markeredgecolor='#666666', markeredgewidth=0.4)

    ax.set_xlim(-1.5, n_cols + 1.8)
    ax.set_ylim(-0.7, n_rows - 0.3)
    ax.set_xticks(range(n_cols))
    ax.set_xticklabels(DISPLAY, fontsize=9, ha='center')
    ax.xaxis.set_ticks_position('top')
    ax.set_yticks([])
    ax.spines[:].set_visible(False)
    ax.tick_params(length=0)

    ax.text(-0.8, n_rows - 0.05, 'n', ha='right', va='bottom', fontsize=8,
            fontstyle='italic', color='#666666')

    handles = [mpatches.Patch(facecolor=c, edgecolor='#666', label=k, linewidth=0.5)
               for k, c in KINGDOM_COLORS.items()]
    ax.legend(handles=handles, loc='lower left', bbox_to_anchor=(-0.05, -0.02),
              fontsize=7, frameon=False, ncol=3, handletextpad=0.3, columnspacing=1.0,
              title='Dominant kingdom', title_fontsize=7.5)

    fig.subplots_adjust(left=0.12, right=0.88, top=0.90, bottom=0.10)

    out = Path(__file__).parent
    fig.savefig(out / 'vector_panel.pdf', bbox_inches='tight')
    fig.savefig(out / 'vector_panel.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f'Saved: vector_panel.pdf / .png  ({n_rows} vectors)')


if __name__ == '__main__':
    main()
