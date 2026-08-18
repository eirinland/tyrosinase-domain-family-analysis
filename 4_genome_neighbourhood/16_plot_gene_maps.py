#!/usr/bin/env python3
"""
plot_gene_maps.py - render genome-neighbourhood (GNA) gene-arrow figures.

Produces the manuscript supplementary genome-neighbourhood figure and the
companion characterised-BGC reference figure as vector PDFs, from a curated
panel spec (TSV). Each panel is a representative gene neighbourhood for one
active-site substitution group; the gene arrows, categories and marker-gene
labels are taken from the spec, which was curated from - and validated
against - the pipeline outputs:

    4_genome_neighbourhood/groups/<GROUP>/neighbourhoods.tsv          (non-canonical panels)
    4_genome_neighbourhood_canonical/bacterial/summary_by_group.tsv  (canonical bacterial panels)
    supplementary_tables/table_gna_highlights.tsv                    (conserved-locus annotations)

Why a spec rather than a direct render of neighbourhoods.tsv: the raw
neighbourhoods are dominated by "hypothetical protein" flanks, and each
published panel is a hand-chosen representative with shortened marker labels.
The spec captures that curation so the figure is reproducible and editable;
every labelled marker gene in the spec traces to a real row in the tables
above.

Dependencies: matplotlib only (no dna_features_viewer / BioPython needed).

Usage:
    python plot_gene_maps.py --spec gna_gene_maps_spec.tsv --out gna_gene_maps.pdf
    python plot_gene_maps.py --spec bgc_gene_maps_spec.tsv --out bgc_gene_maps.pdf --title-mode bgc

Spec TSV columns:
    panel        panel id (A, B, ...); one row of arrows per panel
    species      species name (italicised in the panel header)
    substitution active-site substitution label (e.g. "Gly46>Tyr"); "" for BGC ref panels
    locus        conserved-locus / cluster annotation shown after the substitution
    idx          gene order within the panel (0-based, left to right)
    label        marker-gene label; "" for unlabelled filler genes
    category     one of: ppo, pathway, biosynthetic, regulatory, transport, hypothetical
    strand       "+" (arrow points right) or "-" (points left)
    labelpos     "above", "below" or "none"
"""

import argparse
import csv
from collections import OrderedDict

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Patch
from matplotlib.lines import Line2D

# ---- category -> fill / edge colours (matched to the published figure legend) ----
CATEGORY_COLORS = OrderedDict([
    ("ppo",          ("#A9CCE3", "#2E75B6")),   # light blue, blue edge - the query PPO gene
    ("pathway",      ("#A9DFBF", "#4F9E6A")),   # green  - pathway-associated
    ("biosynthetic", ("#F5CBA7", "#D98A3D")),   # orange - biosynthetic
    ("regulatory",   ("#7F8C8D", "#566573")),   # dark grey - regulatory
    ("transport",    ("#D7BDE2", "#9B6FB5")),   # purple - transport
    ("hypothetical", ("#D5D8DC", "#AEB6BF")),   # light grey - hypothetical/other
])
CATEGORY_LABELS = OrderedDict([
    ("ppo", "PPO"),
    ("pathway", "Pathway-associated"),
    ("biosynthetic", "Biosynthetic"),
    ("regulatory", "Regulatory"),
    ("transport", "Transport"),
    ("hypothetical", "Hypothetical/other"),
])

# ---- geometry (data units; 1 gene slot = GENE_W + GAP) ----
GENE_W = 1.0      # total gene length (block + point)
GENE_H = 0.50     # gene block height
HEAD_W = 0.34     # length of the pointed tip (portion of GENE_W)
GAP = 0.14        # gap between consecutive genes
ROW_DY = 3.1      # vertical spacing between panels


def _read_spec(path):
    panels = OrderedDict()
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            p = row["panel"]
            panels.setdefault(p, {"section": row.get("section", ""),
                                  "species": row["species"],
                                  "substitution": row["substitution"],
                                  "locus": row["locus"],
                                  "genes": []})
            panels[p]["genes"].append(row)
    for p in panels:
        panels[p]["genes"].sort(key=lambda r: int(r["idx"]))
    return panels


def _draw_gene(ax, x0, y, gene):
    """Draw one pentagon ("home-plate") gene arrow, left edge at x0, centred on y.

    Forward genes point right, reverse genes point left; the whole gene is a
    single full-height block with one pointed end (standard gene-map style).
    """
    fill, edge = CATEGORY_COLORS.get(gene["category"], CATEGORY_COLORS["hypothetical"])
    forward = gene["strand"] != "-"
    h = GENE_H / 2.0
    x1 = x0 + GENE_W
    if forward:
        # flat back on the left, point on the right
        pts = [(x0, y - h), (x0, y + h), (x1 - HEAD_W, y + h),
               (x1, y), (x1 - HEAD_W, y - h)]
    else:
        # point on the left, flat back on the right
        pts = [(x1, y - h), (x1, y + h), (x0 + HEAD_W, y + h),
               (x0, y), (x0 + HEAD_W, y - h)]
    poly = Polygon(pts, closed=True, facecolor=fill, edgecolor=edge,
                   linewidth=0.6, joinstyle="miter")
    ax.add_patch(poly)
    lbl = gene["label"].strip()
    if lbl and gene["labelpos"] != "none":
        dy = GENE_H * 0.9 + 0.14 if gene["labelpos"] == "above" else -(GENE_H * 0.9 + 0.14)
        va = "bottom" if gene["labelpos"] == "above" else "top"
        ax.text(x0 + GENE_W / 2, y + dy, lbl, ha="center", va=va,
                fontsize=6.0, fontstyle="italic")


def render(spec_path, out_path, title_mode="gna"):
    panels = _read_spec(spec_path)
    n = len(panels)
    max_genes = max(len(p["genes"]) for p in panels.values())
    fig_w = max(7.0, 0.55 * max_genes)
    fig_h = 0.9 + ROW_DY * n / 3.0 + 0.9
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    xmax = max_genes * (GENE_W + GAP)
    y = 0.0
    current_section = None
    for pid, panel in panels.items():
        genes = panel["genes"]
        # section header (e.g. NON-CANONICAL / CANONICAL BACTERIAL) with divider
        sec = panel.get("section", "").strip()
        if sec and sec != current_section:
            sec_y = y + GENE_H + 1.25
            ax.text(0, sec_y, sec, ha="left", va="bottom", fontsize=8.0,
                    fontweight="bold", color="#333333")
            ax.add_line(Line2D([0, xmax], [sec_y - 0.12, sec_y - 0.12],
                               color="#999999", linewidth=0.7))
            current_section = sec
        # panel header (letter + italic species + substitution - locus)
        header_y = y + GENE_H * 1.7 + 0.55
        sub = panel["substitution"].strip()
        locus = panel["locus"].strip()
        ax.text(0, header_y, pid, ha="left", va="bottom", fontsize=8.5, fontweight="bold")
        # species italic, then the rest roman
        ax.text(0.55, header_y, panel["species"], ha="left", va="bottom",
                fontsize=7.5, fontstyle="italic")
        tail = f"   {sub} - {locus}" if sub else f"   {locus}"
        if tail.strip():
            ax.text(0.55 + 0.16 * len(panel["species"]), header_y, tail,
                    ha="left", va="bottom", fontsize=7.0)
        x = 0.0
        for g in genes:
            _draw_gene(ax, x, y, g)
            x += GENE_W + GAP
        y -= ROW_DY

    # legend
    handles = [Patch(facecolor=CATEGORY_COLORS[c][0], edgecolor=CATEGORY_COLORS[c][1],
                     label=CATEGORY_LABELS[c]) for c in CATEGORY_COLORS]
    ax.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.02),
              ncol=len(handles), frameon=False, fontsize=6.5, handlelength=1.2,
              columnspacing=1.1)

    ax.set_xlim(-0.4, max_genes * (GENE_W + GAP) + 0.4)
    ax.set_ylim(y + ROW_DY - 1.6, GENE_H + 1.7)
    ax.axis("off")
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight")
    fig.savefig(out_path.rsplit(".", 1)[0] + ".png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    return out_path


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--spec", required=True, help="panel spec TSV")
    ap.add_argument("--out", required=True, help="output PDF path (PNG written alongside)")
    ap.add_argument("--title-mode", default="gna", choices=["gna", "bgc"])
    args = ap.parse_args()
    out = render(args.spec, args.out, args.title_mode)
    print("wrote", out, "and", out.rsplit(".", 1)[0] + ".png")


if __name__ == "__main__":
    main()
