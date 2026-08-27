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
    ("ppo", "TDF member"),
    ("pathway", "Pathway-associated"),
    ("biosynthetic", "Biosynthetic"),
    ("regulatory", "Regulatory"),
    ("transport", "Transport"),
    ("hypothetical", "Hypothetical/other"),
])

# ---- output size ----
# Target physical width so the figure embeds 1:1 in the manuscript (no Word
# downscaling, which is what shrinks the fonts). 12.7 cm = 5.0 in = text-column
# width. Data units are kept square (fig_w/xspan == fig_h/yspan) so the pentagon
# arrows are never distorted; height follows from the panel count.
TARGET_WIDTH_IN = 5.0          # 12.7 cm

# ---- geometry (data units; 1 gene slot = GENE_W + GAP) ----
GENE_W = 1.0      # total gene length (block + point)
GENE_H = 0.62     # gene block height
HEAD_W = 0.34     # length of the pointed tip (portion of GENE_W)
GAP = 0.14        # gap between consecutive genes
ROW_DY = 4.6      # vertical spacing between panels (room for stacked label tiers)
SECTION_PAD = 2.4  # extra vertical room above a panel that starts a new section

# ---- font sizes (points; true points because the figure embeds 1:1) ----
FS_LABEL = 7.0     # gene marker labels
FS_PANEL = 9.0     # panel letter (A, B, ...)
FS_SPECIES = 8.0   # species name
FS_TAIL = 7.5      # substitution - locus header tail
FS_SECTION = 10.0  # section header (NON-CANONICAL / CANONICAL BACTERIAL)
FS_LEGEND = 8.0    # legend


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


def _text_data_width(ax, fig, s, **kw):
    """Rendered width of string s in data-x units (text drawn off-canvas)."""
    t = ax.text(0, 0, s, **kw)
    fig.canvas.draw()
    bb = t.get_window_extent(renderer=fig.canvas.get_renderer())
    x0 = ax.transData.inverted().transform((bb.x0, bb.y0))[0]
    x1 = ax.transData.inverted().transform((bb.x1, bb.y0))[0]
    t.remove()
    return x1 - x0


def _place_gene_labels(ax, fig, genes, centres, ybase):
    """Place gene labels in stacked tiers above/below the arrow row so wide
    adjacent labels never overlap. Each label gets a short leader line back to
    its gene. Above/below side comes from the spec's labelpos ('above'/'below').
    """
    TIER_DY = GENE_H * 0.95            # vertical step between tiers
    BASE_GAP = GENE_H * 0.55           # first tier offset from arrow edge
    PAD = 0.10                         # min horizontal gap between labels (data units)
    items = []
    for g, cx in zip(genes, centres):
        lbl = g["label"].strip()
        if lbl and g["labelpos"] != "none":
            side = -1 if g["labelpos"] == "below" else +1
            w = _text_data_width(ax, fig, lbl, fontsize=FS_LABEL, fontstyle="italic")
            items.append({"cx": cx, "lbl": lbl, "side": side, "w": w})
    for side in (+1, -1):
        side_items = sorted([it for it in items if it["side"] == side], key=lambda d: d["cx"])
        occupied = []  # list of (tier, x_right) tracking rightmost extent per tier
        for it in side_items:
            x0 = it["cx"] - it["w"] / 2.0
            x1 = it["cx"] + it["w"] / 2.0
            tier = 0
            # find lowest tier whose current right edge clears this label's left edge
            while any(t == tier and x0 < xr + PAD for t, xr in occupied):
                tier += 1
            occupied.append((tier, x1))
            ytxt = ybase + side * (BASE_GAP + tier * TIER_DY + GENE_H / 2.0)
            va = "bottom" if side > 0 else "top"
            ax.text(it["cx"], ytxt, it["lbl"], ha="center", va=va,
                    fontsize=FS_LABEL, fontstyle="italic")
            # leader line from arrow edge to the label tier
            y_arrow = ybase + side * (GENE_H / 2.0)
            y_lead = ytxt - side * 0.04
            if tier > 0:
                ax.add_line(Line2D([it["cx"], it["cx"]], [y_arrow, y_lead],
                                   color="#bbbbbb", linewidth=0.4, zorder=0))


def _place_text(ax, fig, x, y, s, **kw):
    """Draw text left-aligned at (x, y) in data coords and return the x coord
    just past its right edge, so the next element can be placed without overlap.
    """
    t = ax.text(x, y, s, ha="left", va="bottom", **kw)
    fig.canvas.draw()                       # ensure a renderer exists
    bb = t.get_window_extent(renderer=fig.canvas.get_renderer())
    # convert the text's pixel width to data-x units
    x0d = ax.transData.inverted().transform((bb.x0, bb.y0))[0]
    x1d = ax.transData.inverted().transform((bb.x1, bb.y0))[0]
    return x + (x1d - x0d)


def render(spec_path, out_path, title_mode="gna"):
    panels = _read_spec(spec_path)
    n = len(panels)
    max_genes = max(len(p["genes"]) for p in panels.values())
    xmax = max_genes * (GENE_W + GAP)

    # First pass over layout: assign each panel a baseline y and record where a
    # section starts (so the section band gets its own vertical slot and the
    # divider never overlaps panel A).
    layout = []          # (pid, panel, baseline_y, section_header_or_None)
    y = 0.0
    current_section = None
    for pid, panel in panels.items():
        sec = panel.get("section", "").strip()
        new_sec = sec if (sec and sec != current_section) else None
        if new_sec:
            y -= SECTION_PAD          # reserve room for the band above this panel
            current_section = sec
        layout.append((pid, panel, y, new_sec))
        y -= ROW_DY
    y_bottom = y

    # Figure geometry. keep data units square so pentagons are undistorted, and
    # fix physical width to the manuscript column so it embeds 1:1.
    x_lo, x_hi = -0.4, xmax + 0.4
    y_hi = GENE_H + 3.2
    y_lo = y_bottom + ROW_DY - 3.1
    xspan = x_hi - x_lo
    yspan = y_hi - y_lo
    fig_w = TARGET_WIDTH_IN
    fig_h = fig_w * (yspan / xspan)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    # Fix axis limits and aspect BEFORE any text is drawn, so _place_text's
    # transData measurement uses the final (stable) coordinate transform.
    ax.set_xlim(x_lo, x_hi)
    ax.set_ylim(y_lo, y_hi)
    ax.set_aspect("equal")
    ax.axis("off")

    for pid, panel, ybase, new_sec in layout:
        genes = panel["genes"]
        # section band. title sits in its reserved slot ABOVE the panel header,
        # divider line just under the title, clear of panel A and its label.
        if new_sec:
            title_y = ybase + GENE_H * 1.7 + 2.15
            ax.text(0, title_y, new_sec, ha="left", va="bottom",
                    fontsize=FS_SECTION, fontweight="bold", color="#333333")
            # divider sits just under the title text, with a clear gap above the
            # panel header line below it
            ax.add_line(Line2D([0, xmax], [title_y - 0.12, title_y - 0.12],
                               color="#999999", linewidth=0.9))
        # panel header. place elements sequentially by measuring each one's
        # rendered width, so the panel letter, italic species and roman
        # "substitution - locus" tail never overlap regardless of name length.
        # Sits well above the top label tier so gene labels never reach it.
        header_y = ybase + GENE_H * 1.7 + 1.15
        sub = panel["substitution"].strip()
        locus = panel["locus"].strip()
        tail = f"{sub} - {locus}" if sub else locus
        x_cur = _place_text(ax, fig, 0, header_y, pid + "  ",
                            fontsize=FS_PANEL, fontweight="bold")
        x_cur = _place_text(ax, fig, x_cur, header_y, panel["species"] + "  ",
                            fontsize=FS_SPECIES, fontstyle="italic")
        if tail:
            _place_text(ax, fig, x_cur, header_y, tail, fontsize=FS_TAIL)
        x = 0.0
        centres = []
        for g in genes:
            _draw_gene(ax, x, ybase, g)
            centres.append(x + GENE_W / 2.0)
            x += GENE_W + GAP
        # multi-tier gene labels. labels wider than the gene spacing collide, so
        # stack them into tiers above (and below) the arrows with leader lines,
        # never overwriting a neighbour.
        _place_gene_labels(ax, fig, genes, centres, ybase)

    # legend. anchor at the TRUE horizontal centre of the content (xmax/2 in
    # data coords), not axes-fraction 0.5, so it sits under the panels and the
    # tight bounding box does not pull in empty space on the left.
    handles = [Patch(facecolor=CATEGORY_COLORS[c][0], edgecolor=CATEGORY_COLORS[c][1],
                     label=CATEGORY_LABELS[c]) for c in CATEGORY_COLORS]
    ax.legend(handles=handles, loc="upper center",
              bbox_to_anchor=(xmax / 2.0, y_lo + 1.1),
              bbox_transform=ax.transData, ncol=3, frameon=False,
              fontsize=FS_LEGEND, handlelength=1.4, columnspacing=1.4)

    fig.savefig(out_path, bbox_inches="tight")
    fig.savefig(out_path.rsplit(".", 1)[0] + ".png", dpi=300, bbox_inches="tight")
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
