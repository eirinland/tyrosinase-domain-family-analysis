#!/usr/bin/env python3
"""Combined stacked bar + confident-Cu bar: AF3 vs M3D vs PMM metal predictions.

Panel A: stacked bars — 3 bars per group (AF3 Cu-only; M3D Cu+other; PMM Cu+FECONI+Zn).
Panel B: confident Cu — sites with Cu predicted above confidence threshold
         (AF3 pLDDT >= 70, M3D prob >= 0.7, PMM Cu any).

Colorblind-safe palette: blue (AF3) / warm gold (M3D) / purple (PMM).
"""
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
M3D_COMP = HERE / "af3_vs_metal3d_comparison.tsv"
PMM_CU_TSV = HERE / "pmm_per_site.tsv"
PMM_FE_TSV = HERE / "pmm_feconi_per_site.tsv"
PMM_RAW = NC_DIR / "pinmymetal" / "pmm_nc_results.tsv"
NC_TSV = NC_DIR / "noncanonical_analysis.tsv"
OUT_PNG = HERE / "stacked_combined_figure.png"
OUT_PDF = HERE / "stacked_combined_figure.pdf"

for fname in ["Helvetica Neue", "Helvetica", "Arial"]:
    matches = fm.findSystemFonts()
    if any(fname.lower().replace(" ", "") in f.lower().replace(" ", "") for f in matches):
        plt.rcParams["font.family"] = "sans-serif"
        plt.rcParams["font.sans-serif"] = [fname]
        break

GROUPS = ["3His", "2His+coord", "2His+inert", "1His", "0His"]
GROUP_LABELS = ["3 His\n(canonical)", "2 His +\ncoord.", "2 His +\nnoncoord.", "1 His", "0 His"]

AF3_THRESH = 70.0
M3D_THRESH = 0.7

def lighten(c, f=0.45):
    return tuple(ci + (1 - ci) * f for ci in c)

def darken(c, f=0.3):
    return tuple(ci * (1 - f) for ci in c)

AF3_CU = (0x5B/255, 0x8F/255, 0xAE/255)
M3D_CU = (0xC4/255, 0xA8/255, 0x78/255)
PMM_CU = (0xB4/255, 0x9E/255, 0xBF/255)
M3D_OTHER = lighten(M3D_CU, 0.5)
PMM_FECONI = lighten(PMM_CU, 0.25)
PMM_ZN = lighten(PMM_CU, 0.55)
AF3_EDGE = darken(AF3_CU, 0.35)
M3D_EDGE = darken(M3D_CU, 0.35)
PMM_EDGE = darken(PMM_CU, 0.35)

# ── Load data ────────────────────────────────────────────────────────
acc_site_group = {}
bar = {g: {"total": 0, "af3": 0, "af3_conf": 0,
           "m3d_cu": 0, "m3d_cu_conf": 0, "m3d_other": 0,
           "pmm_cu": 0, "pmm_fe": 0, "pmm_zn": 0} for g in GROUPS}

with open(TSV) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        g = row["chem_group"]
        if g not in bar:
            continue
        bar[g]["total"] += 1
        key = (row["accession"], row["site"])
        acc_site_group[key] = g
        if row["af3_plddt"]:
            bar[g]["af3"] += 1
            if float(row["af3_plddt"]) >= AF3_THRESH:
                bar[g]["af3_conf"] += 1
        if row["m3d_cu_prob"]:
            bar[g]["m3d_cu"] += 1
            if float(row["m3d_cu_prob"]) >= M3D_THRESH:
                bar[g]["m3d_cu_conf"] += 1

# M3D non-Cu element predictions
m3d_site_element = {}
with open(M3D_COMP) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        acc = row["accession"]
        a1 = row.get("af3_cu1_assignment", "").strip()
        a2 = row.get("af3_cu2_assignment", "").strip()
        e1 = row.get("m3d_cu1_closest_any_element", "").strip()
        e2 = row.get("m3d_cu2_closest_any_element", "").strip()
        if a1 in ("CuA", "CuB") and e1:
            m3d_site_element[(acc, a1)] = e1
        if a2 in ("CuA", "CuB") and e2:
            m3d_site_element[(acc, a2)] = e2
        eA = row.get("m3d_CuA_any_element", "").strip()
        eB = row.get("m3d_CuB_any_element", "").strip()
        if eA and (acc, "CuA") not in m3d_site_element:
            m3d_site_element[(acc, "CuA")] = eA
        if eB and (acc, "CuB") not in m3d_site_element:
            m3d_site_element[(acc, "CuB")] = eB

for key, g in acc_site_group.items():
    elem = m3d_site_element.get(key, "")
    if elem and elem != "CU":
        bar[g]["m3d_other"] += 1

# PMM Cu and FECONI
pmm_assigned = set()
with open(PMM_CU_TSV) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        key = (row["accession"], row["site"])
        g = acc_site_group.get(key)
        if g is None:
            continue
        bar[g]["pmm_cu"] += 1
        pmm_assigned.add(key)

with open(PMM_FE_TSV) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        key = (row["accession"], row["site"])
        g = acc_site_group.get(key)
        if g is None:
            continue
        bar[g]["pmm_fe"] += 1
        pmm_assigned.add(key)

# PMM Zinc from raw
CUA_COLS = ["CuA_His1_resnum", "CuA_His2_resnum", "CuA_His3_resnum"]
CUB_COLS = ["CuB_His1_resnum", "CuB_His2_resnum", "CuB_His3_resnum"]
canon_pos = {}
with open(NC_TSV) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        acc = row["accession"]
        cua = {int(row[c]) for c in CUA_COLS if row.get(c, "").strip()}
        cub = {int(row[c]) for c in CUB_COLS if row.get(c, "").strip()}
        canon_pos[acc] = (cua, cub)

def overlap(ligs, refs, tol=1):
    return sum(1 for lr in ligs if any(abs(lr - cr) <= tol for cr in refs))

def assign_site(ligs, cua, cub):
    oa, ob = overlap(ligs, cua), overlap(ligs, cub)
    if oa == 0 and ob == 0:
        return None
    if oa > ob: return "CuA"
    if ob > oa: return "CuB"
    return "CuA" if min(ligs) < (min(cua | cub) + max(cua | cub)) / 2 else "CuB"

with open(PMM_RAW) as f:
    for row in csv.DictReader(f, delimiter="\t"):
        acc = row["accession"]
        if acc not in canon_pos:
            continue
        cua, cub = canon_pos[acc]
        if not cua and not cub:
            continue
        probs = row.get("all_probs", "")
        ligs_raw = row.get("ligands", "")
        if not probs or not ligs_raw:
            continue
        for i, prob_entry in enumerate(probs.split(";")):
            site_ligs = ligs_raw.split(" | ")
            if i >= len(site_ligs):
                break
            metal, prob = prob_entry.rsplit(":", 1)
            if metal.strip() != "Zinc":
                continue
            lig_resnums = []
            for lig in site_ligs[i].split(";"):
                parts = lig.strip().split(",")
                if parts and parts[0].strip().isdigit():
                    lig_resnums.append(int(parts[0].strip()))
            if not lig_resnums:
                continue
            site = assign_site(lig_resnums, cua, cub)
            if site is None:
                continue
            key = (acc, site)
            if key in pmm_assigned:
                continue
            g = acc_site_group.get(key)
            if g is None:
                continue
            bar[g]["pmm_zn"] += 1
            pmm_assigned.add(key)

# ── Print summary ────────────────────────────────────────────────────
print("Group summary:")
for g in GROUPS:
    d = bar[g]
    t = d["total"]
    print(f"  {g} (n={t}): AF3={d['af3']}({100*d['af3']/t:.0f}%) conf={d['af3_conf']}({100*d['af3_conf']/t:.0f}%) "
          f"M3D Cu={d['m3d_cu']}({100*d['m3d_cu']/t:.0f}%) conf={d['m3d_cu_conf']}({100*d['m3d_cu_conf']/t:.0f}%) other={d['m3d_other']}({100*d['m3d_other']/t:.0f}%) "
          f"PMM Cu={d['pmm_cu']}({100*d['pmm_cu']/t:.0f}%)+FE={d['pmm_fe']}({100*d['pmm_fe']/t:.0f}%)+Zn={d['pmm_zn']}({100*d['pmm_zn']/t:.0f}%)")

# ── Figure ───────────────────────────────────────────────────────────
fig, (ax_a, ax_b) = plt.subplots(2, 1, figsize=(8.5, 7),
    gridspec_kw={"height_ratios": [1, 1], "hspace": 0.10})
fig.subplots_adjust(left=0.10, right=0.95, bottom=0.10, top=0.94)

x = np.arange(len(GROUPS))
bw = 0.24

# ── Panel A: stacked bars (all metal types) ──────────────────────────
for i, g in enumerate(GROUPS):
    t = bar[g]["total"]
    af3 = 100 * bar[g]["af3"] / t
    ax_a.bar(i - bw, af3, bw, color=AF3_CU, edgecolor=AF3_EDGE, lw=0.8, zorder=3)

    m_cu = 100 * bar[g]["m3d_cu"] / t
    m_ot = 100 * bar[g]["m3d_other"] / t
    ax_a.bar(i, m_cu, bw, color=M3D_CU, edgecolor=M3D_EDGE, lw=0.8, zorder=3)
    ax_a.bar(i, m_ot, bw, bottom=m_cu, color=M3D_OTHER, edgecolor=M3D_EDGE, lw=0.8, zorder=3)

    p_cu = 100 * bar[g]["pmm_cu"] / t
    p_fe = 100 * bar[g]["pmm_fe"] / t
    p_zn = 100 * bar[g]["pmm_zn"] / t
    ax_a.bar(i + bw, p_cu, bw, color=PMM_CU, edgecolor=PMM_EDGE, lw=0.8, zorder=3)
    ax_a.bar(i + bw, p_fe, bw, bottom=p_cu, color=PMM_FECONI, edgecolor=PMM_EDGE, lw=0.8, zorder=3)
    ax_a.bar(i + bw, p_zn, bw, bottom=p_cu + p_fe, color=PMM_ZN, edgecolor=PMM_EDGE, lw=0.8, zorder=3)

    for dx, val in [(-bw, af3), (0, m_cu + m_ot), (bw, p_cu + p_fe + p_zn)]:
        if val > 0:
            ax_a.text(i + dx, val + 1.5, f"{val:.0f}%", ha="center", va="bottom",
                      fontsize=6.5, fontweight="bold", color="#333333")
    ax_a.text(i, -7, f"n={t}", ha="center", va="top", fontsize=7.5, color="#555555")

ax_a.set_ylim(-12, 115)
ax_a.set_xlim(-0.55, len(GROUPS) - 0.45)
ax_a.set_xticks(x)
ax_a.set_xticklabels([])
ax_a.set_ylabel("Sites with metal predicted (%)", fontsize=10)
ax_a.spines["top"].set_visible(False)
ax_a.spines["right"].set_visible(False)
ax_a.spines["bottom"].set_visible(False)
ax_a.tick_params(bottom=False)

leg_handles = [
    Patch(fc=AF3_CU, ec=AF3_EDGE, lw=0.8, label="AF3 Cu"),
    Patch(fc=M3D_CU, ec=M3D_EDGE, lw=0.8, label="AllM3D Cu"),
    Patch(fc=M3D_OTHER, ec=M3D_EDGE, lw=0.8, label="AllM3D other"),
    Patch(fc=PMM_CU, ec=PMM_EDGE, lw=0.8, label="PMM Cu"),
    Patch(fc=PMM_FECONI, ec=PMM_EDGE, lw=0.8, label="PMM Fe/Co/Ni"),
    Patch(fc=PMM_ZN, ec=PMM_EDGE, lw=0.8, label="PMM Zn"),
]
ax_a.legend(handles=leg_handles, fontsize=7.5, loc="upper right", frameon=False,
            ncol=2, columnspacing=1.0, handlelength=1.2, handleheight=0.9)

# ── Panel B: confident Cu only ───────────────────────────────────────
for i, g in enumerate(GROUPS):
    t = bar[g]["total"]
    af3_c = 100 * bar[g]["af3_conf"] / t
    m3d_c = 100 * bar[g]["m3d_cu_conf"] / t
    pmm_c = 100 * bar[g]["pmm_cu"] / t

    ax_b.bar(i - bw, af3_c, bw, color=AF3_CU, edgecolor=AF3_EDGE, lw=0.8, zorder=3)
    ax_b.bar(i, m3d_c, bw, color=M3D_CU, edgecolor=M3D_EDGE, lw=0.8, zorder=3)
    ax_b.bar(i + bw, pmm_c, bw, color=PMM_CU, edgecolor=PMM_EDGE, lw=0.8, zorder=3)

    for dx, val in [(-bw, af3_c), (0, m3d_c), (bw, pmm_c)]:
        if val > 0:
            ax_b.text(i + dx, val + 1.5, f"{val:.0f}%", ha="center", va="bottom",
                      fontsize=6.5, fontweight="bold", color="#333333")
    ax_b.text(i, -7, f"n={t}", ha="center", va="top", fontsize=7.5, color="#555555")

ax_b.set_ylim(-12, 115)
ax_b.set_xlim(-0.55, len(GROUPS) - 0.45)
ax_b.set_xticks(x)
ax_b.set_xticklabels(GROUP_LABELS, fontsize=9)
ax_b.set_ylabel("Sites with confident Cu (%)", fontsize=10)
ax_b.spines["top"].set_visible(False)
ax_b.spines["right"].set_visible(False)

leg_b = [
    Patch(fc=AF3_CU, ec=AF3_EDGE, lw=0.8, label=f"AF3 (pLDDT ≥ {AF3_THRESH:.0f})"),
    Patch(fc=M3D_CU, ec=M3D_EDGE, lw=0.8, label=f"AllM3D (prob ≥ {M3D_THRESH})"),
    Patch(fc=PMM_CU, ec=PMM_EDGE, lw=0.8, label=f"PMM Cu (prob ≥ {M3D_THRESH})"),
]
ax_b.legend(handles=leg_b, fontsize=7.5, loc="upper right", frameon=False,
            handlelength=1.2, handleheight=0.9)

ax_a.text(-0.07, 1.02, "A", transform=ax_a.transAxes, fontsize=18, fontweight="bold", va="top")
ax_b.text(-0.07, 1.02, "B", transform=ax_b.transAxes, fontsize=18, fontweight="bold", va="top")

plt.savefig(OUT_PNG, dpi=300, bbox_inches="tight")
plt.savefig(OUT_PDF, bbox_inches="tight")
print(f"\nSaved {OUT_PNG.name} and {OUT_PDF.name}")
