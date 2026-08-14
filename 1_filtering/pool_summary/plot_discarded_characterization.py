#!/usr/bin/env python3
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["font.sans-serif"] = ["Helvetica Neue", "Helvetica", "Arial"]
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os, platform

if platform.system() == "Darwin":
    BASE = "/Users/eirinlandsem/Library/Mobile Documents/com~apple~CloudDocs/Proteinkjemi_PhD/Skriving/Thesis/manuscripts:articles/Manuscript_TyrY/New_bioinf/bioinf_redo/Super_reference_pipeline"
else:
    BASE = "/cluster/work/projects/nn1003k/eirin/bioinf/Super_reference_pipeline"

OUTDIR = f"{BASE}/1_filtering/pool_summary"

C_CAN = "#0570B0"
C_NCAN = "#41AB5D"
C_DISC = "#969696"

canon = pd.read_csv(f"{BASE}/canonical_criteria_all_ca.csv")
final = pd.read_csv(f"{BASE}/1_filtering/final_pools/three_pool_assignment_final.csv")
foldseek = pd.read_csv(f"{BASE}/1_filtering/foldseek/foldseek_multiref_best.tsv", sep="\t")
foldseek["accession"] = foldseek["query"].str.replace(r"_taxID_.*", "", regex=True)
taxonomy = pd.read_csv(f"{BASE}/taxonomy_lookup.csv")

df = canon.merge(final[["accession", "pool", "reason"]], on="accession", how="left")
df = df.merge(foldseek[["accession", "qtmscore", "qlen"]], on="accession", how="left")
df = df.merge(taxonomy[["accession", "kingdom"]], on="accession", how="left")
df["pool"] = df["pool"].fillna("discarded")

can = df[df["pool"] == "canonical"]
ncan = df[df["pool"] == "noncanonical"]
disc = df[df["pool"] == "discarded"]

print(f"Canonical: {len(can)}, Non-canonical: {len(ncan)}, Discarded: {len(disc)}")

fig, axes = plt.subplots(2, 2, figsize=(10, 8))
fig.patch.set_facecolor("white")

hist_kw = dict(alpha=0.7, edgecolor="white", linewidth=0.3, density=True)

# A: Sequence length
ax = axes[0, 0]
bins_len = np.linspace(0, 800, 80)
ax.hist(can["qlen"].dropna(), bins=bins_len, color=C_CAN, label=f"Canonical (n={len(can):,})", **hist_kw)
ax.hist(ncan["qlen"].dropna(), bins=bins_len, color=C_NCAN, label=f"Non-canonical (n={len(ncan):,})", **hist_kw)
ax.hist(disc["qlen"].dropna(), bins=bins_len, color=C_DISC, label=f"Discarded (n={len(disc):,})", **hist_kw)
ax.set_xlabel("Sequence length (residues)")
ax.set_ylabel("Density")
ax.set_title("A  Sequence length distribution", loc="left", fontweight="bold")
ax.legend(fontsize=7.5, frameon=False)

# B: qTMscore
ax = axes[0, 1]
bins_qtm = np.linspace(0, 1.0, 60)
ax.hist(can["qtmscore"].dropna(), bins=bins_qtm, color=C_CAN, label="Canonical", **hist_kw)
ax.hist(ncan["qtmscore"].dropna(), bins=bins_qtm, color=C_NCAN, label="Non-canonical", **hist_kw)
ax.hist(disc["qtmscore"].dropna(), bins=bins_qtm, color=C_DISC, label="Discarded", **hist_kw)
ax.set_xlabel("Foldseek qTMscore (best reference)")
ax.set_ylabel("Density")
ax.set_title("B  Structural similarity", loc="left", fontweight="bold")
ax.legend(fontsize=7.5, frameon=False)

# C: Discard reasons
ax = axes[1, 0]
reason_map = {
    "sixHis_fail_geometry": "6 His, fail geometry",
    "corehelix_core_fail": "<6 His, core fail",
    "lt6His_no_foldseek": "<6 His, no Foldseek hit",
}
disc_reasons = disc["reason"].map(reason_map).fillna(disc["reason"])
reason_counts = disc_reasons.value_counts()
colors_r = ["#D94701", "#FD8D3C", "#FDBE85"]
bars = ax.barh(range(len(reason_counts)), reason_counts.values, color=colors_r[:len(reason_counts)],
               edgecolor="white", linewidth=0.3)
ax.set_yticks(range(len(reason_counts)))
ax.set_yticklabels(reason_counts.index, fontsize=9)
ax.set_xlabel("Count")
ax.set_title("C  Discard reasons", loc="left", fontweight="bold")
for bar, val in zip(bars, reason_counts.values):
    ax.text(bar.get_width() + 50, bar.get_y() + bar.get_height()/2,
            f"{val:,} ({100*val/len(disc):.1f}%)", va="center", fontsize=8.5)
ax.set_xlim(0, reason_counts.max() * 1.35)

# D: Taxonomy of discarded
ax = axes[1, 1]
kingdoms = ["Fungi", "Animals", "Plants", "Bacteria", "Oomycota", "Other Eukaryota", "Archaea"]
disc_tax = disc["kingdom"].value_counts()
x_pos = np.arange(len(kingdoms))
width = 0.35
can_tax = can["kingdom"].value_counts()
disc_pct = [100 * disc_tax.get(k, 0) / len(disc) for k in kingdoms]
can_pct = [100 * can_tax.get(k, 0) / len(can) for k in kingdoms]
ax.bar(x_pos - width/2, can_pct, width, color=C_CAN, edgecolor="white", linewidth=0.3, label="Canonical")
ax.bar(x_pos + width/2, disc_pct, width, color=C_DISC, edgecolor="white", linewidth=0.3, label="Discarded")
ax.set_xticks(x_pos)
ax.set_xticklabels(kingdoms, rotation=30, ha="right", fontsize=9)
ax.set_ylabel("Percentage (%)")
ax.set_title("D  Taxonomic distribution", loc="left", fontweight="bold")
ax.legend(fontsize=7.5, frameon=False)

for a in axes.flat:
    a.spines["top"].set_visible(False)
    a.spines["right"].set_visible(False)
    a.tick_params(labelsize=9)

fig.suptitle(
    f"Discarded pool characterization (n={len(disc):,})",
    fontsize=14, fontweight="bold", y=0.99,
)

plt.tight_layout(rect=[0, 0, 1, 0.95])
for fmt in ["pdf", "png"]:
    out = os.path.join(OUTDIR, f"discarded_characterization.{fmt}")
    fig.savefig(out, dpi=300, facecolor="white", bbox_inches="tight")
    print(f"Saved: {out}")
plt.close()

print(f"\n=== Discarded pool statistics ===")
print(f"Total: {len(disc):,}")
print(f"Median length: {disc['qlen'].median():.0f} residues")
print(f"Median qTMscore: {disc['qtmscore'].median():.3f}")
print(f"Reasons: {disc['reason'].value_counts().to_dict()}")
print(f"No Foldseek hit: {disc['qtmscore'].isna().sum()}")
