#!/usr/bin/env python3
"""Rebuild the per-anchor a1-a4 detail under the FINAL (widened) M7 helicity windows.

Why this exists
---------------
The local mirror's `1_filtering/core_helix_filter/core_helix_results.tsv` is the
PRE-widening run: it has 1,036 core_ok, and `core_helix_check.py` still carries the
strict windows HELIX_D3 (4.8, 6.4) / HELIX_D4 (5.4, 7.4). The authoritative pools
(`three_pool_assignment_final.csv`, 1,060 non-canonical) are the post-widening run
with HELIX_D3 (4.0, 6.4) / HELIX_D4 (4.8, 8.2). The widened per-anchor table was
never mirrored off Olivia and Olivia's work area has since been purged.

The widening changed ONLY the helicity predicate, which is a pure function of the
query's own Ca coordinates and the anchor residue index `a*_qres`. The Foldseek
alignment, the reference anchors, and therefore a*_anchor / a*_qres / a*_dist /
a*_plddt are unchanged. So the widened a*_helical / a*_sat / n_helix_ok can be
recomputed exactly from the stored a*_qres plus the AF3 coordinates, with no
Foldseek run.

Self-check: the recomputed n_helix_ok / core_ok must reproduce the authoritative
`three_pool_assignment_final.csv`. Residual mismatches are reported, and are
expected only where widening promotes a DIFFERENT anchor of the same helix to
"best" (the table stores just one candidate per helix), which cannot be recovered
without the alignment.

Reads the 2.55 GB AF3 model archive as a single stream (no extraction; the local
disk has no room for 9,898 CIFs).

Output: core_helix_results_widened.tsv, widening_check.txt
"""
import os
import platform
import tarfile

import numpy as np
import pandas as pd

if platform.system() == "Darwin":
    BASE = ("/Users/eirinlandsem/Library/Mobile Documents/com~apple~CloudDocs/"
            "Proteinkjemi_PhD/Skriving/Thesis/manuscripts:articles/Manuscript_TyrY/"
            "New_bioinf/bioinf_redo/Super_reference_pipeline")
    TARBALL = ("/Users/eirinlandsem/Library/Mobile Documents/com~apple~CloudDocs/"
               "Proteinkjemi_PhD/Skriving/Thesis/manuscripts:articles/Manuscript_TyrY/"
               "New_bioinf/ppo_af3_models.tar.gz")
else:
    BASE = "/cluster/work/projects/nn1003k/eirin/bioinf/Super_reference_pipeline"
    TARBALL = os.environ.get("PPO_AF3_TARBALL", "")

HERE = os.path.dirname(os.path.abspath(__file__))

# FINAL widened windows (2026-06-21 sweep, FP=0 across the whole grid)
HELIX_D3 = (4.0, 6.4)
HELIX_D4 = (4.8, 8.2)
DMAX, PMIN, NEED = 4.0, 70.0, 1


def parse_cif_ca(fh):
    """Ca coordinate + B-factor (pLDDT) maps keyed on label_seq_id. Mirrors
    core_helix_check.parse_cif/ca_map/ca_bfac_map."""
    cols, inA, coll = [], False, False
    ca, bf = {}, {}
    for raw in fh:
        l = raw.decode("utf-8", "replace").strip() if isinstance(raw, bytes) else raw.strip()
        if l == "loop_":
            coll, inA, cols = False, False, []
            continue
        if l.startswith("_atom_site."):
            coll = True
            cols.append(l)
            continue
        if coll and cols:
            coll = False
            inA = True
        if inA:
            if l.startswith("_") or l == "#" or not l:
                break
            p = l.split()
            if len(p) != len(cols):
                continue
            r = dict(zip(cols, p))
            try:
                if r["_atom_site.label_atom_id"].upper() != "CA":
                    continue
                s = r["_atom_site.label_seq_id"]
                if not s.lstrip("-").isdigit():
                    continue
                s = int(s)
                ca[s] = np.array([float(r["_atom_site.Cartn_x"]),
                                  float(r["_atom_site.Cartn_y"]),
                                  float(r["_atom_site.Cartn_z"])])
                bf[s] = float(r.get("_atom_site.B_iso_or_equiv", "0") or 0)
            except (KeyError, ValueError):
                continue
    return ca, bf


def helical_turns(qca, r):
    cnt = 0
    for j in range(r - 4, r + 1):
        if all(k in qca for k in (j, j + 3, j + 4)):
            d3 = float(np.linalg.norm(qca[j] - qca[j + 3]))
            d4 = float(np.linalg.norm(qca[j] - qca[j + 4]))
            if HELIX_D3[0] <= d3 <= HELIX_D3[1] and HELIX_D4[0] <= d4 <= HELIX_D4[1]:
                cnt += 1
    return cnt


def residue_helical(qca, r, need=NEED):
    if r is None:
        return False
    return any(helical_turns(qca, rr) >= need for rr in (r, r - 1, r + 1))


# ---------------------------------------------------------------- inputs
strict = pd.read_csv(f"{BASE}/1_filtering/core_helix_filter/core_helix_results.tsv", sep="\t")
pools = pd.read_csv(f"{BASE}/1_filtering/final_pools/three_pool_assignment_final.csv")

want = {a: i for i, a in enumerate(strict.accession)}
qres = {}
for i, row in strict.iterrows():
    qres[row.accession] = [row[f"a{h}_qres"] for h in (1, 2, 3, 4)]

print(f"streaming {TARBALL} for {len(want):,} scope structures ...")
hel = {}
seen = 0
with tarfile.open(TARBALL, "r|gz") as tf:
    for m in tf:
        if not m.name.endswith("_model.cif"):
            continue
        acc = os.path.basename(m.name).split("_taxID_")[0]
        if acc not in want:
            continue
        f = tf.extractfile(m)
        ca, bf = parse_cif_ca(f)
        out = []
        for q in qres[acc]:
            if pd.isna(q):
                out.append((None, None))
            else:
                q = int(q)
                out.append((residue_helical(ca, q), bf.get(q, 0.0)))
        hel[acc] = out
        seen += 1
        if seen % 1000 == 0:
            print(f"  {seen:,}/{len(want):,}")
print(f"  done: {seen:,} structures read")

# ---------------------------------------------------------------- rebuild table
wid = strict.copy()
n_ok = []
for i, row in wid.iterrows():
    got = hel.get(row.accession)
    nok = 0
    for k, h in enumerate((1, 2, 3, 4)):
        if got is None or got[k][0] is None:
            wid.at[i, f"a{h}_helical"] = 0
            wid.at[i, f"a{h}_sat"] = 0
            continue
        helical, plddt = got[k]
        sat = bool(row[f"a{h}_dist"] <= DMAX) and bool(helical) and float(plddt) >= PMIN
        wid.at[i, f"a{h}_helical"] = int(bool(helical))
        wid.at[i, f"a{h}_sat"] = int(sat)
        nok += int(sat)
    n_ok.append(nok)
wid["n_helix_ok"] = n_ok
wid["core_ok"] = wid.n_helix_ok == 4

# ---------------------------------------------------------------- self-check
chk = wid[["accession", "n_helix_ok", "core_ok"]].merge(
    pools[["accession", "pool", "n_helix_ok", "core_ok"]], on="accession",
    suffixes=("_new", "_auth"))
ok_n = (chk.n_helix_ok_new == chk.n_helix_ok_auth).sum()
ok_c = (chk.core_ok_new == chk.core_ok_auth).sum()
old = strict[["accession", "n_helix_ok", "core_ok"]].merge(
    pools[["accession", "n_helix_ok", "core_ok"]], on="accession", suffixes=("_old", "_auth"))

lines = [
    "Recomputed the M7 per-anchor detail under the FINAL widened windows",
    f"  HELIX_D3 {HELIX_D3}   HELIX_D4 {HELIX_D4}   dmax {DMAX}  pmin {PMIN}  need {NEED}",
    "",
    f"scope rows                              {len(wid):,}",
    f"core_ok (recomputed, widened)           {int(wid.core_ok.sum()):,}   [target 1,060]",
    f"core_ok (mirrored strict-threshold run) {int(strict.core_ok.sum()):,}",
    "",
    "agreement with three_pool_assignment_final.csv (authoritative):",
    f"  n_helix_ok  recomputed {ok_n:,}/{len(chk):,} ({100*ok_n/len(chk):.2f}%)"
    f"   vs mirrored strict {(old.n_helix_ok_old == old.n_helix_ok_auth).sum():,}/{len(old):,}",
    f"  core_ok     recomputed {ok_c:,}/{len(chk):,} ({100*ok_c/len(chk):.2f}%)"
    f"   vs mirrored strict {(old.core_ok_old == old.core_ok_auth).sum():,}/{len(old):,}",
]
bad = chk[chk.n_helix_ok_new != chk.n_helix_ok_auth]
if len(bad):
    lines += ["", f"residual n_helix_ok mismatches: {len(bad)}",
              "  (expected only where widening promotes a different anchor of the same",
              "   helix to 'best'; the mirrored table stores one candidate per helix)",
              bad.head(20).to_string(index=False)]
txt = "\n".join(lines)
print("\n" + txt)

wid.to_csv(f"{HERE}/core_helix_results_widened.tsv", sep="\t", index=False)
with open(f"{HERE}/widening_check.txt", "w") as fh:
    fh.write(txt + "\n")
print(f"\nwrote {HERE}/core_helix_results_widened.tsv, widening_check.txt")
