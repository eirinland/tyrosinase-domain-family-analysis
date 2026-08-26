#!/usr/bin/env python3
"""What if M7 also required the four core-helix anchors to land on four SEPARATE
helices of the query (no helix serving two anchor slots)?

M7 accepts core helix i when at least one of its reference His anchors finds a query
Ca that is close, locally helical and confident. The four helices are tested
independently, so nothing forbids two anchors being satisfied by the SAME piece of
query backbone: a structure that has collapsed two core helices into one, or that is
aligned badly enough for one long helix to serve two slots, can still score 4/4.
This script adds the missing constraint and measures what it would cost.

The rule needs a Ca-only definition of "the same helix", and that definition is
calibrated, not assumed (calibrate_segments.py): a residue counts as helical when it
is covered by >= 3 genuine helical turns (Ca i->i+3 in 4.0-6.4 A and i->i+4 in
4.8-8.2 A, the FINAL widened M7 windows); maximal runs of consecutive helical
residues are the query's helices; an anchor is assigned to the run holding its query
residue or a neighbour (+/-1, M7's own cap fallback). That is the ONLY tested
definition under which all five copper-bearing references give four separate,
correctly bounded core helices. Looser variants (>=1 turn) merge PmTYR's a3 and a4
through the intervening loop and are reported as sensitivity only, flagged by their
reference control.

Scope: the 1,060 retained non-canonical structures -- the rule is an extra AND on
core_ok, so it can only move structures OUT of non-canonical -- plus the 241
hand-labelled benchmark structures, so the recall cost is measured, not assumed.

Run cache_ca.py first. Outputs: separate_helix_nc.tsv, separate_helix_bench.tsv,
summary.txt
"""
import os

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
# pipeline root, two levels up: <root>/1_filtering/separate_helix_check
BASE = os.path.dirname(os.path.dirname(HERE))
# The AF3 model archive is a Zenodo item, not part of the repository. Only cache_ca.py
# needs it, and only when ca_cache.npz is missing.
TARBALL = os.environ.get("PPO_AF3_TARBALL", os.path.join(
    os.path.dirname(os.path.dirname(BASE)), "ppo_af3_models.tar.gz"))

WIDE3, WIDE4 = (4.0, 6.4), (4.8, 8.2)     # FINAL M7 helicity windows
STRICT3, STRICT4 = (4.8, 6.4), (5.4, 7.4)  # pre-widening windows
DMAX, PMIN, NEED = 4.0, 70.0, 1            # M7 anchor gates (unchanged)
PAIRS = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
LOBE = {0: "CuA", 1: "CuA", 2: "CuB", 3: "CuB"}

# key, label, d3 window, d4 window, turns needed per helical residue, gap bridged
DEFS = [
    ("primary", "calibrated: widened windows, >=3 covering turns", WIDE3, WIDE4, 3, 0),
    ("bridge2", "as primary, runs bridged over gaps <=2", WIDE3, WIDE4, 3, 2),
    ("loose", "widened windows, >=1 covering turn", WIDE3, WIDE4, 1, 0),
    ("strict_win", "pre-widening windows, >=1 covering turn", STRICT3, STRICT4, 1, 0),
]

CORE_REF = {
    "ref_PmTYR": [42, 69, 204, 231], "ref_2Y9W_Abisporus": [61, 94, 259, 296],
    "ref_5CE9_Jregia": [87, 117, 243, 273], "ref_1BT3_Ibatatas": [88, 118, 240, 274],
    "ref_1JS8_squid": [2543, 2571, 2671, 2702],
}


# ---------------------------------------------------------------- geometry
def turns(ca, r, w3, w4):
    n = 0
    for j in range(r - 4, r + 1):
        if all(k in ca for k in (j, j + 3, j + 4)):
            d3 = float(np.linalg.norm(ca[j] - ca[j + 3]))
            d4 = float(np.linalg.norm(ca[j] - ca[j + 4]))
            if w3[0] <= d3 <= w3[1] and w4[0] <= d4 <= w4[1]:
                n += 1
    return n


def segments(ca, w3, w4, need, bridge):
    hel = sorted(r for r in ca if turns(ca, r, w3, w4) >= need)
    if not hel:
        return []
    out = [[hel[0], hel[0]]]
    for r in hel[1:]:
        if r - out[-1][1] - 1 <= bridge:
            out[-1][1] = r
        else:
            out.append([r, r])
    return [tuple(x) for x in out]


def seg_of(segs, r):
    if r is None:
        return None
    for rr in (r, r - 1, r + 1):
        for i, (s, e) in enumerate(segs):
            if s <= rr <= e:
                return i
    return None


def parse_pdb_ca(path):
    ca = {}
    for line in open(path):
        if line[:6].strip() in ("ATOM", "HETATM") and line[12:16].strip().upper() == "CA":
            try:
                ca[int(line[22:26])] = np.array([float(line[30:38]), float(line[38:46]),
                                                 float(line[46:54])])
            except ValueError:
                pass
    return ca


def parse_cif_ca(path):
    cols, inA, coll = [], False, False
    ca = {}
    for raw in open(path):
        l = raw.strip()
        if l == "loop_":
            coll, inA, cols = False, False, []
            continue
        if l.startswith("_atom_site."):
            coll = True
            cols.append(l)
            continue
        if coll and cols:
            coll, inA = False, True
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
                ca[int(r["_atom_site.label_seq_id"])] = np.array(
                    [float(r["_atom_site.Cartn_x"]), float(r["_atom_site.Cartn_y"]),
                     float(r["_atom_site.Cartn_z"])])
            except (KeyError, ValueError):
                continue
    return ca


# ---------------------------------------------------------------- inputs
pools = pd.read_csv(f"{BASE}/1_filtering/final_pools/three_pool_assignment_final.csv")
wid = pd.read_csv(f"{BASE}/1_filtering/si_discarded_examples/core_helix_results_widened.tsv",
                  sep="\t")
bench_m7 = pd.read_csv(f"{BASE}/1_filtering/core_helix_filter/core_helix_bench.tsv", sep="\t")
labels = pd.read_csv(f"{BASE}/1_filtering/benchmark/benchmark_results.tsv", sep="\t")

nc = set(pools.loc[pools.pool == "noncanonical", "accession"])
nc_rows = wid[wid.accession.isin(nc)].copy()
assert len(nc_rows) == len(nc) == 1060, (len(nc_rows), len(nc))

bench = bench_m7.merge(labels[["accession", "your_label", "stratum"]], on="accession",
                       how="left")
bench = bench[bench.your_label.isin(["canonical", "noncanonical", "discard"])].copy()

z = np.load(f"{HERE}/ca_cache.npz")
cached = {k.split("|")[0] for k in z.files}


def ca_of(acc):
    if f"{acc}|seq" not in z.files:
        return None, None
    s, x, b = z[f"{acc}|seq"], z[f"{acc}|xyz"], z[f"{acc}|bf"]
    return {int(k): v.astype(float) for k, v in zip(s, x)}, \
        {int(k): float(v) for k, v in zip(s, b)}


# ---------------------------------------------------------------- reference control
ref_ca = {"ref_PmTYR": parse_cif_ca(f"{BASE}/1_filtering/B2ZB02_taxID_1404_model.cif")}
for name in CORE_REF:
    if name != "ref_PmTYR":
        ref_ca[name] = parse_pdb_ca(f"{BASE}/1_filtering/foldseek/{name}.pdb")

ref_control = {}
for key, label, w3, w4, need, bridge in DEFS:
    res = {}
    for name, ca in ref_ca.items():
        segs = segments(ca, w3, w4, need, bridge)
        idx = [seg_of(segs, a) for a in CORE_REF[name]]
        res[name] = len({i for i in idx if i is not None}) == 4
    ref_control[key] = res


# ---------------------------------------------------------------- scoring
def score(df):
    rows = []
    for _, row in df.iterrows():
        q = [None if pd.isna(row[f"a{h}_qres"]) else int(row[f"a{h}_qres"]) for h in (1, 2, 3, 4)]
        rec = {"accession": row.accession, "best_ref": row.best_ref, "best_qtm": row.best_qtm,
               "a1_qres": q[0], "a2_qres": q[1], "a3_qres": q[2], "a4_qres": q[3]}
        ca, bf = ca_of(row.accession)
        if ca is None:
            rows.append(rec)
            continue
        rec["nres"] = len(ca)
        # M7 per-anchor decision under the FINAL widened windows (the bench table on
        # disk is the pre-widening run, so recompute for both scopes identically)
        nok = 0
        for k, h in enumerate((1, 2, 3, 4)):
            if q[k] is None or pd.isna(row[f"a{h}_dist"]):
                continue
            hel = any(turns(ca, rr, WIDE3, WIDE4) >= NEED for rr in (q[k] - 1, q[k], q[k] + 1))
            if float(row[f"a{h}_dist"]) <= DMAX and hel and float(row[f"a{h}_plddt"]) >= PMIN:
                nok += 1
        rec["n_helix_ok_widened"] = nok
        rec["core_ok"] = nok == 4
        sep = [(abs(q[i] - q[j]), i, j) for i, j in PAIRS if q[i] is not None and q[j] is not None]
        if sep:
            s, i, j = min(sep)
            rec["min_idx_sep"] = s
            rec["min_idx_pair"] = f"a{i+1}/a{j+1}"
            rec["min_pair_ca_dist"] = round(float(np.linalg.norm(ca[q[i]] - ca[q[j]])), 2) \
                if q[i] in ca and q[j] in ca else None
        for key, label, w3, w4, need, bridge in DEFS:
            segs = segments(ca, w3, w4, need, bridge)
            idx = [seg_of(segs, x) for x in q]
            coll = [(i, j) for i, j in PAIRS if idx[i] is not None and idx[i] == idx[j]]
            rec[f"{key}_ndistinct"] = len({i for i in idx if i is not None})
            rec[f"{key}_nunassigned"] = sum(1 for i in idx if i is None)
            rec[f"{key}_collision"] = bool(coll)
            rec[f"{key}_pairs"] = ";".join(f"a{i+1}/a{j+1}" for i, j in coll)
            rec[f"{key}_crosslobe"] = any(LOBE[i] != LOBE[j] for i, j in coll)
            rec[f"{key}_segs"] = ";".join("" if i is None else f"{segs[i][0]}-{segs[i][1]}"
                                          for i in idx)
        rows.append(rec)
    return pd.DataFrame(rows)


ncs = score(nc_rows)
bns = score(bench).merge(bench[["accession", "your_label", "stratum"]], on="accession",
                         how="left")

tax = pd.read_csv(f"{BASE}/taxonomy_lookup.csv")
ncs = ncs.merge(tax[["accession", "kingdom", "genus"]], on="accession", how="left")
st3 = pd.read_csv(f"{BASE}/3_noncanonical_analysis/helix_and_gap_filtered_structures.tsv",
                  sep="\t")
ncs = ncs.merge(st3[["accession", "classification"]].rename(
    columns={"classification": "site_tier"}), on="accession", how="left")
ncs["n_his"] = ncs.accession.map(pools.set_index("accession").n_his)

ncs.to_csv(f"{HERE}/separate_helix_nc.tsv", sep="\t", index=False)
bns.to_csv(f"{HERE}/separate_helix_bench.tsv", sep="\t", index=False)

# ---------------------------------------------------------------- report
L = []


def say(*s):
    L.append(" ".join(str(x) for x in s))
    print(*s)


say("EXTRA M7 REQUIREMENT: the four core-helix anchors must sit on four SEPARATE")
say("helices of the query (one helix may not satisfy two anchor slots).")
say(f"M7 gates unchanged: dmax {DMAX} A, pLDDT >= {PMIN}, helicity windows "
    f"d3 {WIDE3} d4 {WIDE4}.")
say("")
say("REFERENCE CONTROL -- does the definition give 4 separate core helices on the")
say("five copper-bearing references (it must, or the rule is measuring an artefact)?")
for key, label, w3, w4, need, bridge in DEFS:
    r = ref_control[key]
    bad = [n for n, ok in r.items() if not ok]
    say(f"  {key:11} {label:52} {sum(r.values())}/5"
        + (f"  MERGED on {', '.join(bad)}" if bad else "  OK"))
say("")
say(f"M7 core_ok recomputed on the retained pool: {int(ncs.core_ok.sum()):,}/{len(ncs):,}"
    " (sanity: must be 1,060)")
say("")
say("EFFECT ON THE POOLS (canonical 21,893 is untouched -- it never uses M7):")
for key, label, w3, w4, need, bridge in DEFS:
    c = ncs[f"{key}_collision"]
    x = ncs[f"{key}_crosslobe"]
    flag = "" if all(ref_control[key].values()) else "   [fails reference control]"
    say(f"  {key:11} collisions {int(c.sum()):4}  ({int((c & x).sum()):3} cross-lobe)"
        f"  ->  non-canonical {len(ncs)-int(c.sum()):5,}   discarded {9116+int(c.sum()):6,}"
        f"{flag}")
    if c.sum():
        say(f"{'':14}colliding pair: {ncs.loc[c, f'{key}_pairs'].value_counts().to_dict()}")
say(f"  {'idx<=10':11} residue-index proxy used by the SI panel-E analysis: "
    f"{int((ncs.min_idx_sep <= 10).sum())}")
say("")
say(f"inter-anchor separation across the retained pool: median "
    f"{ncs.min_idx_sep.median():.0f} residues, median closest-pair Ca-Ca "
    f"{ncs.min_pair_ca_dist.median():.1f} A")
say("")

K = "primary"
drop = ncs[ncs[f"{K}_collision"]].copy()
say(f"=== structures the calibrated rule ({K}) would move to discarded: {len(drop)} ===")
if len(drop):
    for c in ("n_his", "site_tier", "kingdom"):
        say(f"  by {c}: {drop[c].value_counts().sort_index().to_dict()}")
    say(f"  best_qtm median {drop.best_qtm.median():.3f} (pool {ncs.best_qtm.median():.3f});"
        f" length median {drop.nres.median():.0f} (pool {ncs.nres.median():.0f})")
    say(f"  cross-lobe (one helix serving a CuA and a CuB slot): "
        f"{int(drop[f'{K}_crosslobe'].sum())}")
    say("")
    cols = ["accession", "best_ref", "best_qtm", "nres", "n_his", "site_tier",
            f"{K}_pairs", f"{K}_crosslobe", f"{K}_segs", "min_idx_sep",
            "min_pair_ca_dist", "kingdom", "genus"]
    say(drop.sort_values(["min_idx_sep"])[cols].to_string(index=False))
say("")

# ---------------------------------------------------------------- benchmark
say("=== BENCHMARK (241 hand-labelled structures; core present = canonical or "
    "non-canonical) ===")
present = bns.your_label.isin(["canonical", "noncanonical"])


def scoreline(name, call):
    tp = int((call & present).sum()); fp = int((call & ~present).sum())
    fn = int((~call & present).sum()); tn = int((~call & ~present).sum())
    prec = tp / (tp + fp) if tp + fp else float("nan")
    rec = tp / (tp + fn) if tp + fn else float("nan")
    f1 = 2 * prec * rec / (prec + rec) if prec + rec else float("nan")
    say(f"  {name:44} TP={tp:3} FP={fp:2} FN={fn:3} TN={tn:3}  "
        f"Acc={(tp+tn)/len(bns):.3f} Prec={prec:.2f} Rec={rec:.2f} F1={f1:.2f}")


scoreline("M7 as published", bns.core_ok)
for key, label, w3, w4, need, bridge in DEFS:
    scoreline(f"M7 + separate-helix ({key})", bns.core_ok & ~bns[f"{key}_collision"])
lost = bns[bns.core_ok & bns[f"{K}_collision"] & present]
gain = bns[bns.core_ok & bns[f"{K}_collision"] & ~present]
say("")
say(f"  under the calibrated rule: {len(lost)} hand-labelled TRUE cores lost, "
    f"{len(gain)} false positives removed")
if len(lost):
    say("  lost: " + ", ".join(f"{r.accession}({r.your_label[:5]})" for _, r in lost.iterrows()))
if len(gain):
    say("  removed: " + ", ".join(gain.accession))

with open(f"{HERE}/summary.txt", "w") as fh:
    fh.write("\n".join(L) + "\n")
print(f"\nwrote separate_helix_nc.tsv, separate_helix_bench.tsv, summary.txt")
