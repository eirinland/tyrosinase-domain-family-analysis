#!/usr/bin/env python3
"""Does the '504/1,060 complete HMM vectors' set cover the non-canonical cases the
manuscript actually highlights -- and are the H5Pro structures assigned correctly?

Reads hmm_vs_structure_nc.tsv (written by hmm_vector_check.py) and checks, one by
one, every named structure and every group discussed in the non-canonical results
section, plus the H5Pro group in full.

Output: highlighted_cases.tsv, h5pro_cases.tsv, highlighted_summary.txt
"""
import os

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(os.path.dirname(HERE))
POS = ["CuA_His1", "CuA_His2", "CuA_His3", "CuB_His1", "CuB_His2", "CuB_His3"]
LABEL = {p: f"H{i+1}" for i, p in enumerate(POS)}

df = pd.read_csv(f"{HERE}/hmm_vs_structure_nc.tsv", sep="\t").set_index("accession")
# the recommended convention: every position recovered as the same residue (not just the
# same residue type), falling back to type-only if the table was built without residue
# numbers (WITH_RESNUMS=1). Never the rendered-string comparison -- see hmm_vector_check.py.
OK = "full_correct_same_residue" if "full_correct_same_residue" in df.columns else "full_correct_type"
df["full_vector_correct"] = df[OK]
tax = pd.read_csv(f"{BASE}/taxonomy_lookup.csv").set_index("accession")
cls = pd.read_csv(f"{BASE}/3_noncanonical_analysis/helix_and_gap_filtered_structures.tsv",
                  sep="\t").set_index("accession")["classification"]
df["genus"] = tax.genus.reindex(df.index)
df["species"] = cls.index.to_series().map(lambda a: a)  # placeholder, replaced below
sp = pd.read_csv(f"{BASE}/3_noncanonical_analysis/noncanonical_analysis.tsv",
                 sep="\t").set_index("accession")["species"]
df["species"] = sp.reindex(df.index)
df["classification"] = cls.reindex(df.index)

# structures named in the results text (figure representatives and cited proteins)
NAMED = [
    ("A0A9N8JDS5", "Fig 6A binuclear, Aureobasidium vineae"),
    ("A0A9P7QPY3", "Fig 6B binuclear, Colletotrichum scovillei"),
    ("A0A6G0QFF2", "Fig 6C binuclear, Phytophthora fragariae"),
    ("G3JPN7", "Fig 7A mononuclear, Cordyceps militaris, H5Pro representative"),
    ("H1UWR0", "Fig 7B mononuclear, Colletotrichum higginsianum, H6Gln representative"),
    ("A0A0V1MUW3", "Fig 7C mononuclear, Trichinella papuae"),
    ("A0A8J9R8H5", "phomQ1, annotated tyrosine halogenase, H5Pro group"),
    ("A0A142I737", "phomQ1', annotated tyrosine halogenase, H5Pro group"),
    ("H3GEM4", "Fig 8A degenerate, Phytophthora ramorum"),
    ("A0A8H9LF69", "Fig 8B degenerate, Microbispora bryophytorum, PDB 32AE"),
    ("H2KPL1", "Fig 8C degenerate, Clonorchis sinensis"),
]

# groups discussed in the text, defined the same way the text defines them
GROUPS = [
    ("H5Pro (dominant mononuclear)", dict(vector="HHHHPH"), 178),
    ("H6Gln (second mononuclear)", dict(vector="HHHHHQ"), 43),
    ("Trichinella, CuB degenerate", dict(vector="HHHYTQ"), 14),
    ("Phytophthora, CuA replaced", dict(vector="ALNHHH"), 20),
    ("largest degenerate group", dict(vector="YVYQHY"), 34),
    ("Microbispora actinobacteria", dict(genus="Microbispora"), 8),
    ("Clonorchis / Opisthorchis flukes", dict(vector="NHHNHQ"), 6),
    ("binuclear: His at all six positions", dict(vector="HHHHHH"), 69),
]

L = []


def say(*s):
    L.append(" ".join(str(x) for x in s))
    print(*s)


say("Does the complete-vector-correct set cover the highlighted non-canonical cases?")
say(f"pool {len(df)}; complete six-position vectors correct "
    f"{int(df.full_vector_correct.sum())} by convention '{OK}' (manuscript: 504)")
say("")
say("NAMED STRUCTURES")
say(f"{'accession':12} {'struct':8} {'HMM':8} {'ok?':4} {'wrong positions':34} note")
rows = []
for acc, note in NAMED:
    if acc not in df.index:
        say(f"{acc:12} not in the non-canonical pool")
        continue
    r = df.loc[acc]
    bad = []
    for p in POS:
        if not r[f"{p}_agree"]:
            got = r[f"{p}_hmm"]
            bad.append(f"{LABEL[p]} {r[f'{p}_struct1']}->{'gap' if got == '-' else got}")
    say(f"{acc:12} {r.vector_struct:8} {r.vector_hmm:8} "
        f"{'YES' if r.full_vector_correct else 'no':4} {', '.join(bad) or '-':34} {note}")
    rows.append(dict(accession=acc, note=note, vector_struct=r.vector_struct,
                     vector_hmm=r.vector_hmm, full_vector_correct=r.full_vector_correct,
                     wrong_positions="; ".join(bad), classification=r.classification,
                     genus=r.genus, species=r.species))
pd.DataFrame(rows).to_csv(f"{HERE}/highlighted_cases.tsv", sep="\t", index=False)
say("")

say("GROUPS")
say(f"{'group':36} {'n':>4} {'text':>5} {'full ok':>8} {'defining pos ok':>16} "
    f"{'gapped at def':>14}")
for name, sel, expect in GROUPS:
    if "vector" in sel:
        m = df[df.vector_struct == sel["vector"]]
        defpos = [p for i, p in enumerate(POS) if sel["vector"][i] != "H"]
    else:
        m = df[df.genus == sel["genus"]]
        defpos = [p for p in POS]
    n = len(m)
    full = int(m.full_vector_correct.sum())
    if defpos and "vector" in sel:
        okdef = int(m[[f"{p}_agree" for p in defpos]].all(axis=1).sum())
        gapdef = int(m[[f"{p}_gapped" for p in defpos]].any(axis=1).sum())
        d = ",".join(LABEL[p] for p in defpos)
    else:
        okdef = gapdef = -1
        d = "all"
    say(f"{name:36} {n:4} {expect:5} {full:8} "
        f"{(str(okdef) if okdef >= 0 else 'n/a'):>16} "
        f"{(str(gapdef) if gapdef >= 0 else 'n/a'):>14}   defining {d}")
say("")
say("  'text' = the count stated in the manuscript; 'defining pos ok' = HMM recovers")
say("  the substituted position(s) that define the group.")
say("")

# ---------------------------------------------------------------- H5Pro in full
h5 = df[df.vector_struct == "HHHHPH"].copy()
p = "CuB_His2"
h5["hmm_at_H5"] = h5[f"{p}_hmm"]
h5["H5_correct"] = h5[f"{p}_agree"]
h5["H5_gapped"] = h5[f"{p}_gapped"]
h5[["vector_struct", "vector_hmm", "hmm_at_H5", "H5_correct", "H5_gapped",
    "full_vector_correct", "n_agree", "n_gapped", "classification", "genus",
    "species"]].to_csv(f"{HERE}/h5pro_cases.tsv", sep="\t")

say("H5PRO GROUP IN FULL (structural vector H H H H P H)")
say(f"  members: {len(h5)}  (manuscript: 178)")
say(f"  HMM correct at H5 (reads P):        {int(h5.H5_correct.sum())}/{len(h5)}")
say(f"  HMM gapped at H5 (position lost):   {int(h5.H5_gapped.sum())}/{len(h5)}")
say(f"  HMM aligned but wrong residue:      "
    f"{int((~h5.H5_correct & ~h5.H5_gapped).sum())}/{len(h5)}")
if int((~h5.H5_correct & ~h5.H5_gapped).sum()):
    say(f"     residues read instead of P: "
        f"{h5[~h5.H5_correct & ~h5.H5_gapped].hmm_at_H5.value_counts().to_dict()}")
say(f"  complete six-position vector correct: {int(h5.full_vector_correct.sum())}/{len(h5)}")
say(f"  positions correct, distribution: {h5.n_agree.value_counts().sort_index().to_dict()}")
say("  where the other five positions fail (all are structurally His):")
for q in POS:
    if q == p:
        continue
    say(f"    {LABEL[q]:3} correct {int(h5[f'{q}_agree'].sum()):3}/{len(h5)}  "
        f"gapped {int(h5[f'{q}_gapped'].sum()):3}  "
        f"wrong-residue {int((~h5[f'{q}_agree'] & ~h5[f'{q}_gapped']).sum()):3}")
say("  the H5Pro members HMM gets wrong:")
for acc, r in h5[~h5.full_vector_correct].iterrows():
    bad = [f"{LABEL[q]} {r[f'{q}_struct1']}->{'gap' if r[f'{q}_gapped'] else r[f'{q}_hmm']}"
           for q in POS if not r[f"{q}_agree"]]
    say(f"    {acc:12} {r.vector_struct} -> {r.vector_hmm}   {', '.join(bad)}   "
        f"{r.species}")
say("")
say("  the two annotated halogenases:")
for acc in ("A0A8J9R8H5", "A0A142I737"):
    if acc in h5.index:
        r = h5.loc[acc]
        say(f"    {acc}  struct {r.vector_struct}  hmm {r.vector_hmm}  "
            f"H5 {'correct' if r.H5_correct else ('gapped' if r.H5_gapped else 'wrong')}  "
            f"full {'ok' if r.full_vector_correct else 'no'}")

say("")
say("H4 GLUTAMATE group (text: 4 confidently predicted structures, 2 distant species)")
glu = df[df.vector_struct == "HHHEHH"]
say(f"  members {len(glu)}; full vector correct {int(glu.full_vector_correct.sum())}; "
    f"H4 correct {int(glu.CuB_His1_agree.sum())}; H4 gapped {int(glu.CuB_His1_gapped.sum())}")
for acc, r in glu.iterrows():
    say(f"    {acc:12} {r.vector_struct} -> {r.vector_hmm}   {r.species}")
say("")
say("HOW THE COMPLETE-VECTOR SET SPLITS BY SITE CLASS")
for c in ("binuclear", "mononuclear", "no_cu"):
    m = df[df.classification == c]
    say(f"  {c:12} {int(m.full_vector_correct.sum()):4}/{len(m):4} "
        f"({100*m.full_vector_correct.mean():.0f}%)   median positions correct "
        f"{m.n_agree.median():.0f}/6, median gapped {m.n_gapped.median():.0f}")
say("")
say("BY NUMBER OF STRUCTURALLY SUBSTITUTED POSITIONS")
df["n_subst"] = sum(df[f"{p}_substituted"] for p in POS)
for k, m in df.groupby("n_subst"):
    say(f"  {k} substituted: {int(m.full_vector_correct.sum()):4}/{len(m):4} "
        f"({100*m.full_vector_correct.mean():.0f}%)")

with open(f"{HERE}/highlighted_summary.txt", "w") as fh:
    fh.write("\n".join(L) + "\n")
print("\nwrote highlighted_cases.tsv, h5pro_cases.tsv, highlighted_summary.txt")
