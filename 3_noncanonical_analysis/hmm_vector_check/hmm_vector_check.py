#!/usr/bin/env python3
"""Per-structure HMM-vs-structure comparison for the six canonical His positions of
the non-canonical pool -- i.e. the data behind Table S6 and the "504/1,060 complete
vectors" claim, rebuilt so individual structures and groups can be looked up.

The published table (3_noncanonical_analysis/hmm_his_comparison.tsv) is a summary
with no surviving generator, so this reproduces it from the inputs named in the
Methods: hmmalign of the PF00264.26 profile (214 match states) over the query
sequences, with PmTYR/B2ZB02 as the reference that ties match-state columns to the
six canonical His positions. Reproduction of the published per-position counts is
asserted, not assumed -- if they do not match, the script says so.

hmmalign is run through pyhmmer (same HMMER3 code) so no external binary is needed.

Outputs: hmm_vs_structure_nc.tsv (one row per structure, both vectors + per-position
agreement), summary.txt
"""
import gzip
import os

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.dirname(os.path.dirname(HERE))
HMM = f"{BASE}/2_canonical_analysis/hmm/PF00264.hmm"
COLS_TSV = f"{BASE}/2_canonical_analysis/hmm/hmm_match_columns.tsv.gz"
FASTA = next((p for p in (f"{BASE}/2_canonical_analysis/hmm/query.fasta.gz",
                          f"{BASE}/2_canonical_analysis/hmm/query.fasta")
              if os.path.exists(p)), None)

POS = ["CuA_His1", "CuA_His2", "CuA_His3", "CuB_His1", "CuB_His2", "CuB_His3"]
PMTYR_HIS = {"CuA_His1": 42, "CuA_His2": 60, "CuA_His3": 69,
             "CuB_His1": 204, "CuB_His2": 208, "CuB_His3": 231}
AA1 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
       "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
       "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
       "TYR": "Y", "VAL": "V"}

# ---------------------------------------------------------------- inputs
pools = pd.read_csv(f"{BASE}/1_filtering/final_pools/three_pool_assignment_final.csv")
nc = list(pools.loc[pools.pool == "noncanonical", "accession"])
struct = pd.read_csv(f"{BASE}/3_noncanonical_analysis/noncanonical_analysis.tsv", sep="\t")
struct = struct.set_index("accession")
# ---------------------------------------------------------------- alignment
# Primary route: the deposited match-state table + reference map written by
# 2_canonical_analysis/hmm/build_alignment.py. Needs nothing but pandas, so the
# 504-of-1,060 claim is checkable from the deposit without running hmmalign.
if not os.path.exists(COLS_TSV):
    raise SystemExit(
        f"{COLS_TSV} not found -- run 2_canonical_analysis/hmm/build_alignment.py first "
        "(it regenerates the PF00264 alignment from the deposited query.fasta.gz)")

refmap = pd.read_csv(f"{BASE}/2_canonical_analysis/hmm/reference_map.tsv", sep="\t")
refmap = refmap.set_index("resnum")
ms_of = {}
for p, r in PMTYR_HIS.items():
    row = refmap.loc[r]
    assert row.residue == "H", f"B2ZB02 residue {r} is {row.residue}, not H ({p})"
    assert not pd.isna(row.match_state), f"B2ZB02 His{r} sits in an insert column"
    ms_of[p] = int(row.match_state)
print("PmTYR His -> PF00264 match state: "
      + ", ".join(f"{p}(His{PMTYR_HIS[p]})={ms_of[p]}" for p in POS))

want = set(nc)
match_cols = {}
with gzip.open(COLS_TSV, "rt") as fh:
    header = fh.readline()
    for line in fh:
        acc, _, cols_str = line.rstrip("\n").partition("\t")
        if acc in want:
            match_cols[acc] = cols_str
missing = want - set(match_cols)
assert not missing, f"{len(missing)} accessions missing from {os.path.basename(COLS_TSV)}"

# Optional extra: the match-state table carries residue identities but not residue
# NUMBERS, so the stricter "HMM maps the same residue, not just the same residue
# type" variant needs the sequences re-aligned. Off by default; needs pyhmmer.
resnums = {}
if os.environ.get("WITH_RESNUMS", "0") not in ("0", "", "no") and FASTA:
    import pyhmmer
    seqs, name, buf = {}, None, []
    opener = gzip.open if FASTA.endswith(".gz") else open
    with opener(FASTA, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name in want:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].split("|")[0].strip(), []
            else:
                buf.append(line.strip())
    if name in want:
        seqs[name] = "".join(buf)
    alpha = pyhmmer.easel.Alphabet.amino()
    with pyhmmer.plan7.HMMFile(HMM) as f:
        hmm = f.read()
    block = [pyhmmer.easel.TextSequence(name=a.encode(), sequence=s).digitize(alpha)
             for a, s in seqs.items()]
    msa = pyhmmer.hmmalign(hmm, block, trim=False)
    rf = msa.reference
    rf = rf.decode() if isinstance(rf, bytes) else rf
    colset = {i for i, c in enumerate(rf) if c != "."}
    order = sorted(colset)
    for nm, al in zip(msa.names, msa.alignment):
        nm = nm.decode() if isinstance(nm, bytes) else nm
        num, per_col = 0, {}
        for i, c in enumerate(al):
            if c in "-.":
                continue
            num += 1
            if i in colset:
                per_col[i] = num
        resnums[nm] = {k + 1: per_col.get(i) for k, i in enumerate(order)}
    print(f"re-aligned {len(resnums):,} sequences for residue numbers")

# ---------------------------------------------------------------- per-structure
rows = []
for acc in nc:
    v = match_cols[acc]
    rnum = resnums.get(acc, {})
    st = struct.loc[acc]
    rec = {"accession": acc}
    for p in POS:
        h = v[ms_of[p] - 1].upper()
        rec[f"{p}_hmm_resnum"] = rnum.get(ms_of[p]) if rnum else None
        sn = st[f"{p}_resnum"]
        rec[f"{p}_struct_resnum"] = None if pd.isna(sn) else int(sn)
        rec[f"{p}_same_residue"] = (None if not rnum else
                                    (rec[f"{p}_hmm_resnum"] is not None
                                     and rec[f"{p}_hmm_resnum"] == rec[f"{p}_struct_resnum"]))
        s3 = str(st[p])
        rec[f"{p}_struct"] = s3
        rec[f"{p}_hmm"] = "-" if h == "-" else h
        rec[f"{p}_struct1"] = "-" if s3 == "---" else AA1.get(s3, "X")
        rec[f"{p}_agree"] = (h != "-") and (h == rec[f"{p}_struct1"])
        rec[f"{p}_gapped"] = h == "-"
        rec[f"{p}_substituted"] = s3 not in ("HIS", "---")
        rec[f"{p}_mapped"] = s3 != "---"
    rows.append(rec)
df = pd.DataFrame(rows)
df["n_agree"] = sum(df[f"{p}_agree"] for p in POS)
df["n_gapped"] = sum(df[f"{p}_gapped"] for p in POS)
df["vector_struct"] = df.apply(lambda r: "".join(r[f"{p}_struct1"] for p in POS), axis=1)
df["vector_hmm"] = df.apply(lambda r: "".join(r[f"{p}_hmm"] for p in POS), axis=1)
# Conventions, in increasing strictness. The first is retained ONLY to explain the
# discrepancy with the manuscript's 504: comparing the two rendered vector strings
# scores an unmapped structural position ("---" -> "-") against an HMM deletion ("-")
# as a match, although nothing was recovered at that position. That contradicts the
# per-position rule behind Table S6, where a deletion never agrees. Do not quote it.
df["full_vector_string_equal"] = df.vector_hmm == df.vector_struct
df["full_correct_type"] = df[[f"{p}_agree" for p in POS]].all(axis=1)
df["full_correct_mapped_only"] = df.apply(
    lambda r: all(r[f"{p}_agree"] for p in POS if r[f"{p}_mapped"]), axis=1)
df["full_correct_subst_only"] = df.apply(
    lambda r: all(r[f"{p}_agree"] for p in POS if r[f"{p}_substituted"]), axis=1)
HAVE_RESNUMS = bool(resnums)
if HAVE_RESNUMS:
    df["full_correct_same_residue"] = df.apply(
        lambda r: all(r[f"{p}_agree"] and r[f"{p}_same_residue"] for p in POS), axis=1)
df.to_csv(f"{HERE}/hmm_vs_structure_nc.tsv", sep="\t", index=False)

# ---------------------------------------------------------------- report
L = []


def say(*s):
    L.append(" ".join(str(x) for x in s))
    print(*s)


pub = pd.read_csv(f"{BASE}/3_noncanonical_analysis/hmm_his_comparison.tsv", sep="\t")
pub = pub.set_index("Position")
say("REPRODUCING Table S6 (published hmm_his_comparison.tsv) on the 1,060 pool")
say(f"{'position':10} {'n_subst':>8} {'aligned':>8} {'gapped':>7} {'agree':>6} "
    f"{'disagree':>9} {'%aln':>6} {'%incl_gaps':>11}   published")
ok = True
for p in POS:
    m = df[df[f"{p}_substituted"]]
    n = len(m)
    g = int(m[f"{p}_gapped"].sum())
    a = int(m[f"{p}_agree"].sum())
    d = n - g - a
    pa = 100 * a / (n - g) if n > g else float("nan")
    pi = 100 * a / n if n else float("nan")
    r = pub.loc[p]
    same = (n == r.n_substituted and n - g == r.HMM_aligned and g == r.HMM_gapped
            and a == r.Agree and d == r.Disagree)
    ok &= bool(same)
    verdict = "MATCH" if same else (
        f"DIFFERS: published {int(r.n_substituted)}/{int(r.HMM_aligned)}/"
        f"{int(r.HMM_gapped)}/{int(r.Agree)}/{int(r.Disagree)}")
    say(f"{p:10} {n:8} {n-g:8} {g:7} {a:6} {d:9} {pa:6.1f} {pi:11.1f}   {verdict}")
say(f"  -> published per-position counts reproduced exactly: {ok}")
say("")
say("COMPLETE SIX-POSITION VECTORS  (manuscript states 504/1,060)")
n_str = int(df.full_vector_string_equal.sum())
n_typ = int(df.full_correct_type.sum())
say(f"  RECOMMENDED -- all six recovered, same residue: "
    f"{int(df.full_correct_same_residue.sum()) if HAVE_RESNUMS else 'n/a'}/1060")
say(f"  all six recovered, residue TYPE only:           {n_typ}/1060")
say(f"  string equality of the rendered vectors:        {n_str}/1060  "
    f"[NOT quotable: {n_str - n_typ} structures scored an unmapped structural")
say(f"                                                   position against an HMM "
    f"deletion as a match]")
say(f"  all STRUCTURALLY MAPPED positions correct:      {int(df.full_correct_mapped_only.sum())}/1060")
say(f"  all SUBSTITUTED positions correct:              {int(df.full_correct_subst_only.sum())}/1060")
if HAVE_RESNUMS:
    tc = df[df.full_correct_type]   # count only inside the type-correct set
    mis = sum(int((tc[f"{p}_agree"] & ~tc[f"{p}_same_residue"]
                   & tc[f"{p}_hmm_resnum"].notna()).sum()) for p in POS)
    say("")
    say(f"  type-only minus same-residue = {n_typ - int(df.full_correct_same_residue.sum())} "
        f"structures, from {mis} position-instances where the HMM reports the right")
    say("  residue type at the wrong residue number (genuine misalignment).")
else:
    say("  (set WITH_RESNUMS=1 for the same-residue check; needs pyhmmer)")
say(f"  structures with >=1 gapped position: {int((df.n_gapped > 0).sum())}; "
    f"median positions correct: {df.n_agree.median():.0f}/6")

with open(f"{HERE}/summary.txt", "w") as fh:
    fh.write("\n".join(L) + "\n")
print(f"\nwrote hmm_vs_structure_nc.tsv, summary.txt")
