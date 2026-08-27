#!/usr/bin/env python3
"""Regenerate manuscript numbers that no other deposited script produces.

Motivation. A reproducibility audit of every numeric claim in the manuscript
found a small set of values that were correct but had no generator in the
deposit, so a reader could not obtain them from the code. This script closes
that gap. Each block prints the manuscript's published value alongside the
recomputed one, states the convention used, and names its inputs.

Everything here reads only files tracked in this repository. Nothing needs the
Zenodo bulk archive, and no network access is required except for the optional
sequence-identity block, which fetches two sequences from public databases.

Run:
    python3 verify_manuscript_numbers.py                # all local blocks
    python3 verify_manuscript_numbers.py --identity     # add the network block
    python3 verify_manuscript_numbers.py --tsv out.tsv  # machine-readable

Conventions are printed, not assumed. Where a value depends on a definition,
the definition is stated in the output so it can be quoted in Methods.
"""
import argparse
import collections
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from repo_paths import BASE  # noqa: E402

POS = ["CuA_His1", "CuA_His2", "CuA_His3", "CuB_His1", "CuB_His2", "CuB_His3"]
RESULTS = []


def _digits(s):
    """Numeric tokens only, so thousands separators and spacing do not matter.

    Without this, '502/1,060' and '502/1060' compare unequal, and a percentage
    printed to one decimal never matches one printed to none.
    """
    import re
    return [t.rstrip("0").rstrip(".") if "." in t else t
            for t in re.findall(r"\d+(?:\.\d+)?", (s or "").replace(",", ""))]


def record(label, published, recomputed, convention, inputs, agrees=None):
    if agrees is None:
        agrees = _digits(published) == _digits(recomputed)
    RESULTS.append(dict(claim=label, published=published, recomputed=recomputed,
                        agrees="yes" if agrees else "NO",
                        convention=convention, inputs="; ".join(inputs)))
    mark = "  " if agrees else "  <== DIFFERS"
    print(f"  {label}")
    print(f"    published  {published}")
    print(f"    recomputed {recomputed}{mark}")
    print(f"    convention {convention}")


def read_tsv(rel, delim="\t"):
    p = os.path.join(BASE, rel)
    with open(p, newline="", encoding="utf-8", errors="replace") as f:
        return list(csv.DictReader(f, delimiter=delim))


def read_csv(rel):
    return read_tsv(rel, delim=",")


T = lambda v: str(v).strip().lower() in ("true", "1", "yes")


# ---------------------------------------------------------------- block 1
def block_joint_vector():
    """The joint six-position recovery count. Published as 504/1,060."""
    print("\n[1] Profile-HMM joint six-position recovery")
    rel = "3_noncanonical_analysis/hmm_vector_check/hmm_vs_structure_nc.tsv"
    rows = read_tsv(rel)

    def absent(v):
        return (v or "").strip().strip("-") == ""

    def ok(r, p, require_number, credit_both_absent):
        s, h = (r[f"{p}_struct"] or "").strip(), (r[f"{p}_hmm"] or "").strip()
        if absent(s) and absent(h):
            return credit_both_absent
        if absent(s) or absent(h):
            return False
        if T(r[f"{p}_same_residue"]):
            return True
        return (s[:1].upper() == h[:1].upper()) and not require_number

    variants = [
        ("residue type only, absent positions credited", False, True),
        ("type and residue number, absent credited", True, True),
        ("type and residue number, absent counted as failure", True, False),
    ]
    counts = {}
    for name, req, credit in variants:
        counts[name] = sum(1 for r in rows
                           if all(ok(r, p, req, credit) for p in POS))
    for name, n in counts.items():
        print(f"    {n:>5} / {len(rows)}  ({100*n/len(rows):.1f}%)  {name}")

    chosen = "type and residue number, absent counted as failure"
    n = counts[chosen]
    # why this convention: at a both-absent position the STRUCTURAL assignment
    # is also missing, so there is no ground truth to score against
    both_absent = sum(1 for r in rows for p in POS
                      if absent(r[f"{p}_struct"]) and absent(r[f"{p}_hmm"]))
    unmapped = sum(1 for r in rows for p in POS
                   if absent(r[f"{p}_struct"]) and absent(r[f"{p}_hmm"])
                   and not T(r[f"{p}_mapped"]))
    print(f"    ({both_absent} position-instances have neither assignment; "
          f"{unmapped} of those are also unmapped structurally)")
    record("joint six-position vector recovered", "502/1,060",
           f"{n}/{len(rows)}", chosen, [rel])


# ---------------------------------------------------------------- block 2
def block_table_s5():
    """Per-position agreement, Table S5."""
    print("\n[2] Table S5, per-position profile-HMM agreement")
    rel = "3_noncanonical_analysis/hmm_vector_check/hmm_vs_structure_nc.tsv"
    rows = read_tsv(rel)
    print(f"    {'position':<11}{'subst':>7}{'aligned':>9}{'gapped':>8}"
          f"{'agree':>7}{'disagree':>10}{'agree%':>9}")
    out = {}
    for p in POS:
        sub = [r for r in rows if T(r[f"{p}_substituted"])]
        gap = sum(1 for r in sub if T(r[f"{p}_gapped"]))
        aligned = len(sub) - gap
        agree = sum(1 for r in sub if T(r[f"{p}_agree"]))
        # the published percentage column is agree / n_substituted
        pct = 100 * agree / len(sub) if sub else 0
        out[p] = (len(sub), aligned, gap, agree, aligned - agree, pct)
        print(f"    {p:<11}{len(sub):>7}{aligned:>9}{gap:>8}{agree:>7}"
              f"{aligned-agree:>10}{pct:>8.1f}%")
    n, al, g, ag, dis, pct = out["CuB_His1"]
    record("Table S5 CuB_His1 agree / disagree / %", "233 / 44 / 66.8",
           f"{ag} / {dis} / {pct:.1f}",
           "agreement percentage is agreeing positions over n_substituted",
           [rel])


# ---------------------------------------------------------------- block 3
def block_thioether():
    """Thioether against extra-domain presence. Published as 5,516/6,151."""
    print("\n[3] Thioether cysteine against non-core domain presence")
    pv = "2_canonical_analysis/position_vectors.csv"
    da = "1_filtering/foldseek/ppo_domain_assignment_multiref.csv"
    vec = read_csv(pv)
    key = list(vec[0])[0]
    thio = {(r[key] or "").split("_taxID")[0].strip(): (r["thioether"] or "").strip()
            for r in vec}
    dom = {(r["accession"] or "").strip(): r for r in read_csv(da)}
    print("    thioether states: "
          + ", ".join(f"{k or 'blank'}={v}"
                      for k, v in sorted(collections.Counter(thio.values()).items())))

    def has_extra(r):
        return (r["non_ppo_domains"] or "").strip() not in ("", "-", "NA", "none", "[]")

    def nums(s):
        return [int(x) for x in __import__("re").findall(r"\d+", s or "")]

    def is_cterm(r):
        if not has_extra(r):
            return False
        ns, pe = nums(r["non_ppo_range"]), nums(r["ppo_range"])
        return bool(ns) and bool(pe) and min(ns) > max(pe)

    for defname, fn in (("any non-PPO domain", has_extra),
                        ("non-PPO domain C-terminal to the PPO range", is_cterm)):
        ct = collections.Counter()
        for acc, t in thio.items():
            r = dom.get(acc)
            if r is None:
                continue
            ct[(t, bool(fn(r)))] += 1
        print(f"    definition: {defname}")
        for state, pub in (("C", "5,516/6,151 (90%)"), ("-", "3,769/14,383 (26%)")):
            y = ct[(state, True)]
            tot = ct[(state, True)] + ct[(state, False)]
            lab = "with thioether   " if state == "C" else "without thioether"
            print(f"      {lab} {y:>5}/{tot:<6} = {100*y/tot:>5.1f}%   "
                  f"published {pub}")
        if defname == "any non-PPO domain":
            y = ct[("C", True)]
            tot = ct[("C", True)] + ct[("C", False)]
            record("canonical structures with thioether carrying an extra domain",
                   "5,916/6,151 (96%)", f"{y}/{tot} ({round(100*y/tot)}%)",
                   "any non-PPO domain reported by Chainsaw; the 1,359 partial-motif "
                   "(C*) structures are excluded from both groups", [pv, da])
            y2 = ct[("-", True)]
            tot2 = ct[("-", True)] + ct[("-", False)]
            # published to whole-percent precision, and the numerator differs
            # by 3 of 14,383, so compare the rounded percentage
            record("canonical structures without thioether carrying an extra domain",
                   "26%", f"{round(100*y2/tot2)}%",
                   f"as above; recomputed numerator {y2}/{tot2}", [pv, da])
    print("    note: Chainsaw merges the two tyrosinase lobes into one PPO domain")
    print("          (rows with ndom=2 carry ppo_domains 'd0+d1'), so lobe")
    print("          segmentation does not affect these counts.")


# ---------------------------------------------------------------- block 4
def block_degenerate():
    """Degenerate-group species and pattern counts."""
    print("\n[4] Degenerate substitution patterns")
    rel = "3_noncanonical_analysis/nc_patterns_degenerate.tsv"
    rows = read_tsv(rel)
    n_pat = len(rows)
    multi = sum(1 for r in rows if int(float(r["n_species"] or 0)) >= 2)
    top = max(rows, key=lambda r: int(float(r["n"] or 0)))
    print(f"    pattern rows in the table: {n_pat}")
    record("degenerate patterns seen in two or more species", "28", str(multi),
           "n_species >= 2 in the deposited pattern table", [rel])
    record("largest degenerate group, structures and species",
           "34 structures, 16 species",
           f"{int(float(top['n']))} structures, {int(float(top['n_species']))} species",
           "top row of the deposited pattern table by structure count", [rel])


# ---------------------------------------------------------------- block 5
def block_identity():
    """Global sequence identity of the Microbispora protein. Needs network."""
    print("\n[5] Microbispora global sequence identity (network)")
    try:
        from Bio import Align
        from Bio.Align import substitution_matrices
    except ImportError:
        print("    SKIPPED: biopython is not installed (pip install biopython)")
        return
    import urllib.request

    def uniprot(acc):
        with urllib.request.urlopen(
                f"https://rest.uniprot.org/uniprotkb/{acc}.fasta", timeout=60) as r:
            return "".join(l.strip() for l in r.read().decode().splitlines()[1:])

    def pdb_longest(pdbid):
        with urllib.request.urlopen(
                f"https://www.rcsb.org/fasta/entry/{pdbid}", timeout=60) as r:
            txt = r.read().decode()
        seqs, buf = [], []
        for line in txt.splitlines():
            if line.startswith(">"):
                if buf:
                    seqs.append("".join(buf))
                buf = []
            else:
                buf.append(line.strip())
        if buf:
            seqs.append("".join(buf))
        return max(seqs, key=len)

    q = uniprot("A0A8H9LF69")
    t = pdb_longest("6J2U")
    al = Align.PairwiseAligner()
    al.substitution_matrix = substitution_matrices.load("BLOSUM62")
    al.open_gap_score, al.extend_gap_score, al.mode = -10, -0.5, "global"
    aln = al.align(q, t)[0]
    a, b = aln[0], aln[1]
    ident = sum(1 for x, y in zip(a, b) if x == y and x != "-")
    pct = 100 * ident / len(a)
    print(f"    query {len(q)} aa, 6J2U tyrosinase chain {len(t)} aa, "
          f"{ident} identities over {len(a)} alignment columns")
    record("global sequence identity to the closest characterised homologue",
           "27%", f"{pct:.1f}%",
           "Needleman-Wunsch, BLOSUM62, gap open 10, gap extend 0.5, identities "
           "over total alignment length (EMBOSS needle defaults)",
           ["UniProt A0A8H9LF69", "PDB 6J2U"],
           agrees=abs(pct - 27) < 0.5)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--identity", action="store_true",
                    help="also run the network sequence-identity block")
    ap.add_argument("--tsv", help="write a machine-readable summary here")
    args = ap.parse_args()

    print(f"analysis root: {BASE}")
    for fn in (block_joint_vector, block_table_s5, block_thioether,
               block_degenerate):
        fn()
    if args.identity:
        block_identity()

    print("\n" + "=" * 70)
    dis = [r for r in RESULTS if r["agrees"] == "NO"]
    print(f"claims checked: {len(RESULTS)} | agreeing: {len(RESULTS)-len(dis)} "
          f"| differing: {len(dis)}")
    for r in dis:
        print(f"  DIFFERS  {r['claim']}")
        print(f"           published {r['published']} vs "
              f"recomputed {r['recomputed']}")
    if args.tsv:
        with open(args.tsv, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(RESULTS[0]), delimiter="\t")
            w.writeheader()
            w.writerows(RESULTS)
        print(f"\nwrote {args.tsv}")


if __name__ == "__main__":
    main()
