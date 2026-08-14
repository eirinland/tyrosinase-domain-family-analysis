"""
Cross-reference characterized PPOs against the two fungal C-terminal domain
types found inside TM0.8 cluster 3 (rep A0AAJ0MJ41).

types = the two largest TM0.7 sub-clusters of the extracted C-terminal domains
(short ~174 res vs long ~421 res). For each characterized PPO that sits in TM0.8
cluster 3, report which type its C-terminal domain belongs to.
"""
import csv, re, sys
from collections import defaultdict
import openpyxl

PIPE = "/cluster/work/projects/nn1003k/eirin/bioinf/Super_reference_pipeline"
CTERM = f"{PIPE}/foldseek/pools/cterm_domain_similarity/c3_fungi"
CLUST08 = f"{PIPE}/foldseek/pools/results/cluster_cluster.tsv"
SUB07 = f"{CTERM}/sub07/clu_cluster.tsv"
MAN = f"{CTERM}/manifest.csv"
TAX = f"{PIPE}/taxonomy_lookup.csv"
XLSX = "/cluster/work/projects/nn1003k/eirin/bioinf/characterized_PPOs.xlsx"
REP = "A0AAJ0MJ41"


def bare(name):
    return name.split("_taxID_")[0]


def strip_a(name):
    return name[:-2] if name.endswith("_A") else name


# --- TM0.8 cluster 3 membership (bare acc) ---
clust3 = set()
with open(CLUST08) as f:
    for line in f:
        r, m = line.rstrip("\n").split("\t")
        if bare(r) == REP:
            clust3.add(bare(m))

# --- manifest: bare acc -> (full member, cterm_len, status) ---
man = {}
with open(MAN) as f:
    for r in csv.DictReader(f):
        man[r["accession"]] = (r["member"], r["cterm_len"], r["status"])

# --- TM0.7 sub-cluster of the C-terminal domains: bare acc -> sub rep ---
sub_of = {}
members_of = defaultdict(list)
with open(SUB07) as f:
    for line in f:
        r, m = line.rstrip("\n").split("\t")
        rr, mm = strip_a(r), strip_a(m)
        sub_of[mm] = rr
        members_of[rr].append(mm)

# label the two largest sub-clusters by size -> short/long by median length
def med_len(accs):
    L = sorted(int(man[a][1]) for a in accs if a in man and man[a][1].isdigit())
    return L[len(L) // 2] if L else -1

top = sorted(members_of.items(), key=lambda kv: len(kv[1]), reverse=True)[:2]
top_sorted = sorted(top, key=lambda kv: med_len(kv[1]))
SHORT_REP, LONG_REP = top_sorted[0][0], top_sorted[1][0]
label = {SHORT_REP: "SHORT", LONG_REP: "LONG"}

# --- taxonomy ---
tax = {}
with open(TAX) as f:
    for r in csv.DictReader(f):
        tax[r["accession"]] = (r.get("phylum", "?"), r.get("genus", "?"))

# --- characterized PPOs ---
wb = openpyxl.load_workbook(XLSX)
ws = wb["Characterized PPOs"]
chars = {}
for row in list(ws.iter_rows(values_only=True))[1:]:
    if row and row[0]:
        chars[str(row[0]).strip()] = (row[1], row[2])  # activity, uniprot name


def typ(acc):
    if acc not in man or man[acc][2] != "ok":
        return "no_Cterm_domain"
    s = sub_of.get(acc)
    if s in label:
        return label[s]
    return f"minor_sub({med_len(members_of.get(s, []))}res)"


print(f"TM0.8 cluster 3 (rep {REP}): {len(clust3)} members")
print(f"  SHORT type = TM0.7 sub-rep {SHORT_REP}: "
      f"{len(members_of[SHORT_REP])} domains, median {med_len(members_of[SHORT_REP])} res")
print(f"  LONG  type = TM0.7 sub-rep {LONG_REP}: "
      f"{len(members_of[LONG_REP])} domains, median {med_len(members_of[LONG_REP])} res")
print()

# representative example structures (rep + nearest-to-median members)
def examples(rep_acc, n=4):
    accs = [a for a in members_of[rep_acc] if a in man]
    md = med_len(members_of[rep_acc])
    accs.sort(key=lambda a: abs(int(man[a][1]) - md) if man[a][1].isdigit() else 9999)
    out = []
    for a in [rep_acc] + [x for x in accs if x != rep_acc]:
        if a not in man:
            continue
        ph, ge = tax.get(a, ("?", "?"))
        out.append(f"    {man[a][0]}  len={man[a][1]}  {ge} ({ph})")
        if len(out) >= n:
            break
    return out

print("SHORT-type example structures:")
print("\n".join(examples(SHORT_REP)))
print("\nLONG-type example structures:")
print("\n".join(examples(LONG_REP)))
print()

# characterized overlap
in3 = [a for a in chars if a in clust3]
print(f"Characterized PPOs total: {len(chars)} | in TM0.8 cluster 3: {len(in3)}")
print()
print(f"{'acc':<12}{'type':<18}{'len':>5}  {'activity':<6}{'genus/phylum':<26}uniprot_name")
bytype = defaultdict(list)
for a in sorted(in3, key=lambda x: typ(x)):
    t = typ(a)
    bytype[t].append(a)
    act, nm = chars[a]
    ph, ge = tax.get(a, ("?", "?"))
    clen = man[a][1] if a in man else "-"
    print(f"{a:<12}{t:<18}{str(clen):>5}  {str(act):<6}{(ge+'/'+ph):<26}{nm}")

print()
print("counts by type (characterized in cluster 3):")
for t, v in sorted(bytype.items(), key=lambda kv: -len(kv[1])):
    print(f"  {t:<20}: {len(v)}  {sorted(v)}")
