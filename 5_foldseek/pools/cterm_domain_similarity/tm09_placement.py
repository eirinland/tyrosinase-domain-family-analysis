"""
Where do the two fungal C-terminal domain types (defined inside TM0.8 cluster 3)
land in the whole-structure TM0.9 network clustering? Do they separate?
"""
import csv
from collections import defaultdict, Counter

P = "/cluster/work/projects/nn1003k/eirin/bioinf/Super_reference_pipeline/foldseek/pools"
TM08 = f"{P}/results/cluster_cluster.tsv"
TM09 = f"{P}/tm09/results/cluster_cluster.tsv"
SUB07 = f"{P}/cterm_domain_similarity/c3_fungi/sub07/clu_cluster.tsv"
MAN = f"{P}/cterm_domain_similarity/c3_fungi/manifest.csv"
REP = "A0AAJ0MJ41"


def bare(n):
    return n.split("_taxID_")[0]


def strip_a(n):
    return n[:-2] if n.endswith("_A") else n


# TM0.8 cluster-3 membership
clust3 = set()
with open(TM08) as f:
    for line in f:
        r, m = line.rstrip("\n").split("\t")
        if bare(r) == REP:
            clust3.add(bare(m))

# manifest: bare acc -> (cterm_len, status)
man = {}
with open(MAN) as f:
    for r in csv.DictReader(f):
        man[r["accession"]] = (r["cterm_len"], r["status"])

# C-terminal domain TM0.7 sub-cluster -> short/long label
members_of = defaultdict(list)
sub_of = {}
with open(SUB07) as f:
    for line in f:
        r, m = line.rstrip("\n").split("\t")
        members_of[strip_a(r)].append(strip_a(m))
        sub_of[strip_a(m)] = strip_a(r)

def med_len(accs):
    L = sorted(int(man[a][0]) for a in accs if a in man and man[a][0].isdigit())
    return L[len(L) // 2] if L else -1

top2 = sorted(members_of.items(), key=lambda kv: len(kv[1]), reverse=True)[:2]
top2 = sorted(top2, key=lambda kv: med_len(kv[1]))
SHORT, LONG = top2[0][0], top2[1][0]
label = {SHORT: "SHORT", LONG: "LONG"}

def typ(a):
    if a not in man or man[a][1] != "ok":
        return "no_Cterm"
    return label.get(sub_of.get(a), "minor")

# whole-structure TM0.9: bare acc -> tm09 rep
tm09rep = {}
tm09size = Counter()
with open(TM09) as f:
    for line in f:
        r, m = line.rstrip("\n").split("\t")
        tm09rep[bare(m)] = bare(r)
        tm09size[bare(r)] += 1

# map cluster-3 members by type into tm09 clusters
by_t = defaultdict(Counter)          # type -> Counter(tm09rep)
for a in clust3:
    by_t[typ(a)][tm09rep.get(a, "UNCLUSTERED")] += 1

print(f"TM0.8 cluster 3 = {len(clust3)} members "
      f"(SHORT median {med_len(members_of[SHORT])}res, LONG median {med_len(members_of[LONG])}res)")
print()
for t in ("SHORT", "LONG", "minor", "no_Cterm"):
    c = by_t[t]
    tot = sum(c.values())
    if not tot:
        continue
    print(f"== {t}  (n={tot}) -> {len(c)} distinct TM0.9 whole-structure clusters ==")
    for rep, k in c.most_common(6):
        print(f"   tm09 rep {rep:<12} holds {k:>4} of this type "
              f"(rep's full tm09 cluster size = {tm09size[rep]})")
    print()

# separation check: do SHORT and LONG share any TM0.9 cluster?
short_cls = set(by_t["SHORT"])
long_cls = set(by_t["LONG"])
shared = short_cls & long_cls
print(f"SHORT occupies {len(short_cls)} TM0.9 clusters; LONG occupies {len(long_cls)}.")
print(f"TM0.9 clusters containing BOTH SHORT and LONG: {len(shared)} {sorted(shared) if shared else ''}")

# purity of the main long cluster(s)
print("\nType composition of the TM0.9 clusters that hold the LONG type:")
for rep, _ in by_t["LONG"].most_common(5):
    comp = Counter(typ(a) for a in clust3 if tm09rep.get(a) == rep)
    print(f"   tm09 {rep:<12} (size {tm09size[rep]}): {dict(comp)}")
