"""
Find coherent, novel active-site signature groups in the canonical PPO set,
analogous to the 119-structure His230=Y / Fusarium+Pseudomonas convergent group.

Interesting = (a) novel (Hamming >=4 from every characterized PPO vector),
              (b) coherent (a tight shared signature),
              (c) taxonomically striking (cross-kingdom/phylum convergence,
                  or strong genus concentration),
              (d) substitution at a functionally important / conserved position.

Inputs: position_vectors.csv, visualisation/taxonomy_lookup.csv, characterized_PPOs.xlsx
"""
import csv, openpyxl
import numpy as np
from collections import Counter, defaultdict

POS = ['Gly46','Phe65','Trp68','Glu195','Asn205','Arg209','Val218','Ala221','Phe227','His230','thioether']
CANON = {'Gly46':'G','Phe65':'F','Trp68':'W','Glu195':'E','Asn205':'N','Arg209':'R',
         'Val218':'V','Ala221':'A','Phe227':'F','His230':'H','thioether':'C'}
HD_NOVEL = 4

def norm(r):
    if r is None: return '?'
    r = str(r).strip().rstrip('*')
    return '~' if r in ('', '~') else r

def vec_from_string(v):
    parts = v.split('-')
    # thioether is last; empty -> '-'
    parts = [norm(p) for p in parts]
    while len(parts) < 11: parts.append('~')
    return parts[:11]

def hd(a, b):
    """Hamming. Only Gly46 '~' (loop) is a wildcard, matching the manuscript
    convention; '?' (unmapped) is skipped everywhere; '~' at any other
    position is a real (loop) state that differs from a defined residue."""
    d = 0
    for i, (x, y) in enumerate(zip(a, b)):
        if i == 0 and (x == '~' or y == '~'): continue   # Gly46 loop wildcard
        if x == '?' or y == '?': continue                # unmapped
        if x != y: d += 1
    return d

# --- characterized vectors ---
wb = openpyxl.load_workbook("/cluster/work/projects/nn1003k/eirin/bioinf/characterized_PPOs.xlsx", data_only=True)
ws = wb["Characterized PPOs"]; R = list(ws.iter_rows(values_only=True)); Hh = list(R[0])
ci = {p: Hh.index(p) for p in POS}
char_vecs = []
for r in R[1:]:
    if not r[0]: continue
    char_vecs.append([norm(r[ci[p]]) for p in POS])

# --- taxonomy ---
tax = {}
for r in csv.DictReader(open("visualisation/taxonomy_lookup.csv")):
    tax[r['accession']] = (r['kingdom'], r['phylum'], r['genus'])

# --- per-structure vectors ---
structs = []  # (acc, vec, minHD, kingdom, phylum, genus)
for r in csv.DictReader(open("position_vectors.csv")):
    acc = r['accession']
    vstr = r.get('vector','')
    if not vstr: continue
    v = vec_from_string(vstr)
    mh = min(hd(v, cv) for cv in char_vecs)
    k,p,g = tax.get(acc, ('?','?','?'))
    structs.append((acc, tuple(v), mh, k, p, g))
print(f"loaded {len(structs)} structures; {sum(1 for s in structs if s[2]>=HD_NOVEL)} novel (HD>={HD_NOVEL})")

# --- position conservation (whole set) ---
print("\n=== POSITION CONSERVATION (canonical residue %, top variants) ===")
for i,p in enumerate(POS):
    c = Counter(s[1][i] for s in structs if s[1][i] not in '~?')
    tot = sum(c.values())
    canon_pct = 100*c.get(CANON[p],0)/tot if tot else 0
    variants = ', '.join(f'{a}:{n}' for a,n in c.most_common(5) if a!=CANON[p])
    print(f"  {p:10s} {CANON[p]}={canon_pct:5.1f}%   variants: {variants}")

def profile(members, label):
    n = len(members)
    ks = Counter(m[3] for m in members); ps = Counter(m[4] for m in members); gs = Counter(m[5] for m in members)
    medhd = np.median([m[2] for m in members])
    nov = 100*sum(1 for m in members if m[2]>=HD_NOVEL)/n
    # convergence: >=2 kingdoms each >=15%, or fungi+bacteria both present
    kdist = {k: 100*v/n for k,v in ks.items() if k!='?'}
    fb = ('Fungi' in kdist and 'Bacteria' in kdist and kdist.get('Fungi',0)>=15 and kdist.get('Bacteria',0)>=15)
    print(f"\n  [{label}] n={n}  medianHD={medhd:.0f}  novel%={nov:.0f}  genera={len([x for x in gs if x!='?'])}")
    print(f"     kingdoms: {dict(ks.most_common(4))}")
    print(f"     top genera: {dict(gs.most_common(5))}")
    if fb: print(f"     ** CROSS-KINGDOM convergence (Fungi {kdist['Fungi']:.0f}% + Bacteria {kdist['Bacteria']:.0f}%) **")
    return dict(n=n, medhd=medhd, nov=nov, fb=fb, ks=ks, gs=gs)

# --- SCAN 1: novel full-vector groups with cross-kingdom convergence ---
print("\n" + "="*72)
print("SCAN 1: novel full-vector groups (n>=25, HD>=4), ranked by cross-kingdom mix")
print("="*72)
byvec = defaultdict(list)
for s in structs:
    if s[2] >= HD_NOVEL: byvec[s[1]].append(s)
cands = []
for v, mem in byvec.items():
    if len(mem) < 25: continue
    ks = Counter(m[3] for m in mem if m[3]!='?'); tot=sum(ks.values())
    if tot==0: continue
    fungi=100*ks.get('Fungi',0)/tot; bact=100*ks.get('Bacteria',0)/tot
    mix = min(fungi,bact)  # high when both present
    cands.append((mix, len(mem), v, mem, fungi, bact))
cands.sort(reverse=True)
for mix,n,v,mem,fu,ba in cands[:8]:
    profile(mem, '-'.join(v))

# --- SCAN 2: controller-core convergence (reproduces the 119, finds analogues) ---
print("\n" + "="*72)
print("SCAN 2: novel CONTROLLER-CORE groups (Asn205,Arg209,Val218,Ala221,His230)")
print("        n>=25, ranked by cross-kingdom (min(Fungi%,Bacteria%)))")
print("="*72)
CORE_IDX = [POS.index(p) for p in ['Asn205','Arg209','Val218','Ala221','His230']]
bycore = defaultdict(list)
for s in structs:
    if s[2] >= HD_NOVEL:
        bycore[tuple(s[1][i] for i in CORE_IDX)].append(s)
core_cands = []
for core, mem in bycore.items():
    if len(mem) < 25: continue
    ks = Counter(m[3] for m in mem if m[3]!='?'); tot=sum(ks.values())
    if not tot: continue
    fungi=100*ks.get('Fungi',0)/tot; bact=100*ks.get('Bacteria',0)/tot
    oom=100*ks.get('Oomycota',0)/tot; ani=100*ks.get('Animalia',0)/tot
    # convergence = high min across the two biggest distinct-lineage kingdoms
    top2 = sorted([fungi,bact,oom,ani], reverse=True)[:2]
    core_cands.append((top2[1], len(mem), core, mem))
core_cands.sort(reverse=True)
for conv,n,core,mem in core_cands[:10]:
    lbl='Asn205='+core[0]+' Arg209='+core[1]+' Val218='+core[2]+' Ala221='+core[3]+' His230='+core[4]
    profile(mem, lbl)

# --- SCAN 3: rare substitutions at CONSERVED positions (anchor & profile) ---
print("\n" + "="*72)
print("SCAN 3: rare substitutions at conserved positions (all carriers profiled)")
print("="*72)
ANCHORS = [('Phe227','W'),('Phe227','V'),('Phe227','Y'),
           ('Trp68','Y'),('Trp68','A'),
           ('His230','N'),('His230','I'),('His230','F'),
           ('Phe65','T'),('Phe65','M'),('Glu195','D')]
for pos,res in ANCHORS:
    i = POS.index(pos)
    mem = [s for s in structs if s[1][i]==res]
    if len(mem) < 20: 
        print(f"\n  {pos}={res}: n={len(mem)} (too few)"); continue
    # consensus signature
    cons = []
    for j,p in enumerate(POS):
        c = Counter(m[1][j] for m in mem if m[1][j] not in '~?')
        cons.append(c.most_common(1)[0][0] if c else '~')
    profile(mem, f"{pos}={res} | consensus " + '-'.join(cons))
