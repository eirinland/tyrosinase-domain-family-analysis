"""Consolidated summary table for the novel-signature candidate groups.
Encoding convention (kept semantically correct):
  '~' = loop at Gly46 only -> wildcard (skipped in Hamming).
  '-' = absence of the thioether cysteine -> a real state; two '-' match,
        '-' vs 'C' counts as a difference. Thioether is read from its own
        column (C/C*/-), not parsed off the vector-string tail.
  '?' = unmapped/unknown -> skipped.
"""
import csv, statistics, openpyxl
from collections import Counter

POS=['Gly46','Phe65','Trp68','Glu195','Asn205','Arg209','Val218','Ala221','Phe227','His230','thioether']
def res(r):                       # positions 0-9
    if r is None: return '?'
    r=str(r).strip().rstrip('*'); return '?' if r in ('','None') else r
def thio(r):                      # thioether column
    r=str(r).strip()
    return 'C' if r in ('C','C*') else ('-' if r=='-' else '?')
def hd(a,b):
    d=0
    for i,(x,y) in enumerate(zip(a,b)):
        if i==0 and (x=='~' or y=='~'): continue   # Gly46 loop wildcard
        if x=='?' or y=='?': continue               # unmapped
        if x!=y: d+=1
    return d
def fnum(d,k):
    try: return float(d[k])
    except: return None

# characterized (xlsx columns; Gly46 keeps '~', thioether from column)
wb=openpyxl.load_workbook("/cluster/work/projects/nn1003k/eirin/bioinf/characterized_PPOs.xlsx",data_only=True)
ws=wb["Characterized PPOs"]; R=list(ws.iter_rows(values_only=True)); Hh=list(R[0])
ci={p:Hh.index(p) for p in POS}; ai,acti=Hh.index('Accession'),Hh.index('Activity')
charv=[(r[ai],r[acti],[res(r[ci[p]]) for p in POS[:10]]+[thio(r[ci['thioether']])]) for r in R[1:] if r[0]]

tax={r['accession']:(r['kingdom'],r['phylum'],r['genus']) for r in csv.DictReader(open("visualisation/taxonomy_lookup.csv"))}
fr={r['sequence_id'].split('_taxID_')[0]:r for r in csv.DictReader(open('../1_filtering/filter_results.csv'))}
pvr={r['accession']:r for r in csv.DictReader(open('position_vectors.csv')) if r.get('vector')}
V={}
for a,r in pvr.items():
    parts=[res(x) for x in r['vector'].split('-')]
    while len(parts)<10: parts.append('?')
    V[a]=parts[:10]+[thio(r['thioether'])]

def members(pred): return [a for a in V if pred(V[a])]
GROUPS=[
 ("His230=Tyr core (the 119)", "Asn205=D, Arg209=G, Val218=N, His230=Y",
   lambda v: v[4]=='D' and v[5]=='G' and v[6]=='N' and v[9]=='Y', 'His230', 9, 'Y'),
 ("Asn205=Lys", "Asn205=K (cationic residue at the controller pair)",
   lambda v: v[4]=='K', 'Asn205', 4, 'K'),
 ("Phe227=Trp / His230=Leu (black yeasts)", "Phe227=W + His230=L (+Glu195=L in consensus)",
   lambda v: v[8]=='W' and v[9]=='L', 'Phe227', 8, 'W'),
]
def consensus(accs):
    out=[]
    for j in range(11):
        c=Counter(V[a][j] for a in accs if V[a][j]!='?' and not (j==0 and V[a][j]=='~'))
        out.append(c.most_common(1)[0][0] if c else '~')
    return out
def med(xs):
    xs=[x for x in xs if x is not None]; return statistics.median(xs) if xs else float('nan')
def nearest(v):
    b=min(charv,key=lambda c:hd(v,c[2])); return hd(v,b[2]),b
def char_with(idx,r): return [c[0] for c in charv if c[2][idx]==r]

rows=[]
for name,defn,pred,qpos,qidx,subres in GROUPS:
    accs=members(pred); n=len(accs)
    hds=[]; nearc=Counter()
    for a in accs:
        h,b=nearest(V[a]); hds.append(h); nearc[f'{b[0]} ({b[1]})']+=1
    cw=char_with(qidx,subres)
    kc=Counter(tax.get(a,('?','?','?'))[0] for a in accs)
    gc=Counter(tax.get(a,('?','?','?'))[2] for a in accs)
    bigk=[k for k,_ in kc.most_common() if k!='?' and kc[k]>=5][:2]
    conv='n/a (single lineage)'
    if len(bigk)==2:
        c1=consensus([a for a in accs if tax.get(a,('?','?','?'))[0]==bigk[0]])
        c2=consensus([a for a in accs if tax.get(a,('?','?','?'))[0]==bigk[1]])
        conv=f"{bigk[0]} vs {bigk[1]} HD={hd(c1,c2)}"
    plddt=med([fnum(fr[a],'active_site_plddt') for a in accs if a in fr])
    hmm=med([fnum(fr[a],'hmm_coverage_pct') for a in accs if a in fr])
    rmsd=med([fnum(pvr[a],'pmtyr_rmsd') for a in accs])
    cad=med([fnum(pvr[a],f'{qpos}_cadist') for a in accs])
    ssm=[pvr[a].get(f'{qpos}_ss_match') for a in accs]
    ssp=100*sum(1 for x in ssm if x=='True')/len(ssm)
    flagged=sum(1 for a in accs if (a in fr and (fnum(fr[a],'active_site_plddt') or 0)<70)
                or (fnum(pvr[a],f'{qpos}_cadist') or 0)>2.0 or (fnum(pvr[a],'pmtyr_rmsd') or 0)>3.0
                or (a in fr and (fnum(fr[a],'hmm_coverage_pct') or 0)<70))
    rows.append(dict(name=name,defn=defn,cons='-'.join(consensus(accs)),n=n,
        toptax='; '.join(f'{g} {c}' for g,c in gc.most_common(4) if g!='?'),
        conv=conv, medhd=f'{med(hds):.0f} ({min(hds)}-{max(hds)})',
        pct4=f'{100*sum(1 for h in hds if h>=4)/n:.0f}',
        nearest=nearc.most_common(1)[0][0],
        subst='absent (0/84)' if not cw else f'{len(cw)}/84',
        plddt=f'{plddt:.0f}',hmm=f'{hmm:.0f}',rmsd=f'{rmsd:.2f}',cad=f'{cad:.2f}',
        ssp=f'{ssp:.0f}',flagged=f'{flagged}/{n}'))

cols=[('Group','name'),('Defining signature','defn'),('Consensus vector','cons'),('n','n'),
      ('Top genera (n)','toptax'),('Cross-kingdom convergence','conv'),
      ('Median HD (range)','medhd'),('% HD≥4','pct4'),('Nearest charact.','nearest'),
      ('Defining subst. in 84 charact.','subst'),
      ('Active-site pLDDT (med)','plddt'),('HMM cov % (med)','hmm'),('Align RMSD Å (med)','rmsd'),
      ('Subst.-pos CA-dist Å (med)','cad'),('SS-match %','ssp'),('Flagged','flagged')]
with open('candidate_groups_summary.tsv','w',newline='') as f:
    w=csv.writer(f,delimiter='\t'); w.writerow([h for h,_ in cols])
    for r in rows: w.writerow([r[k] for _,k in cols])
with open('candidate_groups_summary.md','w') as f:
    f.write('| '+' | '.join(h for h,_ in cols)+' |\n|'+'|'.join('---' for _ in cols)+'|\n')
    for r in rows: f.write('| '+' | '.join(str(r[k]) for _,k in cols)+' |\n')
print(open('candidate_groups_summary.md').read())
