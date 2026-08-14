"""
QC for the Phe227=V and Phe227=W candidate groups (and reference 119 / baseline).
Two questions:
  (1) Is the region well-modelled?  -> active_site_plddt, pmtyr_rmsd, hmm_coverage
  (2) Is the substitution real (not an alignment artifact)?
        -> CA-mapping distance + secondary-structure-context match at the
           substituted position(s)
"""
import csv
import numpy as np
from collections import Counter

POS = ['Gly46','Phe65','Trp68','Glu195','Asn205','Arg209','Val218','Ala221','Phe227','His230','thioether']

def norm(r):
    r = str(r).strip().rstrip('*'); return '~' if r in ('','~','None') else r
def vfrom(v):
    p=[norm(x) for x in v.split('-')]
    while len(p)<11: p.append('~')
    return p[:11]

# filter_results: active_site_plddt, hmm_coverage
fr = {}
for r in csv.DictReader(open('../1_filtering/filter_results.csv')):
    acc = r['sequence_id'].split('_taxID_')[0]
    fr[acc] = r
# position_vectors: alignment + per-position cadist/ss_match
pv = {}
for r in csv.DictReader(open('position_vectors.csv')):
    pv[r['accession']] = r

def fnum(d, k):
    try: return float(d[k])
    except (KeyError, ValueError, TypeError): return np.nan

# define groups
allacc = [a for a in pv if pv[a].get('vector')]
def vec(a): return vfrom(pv[a]['vector'])
groups = {
  'Phe227=V (cand)':      [a for a in allacc if vec(a)[8]=='V'],
  'Phe227=W (cand)':      [a for a in allacc if vec(a)[8]=='W'],
  '119 His230=Y (ref)':   [a for a in allacc if vec(a)[4]=='D' and vec(a)[5]=='G' and vec(a)[6]=='N' and vec(a)[9]=='Y'],
  'whole dataset':        allacc,
}

def summ(vals):
    v=np.array([x for x in vals if not np.isnan(x)])
    return (np.median(v), np.percentile(v,10), v.min(), v.max(), len(v)) if len(v) else (np.nan,)*4+(0,)

out=open('qc_candidate_groups.txt','w')
def w(s): out.write(s+'\n'); print(s)

for gname, accs in groups.items():
    w('='*72); w(f'{gname}   n={len(accs)}'); w('='*72)
    plddt=[fnum(fr[a],'active_site_plddt') for a in accs if a in fr]
    hmm  =[fnum(fr[a],'hmm_coverage_pct') for a in accs if a in fr]
    rmsd =[fnum(pv[a],'pmtyr_rmsd') for a in accs]
    cad  =[fnum(pv[a],'Phe227_cadist') for a in accs]
    ssm  =[pv[a].get('Phe227_ss_match') for a in accs]
    for lbl,vals in [('active_site_pLDDT',plddt),('hmm_coverage_%',hmm),
                     ('pmtyr_rmsd(Å)',rmsd),('Phe227_CA_dist(Å)',cad)]:
        med,p10,lo,hi,n = summ(vals)
        w(f'  {lbl:20s} median={med:6.2f}  10th_pct={p10:6.2f}  range=[{lo:.2f},{hi:.2f}]  n={n}')
    ss_true = sum(1 for x in ssm if x=='True'); ss_n=sum(1 for x in ssm if x in ('True','False'))
    w(f'  Phe227_ss_match      {ss_true}/{ss_n} True ({100*ss_true/ss_n if ss_n else 0:.0f}%)')
    if gname=='Phe227=W (cand)':  # double substitution: also QC His230=L mapping
        h_cad=[fnum(pv[a],'His230_cadist') for a in accs]
        h_ss=[pv[a].get('His230_ss_match') for a in accs]
        med,p10,lo,hi,n=summ(h_cad)
        w(f'  His230_CA_dist(Å)    median={med:.2f}  10th_pct={p10:.2f}  range=[{lo:.2f},{hi:.2f}]')
        w(f'  His230_ss_match      {sum(1 for x in h_ss if x=="True")}/{len(h_ss)} True')
    # flag bad structures
    bad=[a for a in accs if (a in fr and fnum(fr[a],'active_site_plddt')<70)
         or fnum(pv[a],'Phe227_cadist')>2.0 or fnum(pv[a],'pmtyr_rmsd')>3.0
         or (a in fr and fnum(fr[a],'hmm_coverage_pct')<70)]
    w(f'  FLAGGED (pLDDT<70 | CAdist>2 | rmsd>3 | hmm<70): {len(bad)}/{len(accs)}')
    if bad and gname.endswith('(cand)'):
        for a in bad[:20]:
            w(f'     {a}  plddt={fnum(fr.get(a,{}),"active_site_plddt"):.0f} '
              f'cad={fnum(pv[a],"Phe227_cadist"):.2f} rmsd={fnum(pv[a],"pmtyr_rmsd"):.2f} '
              f'hmm={fnum(fr.get(a,{}),"hmm_coverage_pct"):.0f}')
    w('')
out.close()
print('\nwrote qc_candidate_groups.txt')
