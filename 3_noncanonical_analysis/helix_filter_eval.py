#!/usr/bin/env python3
"""Evaluate helix structural-integrity filter criteria from CACHED per-helix coverage
(no CEAlign re-run). Input = per-helix TSV from helix_ss.py with columns:
  accession, status, helix{1..4}_cov (backbone presence), helix{1..4}_hcov (alpha-helical), min_cov, min_hcov
PmTYR helix lengths: h1=6 (37-42), h2=11 (67-77), h3=5 (204-208), h4=10 (228-238).
EDIT the criteria below to tune; re-running is instant.
Usage: python3 helix_filter_eval.py <perhelix.tsv> [annotate]
  'annotate' cross-references characterised / candidate-group / showcase membership."""
import csv, sys, os, glob

SHORT = ['helix1','helix3']    # short helices: P-SEA unreliable -> require backbone PRESENCE
LONG  = ['helix2','helix4']    # long helices : reliable -> require ALPHA-HELICAL character
HLEN  = {'helix1':6,'helix2':11,'helix3':5,'helix4':10}
def f(r,h,k): return float(r[f'{h}_{k}'])

# ---- tunable criteria (add your own; each takes a row, returns True=KEEP) -------------
def strict(r):                                   # current filter: helical>=0.5 on ALL four
    return all(f(r,h,'hcov')>=0.5 for h in HLEN)
def hybrid(r):                                   # presence on short, helical on long
    return all(f(r,h,'cov')>=0.5 for h in SHORT) and all(f(r,h,'hcov')>=0.5 for h in LONG)
def hybrid_strictpres(r):                        # presence>=0.8 on short, helical>=0.5 on long
    return all(f(r,h,'cov')>=0.8 for h in SHORT) and all(f(r,h,'hcov')>=0.5 for h in LONG)
def minres3(r):                                  # >=3 helical residues in every helix
    return all(f(r,h,'hcov')*HLEN[h]>=3 for h in HLEN)
def plain(r):                                    # backbone present (>=0.5) in all four
    return all(f(r,h,'cov')>=0.5 for h in HLEN)
CRITERIA = {'strict_hcov0.5':strict,'hybrid':hybrid,'hybrid_pres0.8':hybrid_strictpres,
            'min3helical':minres3,'plain_cov0.5':plain}
# --------------------------------------------------------------------------------------

def load_annot():
    BASE="/cluster/work/projects/nn1003k/eirin/bioinf"
    CA=f"{BASE}/bioinf_redo/2_canonical_analysis"
    import openpyxl
    chars={}
    wb=openpyxl.load_workbook(f"{BASE}/characterized_PPOs.xlsx",data_only=True)
    for row in wb.active.iter_rows(min_row=2,values_only=True):
        if row[0]: chars[str(row[0]).strip()]=f"{row[1]}:{row[2]}"
    POS=['Gly46','Phe65','Trp68','Glu195','Asn205','Arg209','Val218','Ala221','Phe227','His230','thioether']
    flg={(r['position'],r['residue']) for r in csv.DictReader(open(f"{CA}/supplementary/flagged_groups.tsv"),delimiter='\t')}
    vec={}
    for r in csv.DictReader(open(f"{CA}/position_vectors.csv")):
        if r.get('vector'):
            p=[ (x.strip().rstrip('*') or '?') for x in r['vector'].split('-')]
            while len(p)<10: p.append('?')
            vec[r['accession']]=p
    def cand(a):
        if a in vec:
            for (pp,rr) in flg:
                if vec[a][POS.index(pp)]==rr: return f"{pp}={rr}"
        return ""
    show={os.path.basename(p).split('_taxID_')[0].split('_',1)[-1] for p in glob.glob(f"{CA}/../3_noncanonical_analysis/showcase_structures/*/*.cif")+glob.glob(f"{CA}/../3_noncanonical_analysis/microbispora/*.cif")}
    return chars,cand,show

def main():
    path=sys.argv[1]; annotate='annotate' in sys.argv[2:]
    rows=[r for r in csv.DictReader(open(path),delimiter='\t') if r['status']=='ok']
    print(f"{len(rows)} structures from {os.path.basename(path)}\n")
    print(f"{'criterion':18} {'KEEP':>7} {'DROP':>7}  recovers-vs-strict")
    base={r['accession'] for r in rows if strict(r)}
    for name,fn in CRITERIA.items():
        keep={r['accession'] for r in rows if fn(r)}
        rec=len(keep-base)
        print(f"  {name:16} {len(keep):>7} {len(rows)-len(keep):>7}  +{rec}")
    if annotate:
        chars,cand,show=load_annot()
        print("\nFalse-removal cost per criterion (important structures it DROPS):")
        for name,fn in CRITERIA.items():
            drop=[r['accession'] for r in rows if not fn(r)]
            nc=sum(a in chars for a in drop); ncand=sum(bool(cand(a)) for a in drop); ns=sum(a in show for a in drop)
            print(f"  {name:16} characterised:{nc:>3}  candidate-members:{ncand:>4}  showcase:{ns:>2}")

if __name__=='__main__': main()
