#!/usr/bin/env python3
"""Consolidated, reproducible pipeline for the PPO active-site novelty analysis.

Produces, from one entry point, every result that defines the groups flagged for
characterisation, with the FINAL agreed criteria:

  * placement (Ca-distance to the canonical position) is the disqualifying filter;
  * coherence and conservation are DESCRIPTORS for prioritisation, not gates;
  * no "erosion" language -- a displaced residue is "not confidently placed", not degraded;
  * Glu195 excluded (modelled imprecisely); Gly46 loops wildcarded; the thioether
    cysteine that maps to the Gly46 helix excluded (copper-verified, not a novelty);
  * the one displaced group kept (oomycete Gly46=Glu) is justified by COPPER distance,
    applied uniformly to all displaced candidate groups (stage E), not by a heuristic.

Stages
  A  coverage ladder (descriptive scale)
  B  activity inseparability among the characterised (gate Phe, controller pair, 96% ceiling)
  C  o-aminophenol-oxidase set (rarest signature; diagnostic Asn on diverse backgrounds)
  D  novelty enumeration + classification -> supplementary/novelty_enumeration.tsv
  E  copper-distance test of displaced candidate groups + thioether-Cys   [needs --cifs]
  F  single-substitution (one step from a characterised active site) leads
  G  HMM cross-validation of the flagged groups              [needs --afa or --hmm-cols]
  H  aromatic-ring compensation (superpose on six His; probe from reference ring)  [needs --cifs]
  I  HMM per-position agreement (structure mapping vs profile column; Table SX)
                                                             [needs --afa or --hmm-cols]
  J  within-group geometric consistency (Ca RMSD of the 10 variable positions)      [needs --cifs]
  K  confidence at the mapped positions (per-position pLDDT, random sample)          [needs --cifs]
  -> supplementary/flagged_groups.tsv, supplementary/hmm_agreement.tsv

Run everything via run_novelty_pipeline.sh (mounts the squashfs, supplies --cifs/--afa).
Without --cifs/--afa the vector-level stages (A-D, F) still run fully.

The 782 MB all_hmmalign.afa that --afa wants is a build artifact and is not
deposited. Regenerate it, or better the 2 MB match-state table that carries the
same information, with hmm/build_alignment.py, then pass
--hmm-cols hmm/hmm_match_columns.tsv.gz. Both routes give identical stage G/I output.
"""
import argparse, csv, os, math, statistics, sys
from itertools import combinations
from collections import Counter, defaultdict
import openpyxl

# ----------------------------------------------------------------------------- config
ROOT = os.environ.get("PPO_ROOT", "/cluster/work/projects/nn1003k/eirin/bioinf")


def _find_xlsx():
    """characterized_PPOs.xlsx lives one level above the pipeline on the cluster and
    at the repository root in the deposit; check both (and $PPO_XLSX) so the script
    runs from either without editing."""
    here = os.path.dirname(os.path.abspath(__file__))
    for cand in (os.environ.get("PPO_XLSX"),
                 os.path.join(here, "characterized_PPOs.xlsx"),
                 os.path.join(here, os.pardir, "characterized_PPOs.xlsx"),
                 os.path.join(here, os.pardir, os.pardir, "characterized_PPOs.xlsx"),
                 f"{ROOT}/characterized_PPOs.xlsx"):
        if cand and os.path.exists(cand):
            return os.path.normpath(cand)
    sys.exit("characterized_PPOs.xlsx not found; set $PPO_XLSX to its path")


XLSX = _find_xlsx()
PVEC = os.environ.get("PVEC", "position_vectors.csv")   # override to re-run with extra structures
TAXf = "visualisation/taxonomy_lookup.csv"
POS  = ['Gly46','Phe65','Trp68','Glu195','Asn205','Arg209','Val218','Ala221','Phe227','His230','thioether']
PMTYR_NUM = {'Gly46':46,'Phe65':65,'Trp68':68,'Glu195':195,'Asn205':205,'Arg209':209,
             'Val218':218,'Ala221':221,'Phe227':227,'His230':230}   # reference residue numbers
ENV  = 1.2          # placement envelope (A); calibrated vs correctly-placed residues
MINN = 20           # minimum carriers
ENVELOPE_POS = None                           # (Glu195 no longer excluded; judged on its own placement baseline)
THIOETHER_CYS = ('Gly46','C')                 # the crosslink cysteine on the Gly46 helix (Cu-verified)
CU_POCKET = 9.0     # A: a residue side chain within this of a copper is "in the pocket"
SKIP = {'?','~'}
AA3 = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H',
       'ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
FUNC = {'C':['SG'],'E':['OE1','OE2'],'D':['OD1','OD2'],'Y':['OH'],'K':['NZ'],'R':['NH1','NH2','NE'],
        'H':['ND1','NE2'],'W':['NE1'],'F':['CZ'],'S':['OG'],'T':['OG1'],'N':['OD1','ND2'],'Q':['OE1','NE2']}
# stage H (aromatic-ring compensation). His/imidazole is NOT treated as aromatic here.
RING_ATOMS = {'PHE':['CG','CD1','CD2','CE1','CE2','CZ'],'TYR':['CG','CD1','CD2','CE1','CE2','CZ'],
              'TRP':['CG','CD1','CD2','NE1','CE2','CE3','CZ2','CZ3','CH2']}
AROM     = {'F','Y','W'}                       # true aromatic one-letter codes (His excluded)
AROM_POS = ['Phe65','Trp68','Phe227']          # the vector's aromatic positions (His230 is not aromatic)
RING_DIST_TOL = 2.0 # A: a query aromatic ring centroid this close to the reference ring centroid
                    # (after six-His superposition) counts as occupying the reference site

# ----------------------------------------------------------------------------- helpers
def res(x):
    if x is None: return '?'
    x = str(x).strip().rstrip('*'); return '?' if x in ('','None') else x
def thio(x):
    x = str(x).strip(); return 'C' if x in ('C','C*') else ('-' if x=='-' else '?')
def hd(a,b):
    d=0
    for i,(x,y) in enumerate(zip(a,b)):
        if i==0 and (x=='~' or y=='~'): continue
        if x=='?' or y=='?': continue
        if x!=y: d+=1
    return d
def diffpos(a,b):
    out=[]
    for i,(x,y) in enumerate(zip(a,b)):
        if i==0 and (x=='~' or y=='~'): continue
        if x=='?' or y=='?': continue
        if x!=y: out.append(i)
    return out
def med(xs):
    xs=[x for x in xs if x is not None]; return statistics.median(xs) if xs else float('nan')
def pct(xs,q):
    xs=sorted(xs); return xs[int(q*(len(xs)-1))]
def hdr(t): print('\n'+'='*78+f'\n{t}\n'+'='*78)

# ----------------------------------------------------------------------------- load
def load():
    wb=openpyxl.load_workbook(XLSX,data_only=True); ws=wb["Characterized PPOs"]
    R=list(ws.iter_rows(values_only=True)); H=list(R[0])
    ai,acti=H.index('Accession'),H.index('Activity')
    act={r[ai]:r[acti] for r in R[1:] if r[0]}
    tax={r['accession']:(r['kingdom'],r['phylum'],r['genus']) for r in csv.DictReader(open(TAXf))}
    V={}; CAD={}; SSM={}
    for r in csv.DictReader(open(PVEC)):
        if not r.get('vector'): continue
        parts=[res(x) for x in r['vector'].split('-')]
        while len(parts)<10: parts.append('?')
        a=r['accession']
        V[a]=parts[:10]+[thio(r['thioether'])]
        CAD[a]={p:r.get(f'{p}_cadist') for p in POS[:10]}
        SSM[a]={p:(r.get(f'{p}_ss_match')=='True') for p in POS[:10]}
    # characterised vectors built identically to the dataset (all 84 are in V)
    charv=[(a,act[a],V[a]) for a in act if a in V]
    return act,tax,V,CAD,SSM,charv

def cad(CAD,a,p):
    try: return float(CAD[a].get(p))
    except: return None
def consensus(V,accs):
    out=[]
    for j in range(11):
        c=Counter(V[a][j] for a in accs if V[a][j] not in SKIP); out.append(c.most_common(1)[0][0] if c else '~')
    return out
def calibrate_envelopes(V,CAD):
    """Per-position placement envelope = p90 of the canonical (majority) residue's Ca-displacement.
    Judges each position against its own modelling baseline instead of one global cutoff, so a floppy
    position (Glu195) is held to a fair standard and nothing needs to be excluded by hand."""
    env={}
    for j,p in enumerate(POS[:10]):
        cnt=Counter(V[a][j] for a in V if V[a][j] not in SKIP)
        canon=cnt.most_common(1)[0][0]
        ds=sorted(d for d in (cad(CAD,a,p) for a in V if V[a][j]==canon) if d is not None)
        env[p]=round(ds[int(0.9*(len(ds)-1))],2) if ds else ENV
    return env

# ----------------------------------------------------------------------------- stage A
def stage_A(V,charv):
    hdr("[A] COVERAGE LADDER (descriptive scale)")
    N=len(V); ACC=list(V)
    sigs=len(set(tuple(V[a]) for a in ACC))
    ucv=len(set(tuple(c[2]) for c in charv))
    nloop=sum(1 for a in ACC if V[a][0]=='~')
    dist=Counter(); actmatch=Counter()
    for a in ACC:
        best=99; ba=None
        for _,ac,cv in charv:
            h=hd(V[a],cv)
            if h<best: best=h; ba=ac
        dist[min(best,3)]+=1
        if best==0: actmatch[ba]+=1
    c0,c1,c2=dist[0],dist[0]+dist[1],dist[0]+dist[1]+dist[2]
    print(f"dataset N={N}; distinct signatures={sigs}; characterised unique vectors={ucv}")
    print(f"Gly46 on a loop (wildcarded): {nloop} ({100*nloop/N:.1f}%); 69/84 characterised also loop")
    print(f"  identical (HD0)       : {c0} ({100*c0/N:.1f}%)  [TYR {actmatch['TYR']}]")
    print(f"  within 1 substitution : {c1} ({100*c1/N:.1f}%)")
    print(f"  within 2 substitutions: {c2} ({100*c2/N:.1f}%)")
    print(f"  divergent (HD>=3)     : {N-c2} ({100*(N-c2)/N:.1f}%)   [descriptive only; not the novelty definition]")

# ----------------------------------------------------------------------------- stage B
def stage_B(charv):
    hdr("[B] ACTIVITY NOT SEPARABLE BY RESIDUE SIGNATURE (characterised set)")
    cv={a:(ac,v) for a,ac,v in charv}; accs=list(cv)
    idx={p:i for i,p in enumerate(POS)}
    def col(p,A): return Counter(cv[a][1][idx[p]] for a in accs if cv[a][0]==A)
    tyrF=col('Val218','TYR')['F']; caoxF=col('Val218','CaOx')['F']; ausF=col('Val218','AUS')['F']
    nT=sum(1 for a in accs if cv[a][0]=='TYR'); nC=sum(1 for a in accs if cv[a][0]=='CaOx'); nA=sum(1 for a in accs if cv[a][0]=='AUS')
    print(f"gate Phe (Val218=F): CaOx {caoxF}/{nC}, TYR {tyrF}/{nT}, AUS {ausF}/{nA}  -> does not mark catechol oxidase")
    # exhaustive ceiling
    best={}
    for k in range(1,4):
        bg=0
        for sub in combinations(POS,k):
            sa=defaultdict(set)
            for a in accs: sa[tuple(cv[a][1][idx[p]] for p in sub)].add(cv[a][0])
            g=sum(1 for a in accs if len(sa[tuple(cv[a][1][idx[p]] for p in sub)])==1)/len(accs)
            bg=max(bg,g)
        best[k]=bg
    print(f"best single-activity separation: 1-pos {best[1]:.0%}, 2-pos {best[2]:.0%}, 3-pos {best[3]:.0%}  (ceiling ~96%, TYR/CaOx overlap)")

# ----------------------------------------------------------------------------- stage C
def stage_C(act,tax,V,charv):
    hdr("[C] O-AMINOPHENOL OXIDASE SET (rarest signature; diagnostic Asn on diverse backgrounds)")
    actvecs=defaultdict(set)
    for a,A,v in charv: actvecs[A].add(tuple(v))
    print("dataset structures exactly sharing each activity's vectors:")
    for A in sorted(actvecs):
        n=sum(1 for a in V if a not in act and tuple(V[a]) in actvecs[A])
        print(f"  {A:11} {n}")
    gN=[a for a in V if V[a][0]=='N' and a not in act]
    onAPO=sum(1 for a in gN if tuple(V[a]) in actvecs['oAPO'])
    print(f"\nGly46=Asn (helix) carriers: {len(gN)}; on a non-oAPO background: {len(gN)-onAPO}; "
          f"distinct vectors {len(set(tuple(V[a]) for a in gN))}")
    print(f"  kingdoms: {dict(Counter(tax.get(a,('?',))[0] for a in gN).most_common())}")
    print("  -> candidates to test for nitroso-forming activity (single demonstrated single-residue switch)")

# ----------------------------------------------------------------------------- stage D
def characterise(V,CAD,SSM,tax,charv,j,resi,envelope=ENV):
    accs=[a for a in V if V[a][j]==resi]; n=len(accs)
    sa=accs[::max(1,n//200)] if n>200 else accs
    cohere=med([hd(V[x],V[y]) for x,y in combinations(sa,2)])
    ds=[cad(CAD,a,POS[j]) for a in accs]; ds=[d for d in ds if d is not None]
    mqc=statistics.median(ds); within=sum(d<=envelope for d in ds)/len(ds)
    nd=[min(hd(V[a],cv) for _,_,cv in charv) for a in accs]
    kacc=defaultdict(list)
    for a in accs: kacc[tax.get(a,('?',))[0]].append(a)
    kcons={k:consensus(V,v) for k,v in kacc.items() if k!='?' and len(v)>=5}
    maxconv=max((hd(x,y) for x,y in combinations(kcons.values(),2)),default=0)
    gc=Counter(tax.get(a,('?','?','?'))[2] for a in accs)
    ssfrac=sum(1 for a in accs if SSM.get(a,{}).get(POS[j]))/n
    return dict(accs=accs,n=n,cohere=cohere,mqc=mqc,within=within,nearest=med(nd),
                nking=len(kcons),maxconv=maxconv,gc=gc,ssfrac=ssfrac,cons=consensus(V,accs))

def stage_D(V,CAD,SSM,tax,charv,env_pos):
    hdr(f"[D] NOVELTY ENUMERATION (novel residue, >= {MINN} carriers; per-position placement envelope)")
    cons_cnt=[Counter(V[a][j] for a in V if V[a][j] not in SKIP) for j in range(11)]
    char_sets=[set(c[2][j] for c in charv)-SKIP for j in range(11)]
    def conservation(j):
        c=cons_cnt[j]; return c.most_common(1)[0][1]/sum(c.values())
    rows=[]
    for j,p in enumerate(POS):
        if p=='thioether': continue                       # binary present/absent, not a residue substitution
        env=env_pos[p]
        consv=conservation(j); tier='high' if consv>=0.8 else ('intermediate' if consv>=0.4 else 'variable')
        for resi,c in cons_cnt[j].items():
            if resi in char_sets[j] or c<MINN: continue
            if (p,resi)==THIOETHER_CYS: continue
            d=characterise(V,CAD,SSM,tax,charv,j,resi,envelope=env)
            placed=(d['mqc']<=env and d['within']>=0.5)
            # placement (vs the position's own envelope) is the only DISQUALIFIER; coherence/conservation are descriptors
            verdict='candidate' if placed else 'displaced'
            rows.append(dict(position=p,residue=resi,tier=tier,cons=round(consv,2),n=d['n'],
                coherence_HD=d['cohere'],median_cadist=round(d['mqc'],2),pct_within_env=round(100*d['within']),
                ss_match_pct=round(100*d['ssfrac']),nearest_char_HD=d['nearest'],
                converging_kingdoms=d['nking'],max_xking_HD=d['maxconv'],
                top_genera='; '.join(f"{g}:{m}" for g,m in d['gc'].most_common(3) if g!='?'),verdict=verdict))
    rank={'high':0,'intermediate':1,'variable':2}
    rows.sort(key=lambda r:(0 if r['verdict']=='candidate' else 1, rank[r['tier']], r['coherence_HD'], -r['n']))
    os.makedirs('supplementary',exist_ok=True)
    cols=['position','residue','tier','cons','n','coherence_HD','median_cadist','pct_within_env',
          'ss_match_pct','nearest_char_HD','converging_kingdoms','max_xking_HD','top_genera','verdict']
    with open('supplementary/novelty_enumeration.tsv','w') as f:
        f.write('\t'.join(cols)+'\n')
        for r in rows: f.write('\t'.join(str(r[c]) for c in cols)+'\n')
    nc=sum(r['verdict']=='candidate' for r in rows); nd=sum(r['verdict']=='displaced' for r in rows)
    print(f"{len(rows)} novel states; placed candidates {nc}; displaced (set aside / Cu-test) {nd}")
    print(f"{'state':12} {'tier':12} {'n':>4} {'cohHD':>5} {'medCa':>5} {'inEnv':>6} {'verdict':>10}")
    for r in rows:
        print(f"{r['position']+'='+r['residue']:12} {r['tier']:12} {r['n']:>4} {r['coherence_HD']:>5.1f} "
              f"{r['median_cadist']:>5.2f} {r['pct_within_env']:>5}% {r['verdict']:>10}")
    return rows

# ----------------------------------------------------------------------------- stage F
def stage_F(V,CAD,charv,env_pos):
    hdr("[F] SINGLE-SUBSTITUTION LEADS (HD=1 from a characterised vector; novel, placed residue; ALL positions)")
    char_sets=[set(c[2][j] for c in charv)-SKIP for j in range(11)]
    cons_cnt=[Counter(V[a][j] for a in V if V[a][j] not in SKIP) for j in range(11)]
    def tier(j):
        c=cons_cnt[j]; f=c.most_common(1)[0][1]/sum(c.values())
        return 'high' if f>=0.8 else ('intermediate' if f>=0.4 else 'variable')
    tally=Counter(); adj=defaultdict(Counter); tiers={}
    for a in V:
        if any(a==c[0] for c in charv): continue
        best=99; ref=None
        for _,ac,cv in charv:
            dp=diffpos(V[a],cv)
            if len(dp)<best: best=len(dp); ref=(dp,ac)
        if best!=1: continue                                   # exactly one position differs
        pidx=ref[0][0]; p=POS[pidx]; r=V[a][pidx]
        if p=='thioether' or r in char_sets[pidx] or r in SKIP: continue   # need a novel residue (not the binary thioether)
        c=cad(CAD,a,p)
        if c is not None and c<env_pos[p]:                     # placed within the position's own envelope
            tally[(p,r)]+=1; adj[(p,r)][ref[1]]+=1; tiers[(p,r)]=tier(pidx)
    tot=sum(tally.values()); hi=sum(n for k,n in tally.items() if tiers[k]=='high')
    print(f"{tot} structures one confident substitution from a characterised active site (all positions; {hi} at highly conserved positions):")
    print(f"{'state':12} {'n':>4} {'tier':>13}  one step from")
    for (p,r),n in tally.most_common():
        print(f"  {p+'='+r:11} {n:>3} {tiers[(p,r)]:>13}  {', '.join(f'{a}:{c}' for a,c in adj[(p,r)].most_common(3))}")
    with open('supplementary/single_substitution.tsv','w') as f:
        f.write("position\tresidue\ttier\tn\tone_step_from\n")
        for (p,r),n in tally.most_common():
            f.write(f"{p}\t{r}\t{tiers[(p,r)]}\t{n}\t{'; '.join(f'{a}:{c}' for a,c in adj[(p,r)].most_common())}\n")
    print("  -> supplementary/single_substitution.tsv")
    return tally

# ----------------------------------------------------------------------------- stage E (copper distance; needs CIFs)
def _parse_atoms(path):
    atoms=[]; cols=[]; collecting=False; in_atom=False
    for raw in open(path):
        line=raw.strip()
        if line=='loop_': collecting=False; in_atom=False; cols=[]; continue
        if line.startswith('_atom_site.'): collecting=True; cols.append(line); continue
        if collecting and cols: collecting=False; in_atom=True
        if in_atom:
            if line.startswith('_') or line=='#' or not line: break
            p=line.split()
            if len(p)!=len(cols): continue
            row=dict(zip(cols,p))
            try:
                atoms.append({'elem':row.get('_atom_site.type_symbol','').upper(),
                    'atom':row.get('_atom_site.label_atom_id','').upper(),
                    'resn':row.get('_atom_site.label_comp_id',''),'seq':row.get('_atom_site.label_seq_id',''),
                    'x':float(row['_atom_site.Cartn_x']),'y':float(row['_atom_site.Cartn_y']),'z':float(row['_atom_site.Cartn_z']),
                    'bfactor':float(row.get('_atom_site.B_iso_or_equiv','0') or 0)})
            except (KeyError,ValueError): continue
    return atoms
def _d(a,b): return math.sqrt((a['x']-b['x'])**2+(a['y']-b['y'])**2+(a['z']-b['z'])**2)
def _cu(atoms):
    cu=[a for a in atoms if a['elem']=='CU']; return cu if len(cu)>=2 else None
def _coordhis(atoms,cu,cut=3.0):
    h=[a for a in atoms if a['resn']=='HIS' and a['atom']=='NE2']
    c=[a for a in h if min(_d(a,cu[0]),_d(a,cu[1]))<=cut]; c.sort(key=lambda a:int(a['seq'])); return c
def _kabsch(P,Q):
    import numpy as np
    cP,cQ=P.mean(0),Q.mean(0); H=(P-cP).T@(Q-cQ); U,S,Vt=np.linalg.svd(H)
    d=np.linalg.det(Vt.T@U.T); R=Vt.T@np.diag([1,1,d])@U.T; return R,cQ-R@cP
# stage-H parallel workers: one CIF -> nearest aromatic ring to the reference probe (six-His frame)
def _occ_init(rne2_list):
    global _OCC_RNE2
    import numpy as np
    _OCC_RNE2=np.array(rne2_list)
def _occ_worker(task):
    import numpy as np
    path,target=task
    if not os.path.exists(path): return None
    at=_parse_atoms(path); cu=_cu(at)
    if cu is None: return None
    ch=_coordhis(at,cu)
    if len(ch)<6: return None
    R,t=_kabsch(np.array([[x['x'],x['y'],x['z']] for x in ch[:6]]),_OCC_RNE2)
    tgt=np.array(target); best=None
    for aa,seq,c in _rings(at):
        v=R@np.array([c['x'],c['y'],c['z']])+t; dd=float(np.linalg.norm(v-tgt))
        if best is None or dd<best[0]: best=(dd,aa)
    return ('t',best[0],best[1]) if best is not None else ('t',9e9,'?')

def stage_E(rows, V, tax, cifdir, seqid):
    import numpy as np, glob
    hdr("[E] COPPER-DISTANCE TEST of displaced groups (uniform Cu ruler) + thioether-Cys")
    ref_cif=glob.glob(f"{cifdir}/B2ZB02_taxID_*_model.cif")
    if not ref_cif:
        print("  PmTYR reference CIF not found in --cifs; skipping stage E"); return {}
    ratoms=_parse_atoms(ref_cif[0]); rcu=_cu(ratoms); rne2=np.array([[a['x'],a['y'],a['z']] for a in _coordhis(ratoms,rcu)])
    rca={int(a['seq']):np.array([a['x'],a['y'],a['z']]) for a in ratoms if a['atom']=='CA' and a['seq'].isdigit()}
    def res_to_cu(p,resi,n_sample=25):
        accs=[a for a in V if V[a][POS.index(p)]==resi]; accs=accs[:n_sample]
        dists=[]
        for a in accs:
            mid=seqid.get(a)
            if not mid: continue
            path=f"{cifdir}/{mid}.cif"
            if not os.path.exists(path): continue
            atoms=_parse_atoms(path); cu=_cu(atoms)
            if cu is None: continue
            ch=_coordhis(atoms,cu)
            if len(ch)<6: continue
            qne2=np.array([[x['x'],x['y'],x['z']] for x in ch[:6]]); Rk,t=_kabsch(qne2,rne2)
            qca={int(x['seq']):Rk@np.array([x['x'],x['y'],x['z']])+t for x in atoms if x['atom']=='CA' and x['seq'].isdigit()}
            if PMTYR_NUM[p] not in rca: continue
            seq46=min(qca,key=lambda s:np.linalg.norm(qca[s]-rca[PMTYR_NUM[p]]))
            ratoms_=[x for x in atoms if x['seq']==str(seq46)]
            fa=[x for x in ratoms_ if x['atom'] in FUNC.get(resi,[])]
            if fa: dists.append(min(min(_d(x,cu[0]),_d(x,cu[1])) for x in fa))
        return dists
    rescued={}
    disp=[r for r in rows if r['verdict']=='displaced' and r['n']>=MINN]
    print(f"testing {len(disp)} displaced groups (sampled CIFs); 'pocket' = side chain within {CU_POCKET} A of a copper")
    for r in disp:
        ds=res_to_cu(r['position'],r['residue'])
        if not ds: continue
        m=statistics.median(ds); frac=sum(d<CU_POCKET for d in ds)/len(ds)
        tag='POCKET (structural value-add)' if (m<CU_POCKET and frac>=0.6) else 'not in pocket'
        print(f"  {r['position']+'='+r['residue']:11} n_tested={len(ds):>2} median Cu-dist={m:.1f} A  in-pocket {frac:.0%}  -> {tag}")
        if m<CU_POCKET and frac>=0.6: rescued[(r['position'],r['residue'])]=round(m,1)
    return rescued

# ----------------------------------------------------------------------------- stage H (aromatic ring-slot compensation; needs CIFs)
def _rings(atoms):
    """centroid of every aromatic side-chain ring: list of (one_letter, seq, centroid)."""
    g=defaultdict(list)
    for a in atoms:
        if a['resn'] in RING_ATOMS and a['atom'] in RING_ATOMS[a['resn']]:
            g[(a['seq'],a['resn'])].append(a)
    out=[]
    for (seq,resn),ats in g.items():
        if len(ats)<4: continue
        c={'x':sum(x['x'] for x in ats)/len(ats),'y':sum(x['y'] for x in ats)/len(ats),'z':sum(x['z'] for x in ats)/len(ats)}
        out.append((AA3.get(resn,'?'),seq,c))
    return out
def stage_H(V, cifdir, seqid):
    import numpy as np, glob
    hdr("[H] AROMATIC-RING COMPENSATION (superpose on the six His; probe from the reference ring)\n"
        "    when a conserved aromatic is replaced by a non-aromatic residue, does another\n"
        "    aromatic ring land where the reference ring sat (3D distance in the common frame)?")
    ref_cif=glob.glob(f"{cifdir}/B2ZB02_taxID_*_model.cif")
    if not ref_cif:
        print("  PmTYR reference CIF not found in --cifs; skipping stage H"); return {}
    ra=_parse_atoms(ref_cif[0]); rcu=_cu(ra)
    rne2=np.array([[a['x'],a['y'],a['z']] for a in _coordhis(ra,rcu)[:6]])
    refc={}                                          # reference aromatic ring centroid per position (ref frame)
    for aa,seq,c in _rings(ra):
        for p in AROM_POS:
            if seq==str(PMTYR_NUM[p]): refc[p]=np.array([c['x'],c['y'],c['z']])
    from multiprocessing import Pool
    nproc=int(os.environ.get("SLURM_CPUS_PER_TASK","4"))
    pool=Pool(nproc, initializer=_occ_init, initargs=(rne2.tolist(),))
    def occupancy(accs,p):                                  # ALL carriers, parallel CIF parse
        target=refc[p].tolist()
        tasks=[(f"{cifdir}/{seqid[a]}.cif",target) for a in accs if seqid.get(a)]
        tested=0; hits=0; ident=Counter()
        for r in pool.imap_unordered(_occ_worker,tasks,chunksize=16):
            if r is None: continue
            tested+=1
            if r[1]<=RING_DIST_TOL: hits+=1; ident[r[2]]+=1
        return tested,hits,ident
    out=[]; comp={}
    for p in AROM_POS:
        if p not in refc:
            print(f"  {p}: reference ring not found in PmTYR; skipping"); continue
        pidx=POS.index(p)
        ct,cbase,_=occupancy([a for a in V if V[a][pidx] in AROM],p)    # baseline: aromatic retained
        print(f"\n  {p}  (probe = reference ring centroid; aromatic ring within {RING_DIST_TOL} A after six-His superposition)")
        print(f"    canonical aromatic retained: ring at reference site {100*cbase/max(ct,1):>3.0f}% (n={ct})  [sanity baseline]")
        loss=Counter(V[a][pidx] for a in V if V[a][pidx] not in AROM and V[a][pidx] not in SKIP)
        for r,n in loss.most_common():
            if n<10: continue                                          # only substitutions with >=10 carriers
            t,h,ident=occupancy([a for a in V if V[a][pidx]==r],p)
            if t==0: continue
            frac=h/t; idents=', '.join(f'{k}:{v}' for k,v in ident.most_common(3))
            tag='COMPENSATED' if frac>=0.7 else ('partial' if frac>=0.4 else 'not compensated')
            print(f"    {p+'='+r:9} n={n:>4} tested={t:>2}  aromatic at reference site {100*frac:>3.0f}%  by [{idents}]  -> {tag}")
            out.append((p,r,n,t,round(100*frac),idents,tag))
            if tag=='COMPENSATED': comp[(p,r)]=(round(100*frac),idents)
    pool.close(); pool.join()
    with open('supplementary/aromatic_compensation.tsv','w') as f:
        f.write("position\tnonaromatic_residue\tn_carriers\tn_tested\taromatic_at_ref_site_pct\tcompensating_rings\tverdict\n")
        for p,r,n,t,pc,idents,tag in out:
            f.write(f"{p}\t{r}\t{n}\t{t}\t{pc}\t{idents}\t{tag}\n")
    print("\n  -> supplementary/aromatic_compensation.tsv")
    return comp     # {(position,non-aromatic residue): (aromatic_at_ref_site_pct, compensating rings)}

# ----------------------------------------------------------------------------- stages J,K (geometry + pLDDT; needs CIFs)
def _ref_variable(cifdir):
    """PmTYR reference: six-His NE2 frame + Ca coords at the 10 variable positions, or None."""
    import glob
    ref_cif=glob.glob(f"{cifdir}/B2ZB02_taxID_*_model.cif")
    if not ref_cif: return None
    ra=_parse_atoms(ref_cif[0]); rcu=_cu(ra)
    if rcu is None: return None
    rne2=[[a['x'],a['y'],a['z']] for a in _coordhis(ra,rcu)[:6]]
    rca={int(a['seq']):[a['x'],a['y'],a['z']] for a in ra if a['atom']=='CA' and a['seq'].isdigit()}
    refpos=[rca[PMTYR_NUM[p]] for p in POS[:10] if PMTYR_NUM[p] in rca]
    return (rne2,refpos) if len(refpos)==10 else None
def _map_init(rne2,refpos):
    global _MAP_RNE2,_MAP_REFPOS
    import numpy as np
    _MAP_RNE2=np.array(rne2); _MAP_REFPOS=[np.array(c) for c in refpos]
def _superpose(at):
    """Six-His Kabsch onto the reference -> transformed query Ca {seq: xyz}, or None."""
    import numpy as np
    cu=_cu(at)
    if cu is None: return None
    ch=_coordhis(at,cu)
    if len(ch)<6: return None
    R,t=_kabsch(np.array([[x['x'],x['y'],x['z']] for x in ch[:6]]),_MAP_RNE2)
    return {int(x['seq']):R@np.array([x['x'],x['y'],x['z']])+t
            for x in at if x['atom']=='CA' and x['seq'].isdigit()}
def _rmsd_worker(task):
    import numpy as np
    path,acc=task
    if not os.path.exists(path): return None
    qca=_superpose(_parse_atoms(path))
    if not qca: return None
    return (acc,[min(qca.values(),key=lambda v:np.linalg.norm(v-refc)).tolist() for refc in _MAP_REFPOS])
def _plddt_worker(task):
    import numpy as np
    path,acc=task
    if not os.path.exists(path): return None
    at=_parse_atoms(path); qca=_superpose(at)
    if not qca: return None
    pl={int(x['seq']):x['bfactor'] for x in at if x['atom']=='CA' and x['seq'].isdigit()}
    return [(j,pl[bs]) for j,refc in enumerate(_MAP_REFPOS)
            for bs in [min(qca,key=lambda s:np.linalg.norm(qca[s]-refc))] if bs in pl]

def stage_J(V, cifdir, seqid, min_group=5):
    import numpy as np
    from multiprocessing import Pool
    hdr(f"[J] WITHIN-GROUP GEOMETRIC CONSISTENCY (Ca RMSD of the 10 variable positions to the group\n"
        f"    centroid; identical-vector groups with >= {min_group} members)")
    ref=_ref_variable(cifdir)
    if ref is None: print("  reference CIF/positions unavailable; skipping"); return
    rne2,refpos=ref
    vec={r['accession']:r['vector'] for r in csv.DictReader(open(PVEC)) if r.get('vector') and not r.get('error')}
    groups=defaultdict(list)
    for a,v in vec.items(): groups[v].append(a)
    groups={v:a for v,a in groups.items() if len(a)>=min_group}
    members=[a for accs in groups.values() for a in accs if a in seqid]
    print(f"  vector groups >= {min_group}: {len(groups)}; structures to superpose: {len(members)}")
    tasks=[(f"{cifdir}/{seqid[a]}.cif",a) for a in members]
    nproc=int(os.environ.get("SLURM_CPUS_PER_TASK","4"))
    coords={}
    with Pool(nproc,initializer=_map_init,initargs=(rne2,refpos)) as pool:
        for r in pool.imap_unordered(_rmsd_worker,tasks,chunksize=50):
            if r: coords[r[0]]=np.array(r[1])
    gm=[]
    for v,accs in groups.items():
        mc=[coords[a] for a in accs if a in coords]
        if len(mc)<2: continue
        cen=np.array(mc).mean(0)
        gm.append(float(np.mean([np.sqrt(np.mean(np.sum((c-cen)**2,axis=1))) for c in mc])))
    gm=np.array(gm)
    print(f"  structures mapped: {len(coords)}; groups with >= 2 mapped: {len(gm)}")
    print(f"  mean within-group RMSD: {gm.mean():.3f} A (median {np.median(gm):.3f} A)")

def stage_K(V, cifdir, seqid, n_sample=2000):
    import numpy as np, random
    from multiprocessing import Pool
    hdr(f"[K] CONFIDENCE AT THE MAPPED POSITIONS (per-position Ca pLDDT; random sample n={n_sample})")
    ref=_ref_variable(cifdir)
    if ref is None: print("  reference CIF/positions unavailable; skipping"); return
    rne2,refpos=ref
    accs=[a for a in V if a in seqid]; random.seed(0)
    sample=random.sample(accs,min(n_sample,len(accs)))
    tasks=[(f"{cifdir}/{seqid[a]}.cif",a) for a in sample]
    nproc=int(os.environ.get("SLURM_CPUS_PER_TASK","4"))
    pv=defaultdict(list)
    with Pool(nproc,initializer=_map_init,initargs=(rne2,refpos)) as pool:
        for r in pool.imap_unordered(_plddt_worker,tasks,chunksize=50):
            if r:
                for j,v in r: pv[j].append(v)
    print(f"  {'position':9} {'mean':>6} {'median':>7} {'min':>6} {'%<70':>6}")
    means=[]
    for j,p in enumerate(POS[:10]):
        vals=pv.get(j,[])
        if not vals: continue
        arr=np.array(vals); means.append(arr.mean())
        print(f"  {p:9} {arr.mean():6.1f} {np.median(arr):7.1f} {arr.min():6.1f} {100*np.mean(arr<70):5.1f}%")
    if means: print(f"\n  variable-position mean pLDDT range: {min(means):.1f}-{max(means):.1f}")

# ----------------------------------------------------------------------------- HMM (stages G + I; needs --afa)
HMM_COL={'Gly46':15,'Phe65':30,'Trp68':33,'Glu195':163,'Asn205':175,'Arg209':179,
         'Val218':188,'Ala221':191,'Phe227':197,'His230':200}   # match-state column per position
def _hmm_positions(afa):
    """Parse the HMM .afa once -> {accession: {position: aligned residue}} at the HMM_COL match columns."""
    MAXC=max(HMM_COL.values()); want={v:k for k,v in HMM_COL.items()}
    def extract(al):
        mc=0; out={}
        for ch in al:
            if ch=='-' or ch.isupper(): mc+=1
            else: continue
            if mc in want: out[want[mc]]=ch
            if mc>=MAXC: break
        return out
    print("  parsing HMM alignment...", file=sys.stderr)
    hp={}; name=""; buf=[]
    for line in open(afa):
        line=line.rstrip('\n')
        if line.startswith('>'):
            if name: hp[name]=extract("".join(buf))
            name=line[1:].split()[0].split('|')[0]; buf=[]
        else: buf.append(line)
    if name: hp[name]=extract("".join(buf))
    return hp

def _hmm_positions_from_cols(path):
    """Same {accession: {position: residue}} map as _hmm_positions, read from the
    compact match-state table written by hmm/build_alignment.py (one row per
    sequence, one character per PF00264 match state). Equivalent to parsing the
    .afa, minus the 782 MB."""
    import gzip
    opener = gzip.open if str(path).endswith('.gz') else open
    print("  reading HMM match-state columns...", file=sys.stderr)
    want = {v: k for k, v in HMM_COL.items()}
    hp = {}
    with opener(path, 'rt') as fh:
        header = fh.readline()
        if not header.startswith('accession'):
            fh.seek(0)
        for line in fh:
            acc, _, cols = line.rstrip('\n').partition('\t')
            if not cols:
                continue
            hp[acc] = {p: cols[c - 1] for c, p in want.items() if c <= len(cols)}
    return hp


def stage_G(V, hp, flagged):
    hdr("[G] HMM CROSS-REFERENCE of the flagged groups (descriptive only; the structure is primary)\n"
        "    high agreement = sequence-recoverable; low agreement = a structure-specific finding")
    print(f"{'group':14} {'n':>5} {'inHMM':>6} {'agree%':>7}  (HMM finds the same residue?)")
    for p,r in flagged:
        if p not in HMM_COL:
            print(f"  {p+'='+r:12} (Gly46 not on HMM grid -> structural value-add, see stage E)"); continue
        mem=[a for a in V if V[a][POS.index(p)]==r]
        inh=[a for a in mem if a in hp and p in hp[a]]
        ag=sum(1 for a in inh if hp[a][p]==r)
        print(f"  {p+'='+r:12} {len(mem):>5} {len(inh):>6} {100*ag/max(len(inh),1):>6.1f}%")

def stage_I(V, hp):
    hdr("[I] HMM PER-POSITION AGREEMENT (independent profile alignment vs the structural mapping)\n"
        "    Table SX: fraction of structures whose HMM-column residue matches the structural assignment")
    cons_cnt=[Counter(V[a][j] for a in V if V[a][j] not in SKIP) for j in range(10)]
    rows=[]
    print(f"  {'position':9} {'tier':13} {'cons%':>6} {'inHMM':>7} {'agree%':>7}")
    for j,p in enumerate(POS[:10]):
        if p not in HMM_COL: continue
        c=cons_cnt[j]; f=c.most_common(1)[0][1]/sum(c.values())
        tier='high' if f>=0.8 else ('intermediate' if f>=0.4 else 'variable')
        inh=[a for a in V if a in hp and p in hp[a] and V[a][j] not in SKIP]
        ag=sum(1 for a in inh if hp[a][p]==V[a][j]); agree=100*ag/max(len(inh),1)
        rows.append((p,tier,100*f,len(inh),agree))
        print(f"  {p:9} {tier:13} {100*f:5.0f}% {len(inh):>7} {agree:6.1f}%")
    hi=[r for r in rows if r[1]=='high']
    if hi:
        print(f"\n  conserved positions ({', '.join(r[0] for r in hi)}): agreement "
              f"{min(r[4] for r in hi):.1f}-{max(r[4] for r in hi):.1f}%")
    os.makedirs('supplementary',exist_ok=True)
    with open('supplementary/hmm_agreement.tsv','w') as fo:
        fo.write("position\ttier\tconservation_pct\tn_in_hmm\tagreement_pct\n")
        for p,t,f,n,a in rows: fo.write(f"{p}\t{t}\t{f:.1f}\t{n}\t{a:.1f}\n")
    print("  -> supplementary/hmm_agreement.tsv")

# ----------------------------------------------------------------------------- main
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--cifs', help="mounted squashfs dir with *_model.cif (enables stages E,H,J,K)")
    ap.add_argument('--afa',  help="all_hmmalign.afa (enables stages G,I)")
    ap.add_argument('--hmm-cols', dest='hmm_cols',
                    help="hmm/hmm_match_columns.tsv.gz from build_alignment.py "
                         "(deposited stand-in for --afa; enables stages G,I)")
    ap.add_argument('--plddt-sample', type=int, default=2000,
                    help="random sample size for the stage-K pLDDT estimate")
    args=ap.parse_args()
    act,tax,V,CAD,SSM,charv=load()
    print(f"loaded {len(V)} structures, {len(charv)} characterised")
    seqid={}
    if args.cifs:
        for fn in os.listdir(args.cifs):
            if fn.endswith('_model.cif'):
                stem=fn[:-4]
                seqid[stem.split('_taxID_')[0]]=stem

    env_pos=calibrate_envelopes(V,CAD)
    print("per-position placement envelopes (p90 of canonical residue):",
          {p:env_pos[p] for p in ['Phe227','Phe65','Trp68','His230','Glu195','Asn205','Ala221','Gly46','Val218','Arg209']})
    stage_A(V,charv)
    stage_B(charv)
    stage_C(act,tax,V,charv)
    rows=stage_D(V,CAD,SSM,tax,charv,env_pos)
    stage_F(V,CAD,charv,env_pos)
    rescued=stage_E(rows,V,tax,args.cifs,seqid) if args.cifs else {}
    comp=stage_H(V,args.cifs,seqid) if args.cifs else {}
    if args.cifs:
        stage_J(V,args.cifs,seqid)
        stage_K(V,args.cifs,seqid,n_sample=args.plddt_sample)

    # ---- assemble the flagged shortlist ----
    # Placement DISQUALIFIES (displaced -> out, unless Cu-rescued). Among placed states,
    # conservation + coherence + size PRIORITISE: the shortlist is the placed, conserved,
    # coherent, sized groups, plus the controller-charge theme, plus Cu-rescued displaced.
    # (The prose then features a chemically-curated subset of this shortlist, e.g. the
    #  His230 aromatics as a group; that final selection is a documented biological judgment.)
    HEADLINE_COH=5; HEADLINE_N=40
    hdr("FLAGGED FOR CHARACTERISATION (shortlist; full ranked set in novelty_enumeration.tsv)")
    flagged=[]; demoted=[]
    for r in rows:
        p,resi=r['position'],r['residue']
        if r['verdict']=='candidate':
            if (p,resi) in comp:     # aromatic-loss whose Cu-relative ring slot stays filled (stage H) -> not a true loss
                demoted.append((p,resi,r['n'],comp[(p,resi)])); continue
            if r['tier']=='high' and r['coherence_HD']<=HEADLINE_COH and r['n']>=HEADLINE_N:
                flagged.append((p,resi,r['tier'],r['n'],f"conserved position, placed, coherent (HD{r['coherence_HD']:.0f})"))
            elif p=='Asn205' and resi in {'K','R','E'}:
                flagged.append((p,resi,r['tier'],r['n'],"charge introduced at the Asn205 controller (never charged in characterised)"))
            elif p=='Glu195' and resi in {'R','K','H'}:
                flagged.append((p,resi,r['tier'],r['n'],"positive charge at the Glu195 catalytic residue (charge reversal; never in characterised)"))
        elif (p,resi) in rescued:
            flagged.append((p,resi,r['tier'],r['n'],
                            f"displaced but Cu-confirmed pocket ({rescued[(p,resi)]} A) [structural value-add]"))
    flagged.sort(key=lambda x:({'high':0,'intermediate':1,'variable':2}[x[2]],-x[3]))
    with open('supplementary/flagged_groups.tsv','w') as f:
        f.write("position\tresidue\ttier\tn\tbasis\n")
        for p,r,t,n,b in flagged: f.write(f"{p}\t{r}\t{t}\t{n}\t{b}\n")
    print(f"{'group':12} {'tier':12} {'n':>5}  basis")
    for p,r,t,n,b in flagged: print(f"  {p+'='+r:11} {t:12} {n:>5}  {b}")
    if demoted:
        print("\nDEMOTED by stage H (aromatic-loss but ring slot stays filled -> not a true aromatic loss):")
        for p,r,n,(pc,idents) in demoted:
            print(f"  {p+'='+r:11} n={n:>4}  slot filled {pc}% by [{idents}]")
    print("\n-> supplementary/flagged_groups.tsv ; supplementary/novelty_enumeration.tsv")

    hp=None
    if args.afa:
        hp=_hmm_positions(args.afa)
    elif args.hmm_cols:
        hp=_hmm_positions_from_cols(args.hmm_cols)
    if hp:
        stage_I(V,hp)
        stage_G(V,hp,[(p,r) for p,r,_,_,_ in flagged])

if __name__=='__main__':
    main()
