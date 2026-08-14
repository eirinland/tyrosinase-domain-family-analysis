"""
Stage-3 non-canonical extraction + active-site mapping (v3).

Replaces the old CEAligner whole-structure superposition with an ICP-refined,
helix-anchored superposition onto B2ZB02 (PmTYR AF3 model).  Validated on 800
canonical structures (ground-truth 6 His): 97.9% all-6 recovery, fold-independent
(holds in the qTM<0.5 divergent bin), because the di-Cu 6-His active site is
structurally conserved across all PPO references — so a single-reference local
superposition suffices and per-reference His mapping is unnecessary.

Per structure:
  - seed helix-Kabsch from the query-vs-B2ZB02 Foldseek cigar (core helices only)
  - ICP refinement on B2ZB02 core-helix CA (correspondence-free, geometry-driven)
  - read the residue at each of the 6 canonical His positions (nearest query CA)
  - Cu assignment to CuA/CuB (query Cu transformed into B2ZB02 frame)
  - 5 A coordination shell per Cu, global pLDDT, per-Cu pLDDT, taxonomy

Output noncanonical_analysis.tsv mirrors the old schema so classify_sites.py
runs unchanged.
"""
import argparse, csv, os, re, math
from collections import defaultdict, Counter
import numpy as np

AA3 = {'ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS',
       'MET','PHE','PRO','SER','THR','TRP','TYR','VAL'}
IMIDAZOLE = {'ND1','NE2','CD2','CE1'}

# B2ZB02 core-helix residue ranges (a1-a4) and the 6 canonical His anchors
HELIX_RANGES = [(34,46),(65,83),(203,211),(226,244)]
ANCHORS = [42,60,69,204,208,231]

COORD_ATOMS = {
    'GLU': ['OE1', 'OE2'], 'ASP': ['OD1', 'OD2'],
    'CYS': ['SG'], 'MET': ['SD'], 'TYR': ['OH'],
}


def parse_atoms(path):
    atoms=[]; cols=[]; inA=False; coll=False
    with open(path) as f:
        for raw in f:
            l=raw.strip()
            if l=='loop_': coll=False; inA=False; cols=[]; continue
            if l.startswith('_atom_site.'): coll=True; cols.append(l); continue
            if coll and cols: coll=False; inA=True
            if inA:
                if l.startswith('_') or l=='#' or not l: break
                p=l.split()
                if len(p)!=len(cols): continue
                r=dict(zip(cols,p))
                try:
                    atoms.append({'elem':r['_atom_site.type_symbol'].upper(),
                        'atom':r['_atom_site.label_atom_id'].upper(),
                        'resn':r['_atom_site.label_comp_id'].upper(),
                        'seq':r['_atom_site.label_seq_id'],
                        'bfac':float(r.get('_atom_site.B_iso_or_equiv','0') or 0),
                        'x':float(r['_atom_site.Cartn_x']),
                        'y':float(r['_atom_site.Cartn_y']),
                        'z':float(r['_atom_site.Cartn_z'])})
                except (KeyError,ValueError): continue
    return atoms


def d3(a,b): return math.sqrt((a['x']-b['x'])**2+(a['y']-b['y'])**2+(a['z']-b['z'])**2)

def vec(a): return np.array([a['x'],a['y'],a['z']])

def get_cu(atoms):
    cu=[a for a in atoms if a['elem']=='CU']
    cu.sort(key=lambda a:a.get('seq',''))
    return cu

def coord_his_ne2(atoms,cus,cutoff=3.5):
    if len(cus)<2: return []
    ne2=[a for a in atoms if a['resn']=='HIS' and a['atom']=='NE2']
    c=[a for a in ne2 if min(d3(a,cus[0]),d3(a,cus[1]))<=cutoff]
    c.sort(key=lambda a:int(a['seq']))
    return c

def ca_map(atoms):
    m={}
    for a in atoms:
        if a['atom']=='CA' and a['seq'].isdigit():
            m[int(a['seq'])]=vec(a)
    return m

def resn3_map(atoms):
    m={}
    for a in atoms:
        if a['seq'].isdigit() and a['resn'] in AA3:
            s=int(a['seq'])
            if s not in m: m[s]=a['resn']
    return m

def atom_xyz_map(atoms):
    m=defaultdict(dict)
    for a in atoms:
        if a['seq'].isdigit():
            m[int(a['seq'])][a['atom']]=np.array([a['x'],a['y'],a['z']])
    return m

def ordered_seqs(atoms):
    return sorted({int(a['seq']) for a in atoms if a['atom']=='CA' and a['seq'].isdigit()})

def kabsch(P,Q):
    cP,cQ=P.mean(0),Q.mean(0)
    H=(P-cP).T@(Q-cQ)
    U,S,Vt=np.linalg.svd(H)
    d=np.linalg.det(Vt.T@U.T)
    R=Vt.T@np.diag([1,1,d])@U.T
    return R, cQ-R@cP

def parse_cigar(c): return [(op,int(n)) for n,op in re.findall(r'(\d+)([MID])',c)]

def cigar_pairs(qstart,tstart,cigar,qseqs,tseqs):
    qi=qstart-1; ti=tstart-1; pairs=[]
    for op,L in parse_cigar(cigar):
        if op=='M':
            for _ in range(L):
                if 0<=qi<len(qseqs) and 0<=ti<len(tseqs):
                    pairs.append((qseqs[qi],tseqs[ti]))
                qi+=1; ti+=1
        elif op=='I': qi+=L
        elif op=='D': ti+=L
    return pairs

def in_helix(t): return any(s<=t<=e for s,e in HELIX_RANGES)

def ref_core_seqs(ref_ca):
    return [s for s in ref_ca if in_helix(s)]

def helix_seed(qca,qseqs,qstart,tstart,cigar,ref_ca,ref_seqs):
    pairs=cigar_pairs(qstart,tstart,cigar,qseqs,ref_seqs)
    P=[]; Q=[]
    for q,t in pairs:
        if in_helix(t) and q in qca and t in ref_ca:
            P.append(qca[q]); Q.append(ref_ca[t])
    if len(P)<8: return None
    R,t=kabsch(np.array(P),np.array(Q))
    return R,t,len(P)

def icp_refine(qca,R,t,ref_ca,rcs,cutoff=4.0,iters=6):
    qseqs=list(qca.keys()); Qbase=np.array([qca[s] for s in qseqs])
    ref_core=np.array([ref_ca[s] for s in rcs])
    last=999.0; n=0
    for _ in range(iters):
        Qt=(R@Qbase.T).T+t
        P=[]; Q=[]
        for rc in ref_core:
            dd=np.linalg.norm(Qt-rc,axis=1); j=int(dd.argmin())
            if dd[j]<=cutoff: P.append(Qbase[j]); Q.append(rc)
        if len(P)<8: break
        R,t=kabsch(np.array(P),np.array(Q))
        Pt=(R@np.array(P).T).T+t
        last=float(np.sqrt(np.mean(np.sum((Pt-np.array(Q))**2,1)))); n=len(P)
    return R,t,last,n

def nearest_res(anchor_seq,ref_ca,qca_t,qres3,cutoff):
    rc=ref_ca[anchor_seq]; bd=999.0; bs=None
    for s,c in qca_t.items():
        dd=float(np.linalg.norm(c-rc))
        if dd<bd: bd=dd; bs=s
    if bs is not None and bd<=cutoff:
        return qres3.get(bs,'UNK'), bd, bs
    return '---', bd, None


def get_coordination(atoms,cu,cutoff):
    cores=defaultdict(set)
    cuv=vec(cu)
    for a in atoms:
        if a['resn'] not in AA3: continue
        if not a['seq'].isdigit(): continue
        if np.linalg.norm(vec(a)-cuv)<=cutoff:
            cores[a['resn']].add(int(a['seq']))
    return {k:len(v) for k,v in sorted(cores.items())}

def fmt_coord(d):
    return ','.join(f"{k}:{v}" for k,v in sorted(d.items())) if d else ''


def load_taxdump(d):
    nodes={};
    with open(os.path.join(d,'nodes.dmp')) as f:
        for line in f:
            p=[x.strip() for x in line.split('|')]
            nodes[int(p[0])]=(int(p[1]),p[2])
    names={}
    with open(os.path.join(d,'names.dmp')) as f:
        for line in f:
            p=[x.strip() for x in line.split('|')]
            if p[3]=='scientific name': names[int(p[0])]=p[1]
    merged={}
    mp=os.path.join(d,'merged.dmp')
    if os.path.exists(mp):
        with open(mp) as f:
            for line in f:
                p=[x.strip() for x in line.split('|')]
                merged[int(p[0])]=int(p[1])
    return nodes,names,merged

def lineage(taxid,nodes,names,merged):
    if taxid in merged: taxid=merged[taxid]
    out={}; cur=taxid; seen=set()
    while cur and cur not in seen and cur!=1:
        seen.add(cur)
        if cur not in nodes: break
        parent,rank=nodes[cur]
        if rank!='no rank' and cur in names: out[rank]=names[cur]
        cur=parent
    return out

def extract_taxid(fn):
    m=re.search(r'taxID_(\d+)',fn)
    return int(m.group(1)) if m else None

def assign_cu(da,db,cutoff):
    if da<=cutoff and db<=cutoff: return 'CuA' if da<db else 'CuB'
    if da<=cutoff: return 'CuA'
    if db<=cutoff: return 'CuB'
    return 'neither'


def build_ref_anchor_sites(ref):
    """Split the 6 B2ZB02 anchors into ordered CuA/CuB His positions.
    CuA = site whose anchors are N-terminal (lower seq); ordered by residue number (sequential)."""
    cus=get_cu(ref)
    ref_ca=ca_map(ref)
    # nearest Cu (by imidazole) for each anchor His
    imid=defaultdict(list)
    for a in ref:
        if a['resn']=='HIS' and a['atom'] in IMIDAZOLE and a['seq'].isdigit():
            imid[int(a['seq'])].append(a)
    site={}
    for p in ANCHORS:
        da=min(d3(x,cus[0]) for x in imid[p]); db=min(d3(x,cus[1]) for x in imid[p])
        site[p]=(0 if da<db else 1, min(da,db))
    groups={0:[],1:[]}
    for p in ANCHORS: groups[site[p][0]].append(p)
    # label the group with smaller min-seq as CuA
    g0min=min(groups[0]); g1min=min(groups[1])
    cuA_g,cuB_g=(0,1) if g0min<g1min else (1,0)
    cuA=sorted(groups[cuA_g])
    cuB=sorted(groups[cuB_g])
    return cus,ref_ca,cuA,cuB


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--cif-dir',required=True)
    ap.add_argument('--pool-csv',required=True)   # accession,failed_step,n_his,cu_dist
    ap.add_argument('--pmtyr',required=True)       # B2ZB02 AF3 model CIF (2 Cu)
    ap.add_argument('--fs-tsv',required=True)      # query vs B2ZB02 foldseek w/ cigar
    ap.add_argument('--taxdump-dir',required=True)
    ap.add_argument('--output',required=True)
    ap.add_argument('--coord-cutoff',type=float,default=5.0)
    ap.add_argument('--assign-cutoff',type=float,default=5.0)
    ap.add_argument('--ca-cutoff',type=float,default=3.0)
    a=ap.parse_args()

    print('Indexing CIFs...',flush=True)
    acc2cif={}
    for fn in os.listdir(a.cif_dir):
        if fn.endswith('.cif'):
            acc2cif[fn.split('_taxID_')[0]]=fn
    print(f'  {len(acc2cif)} CIFs',flush=True)

    pool=list(csv.DictReader(open(a.pool_csv)))
    print(f'Non-canonical pool: {len(pool)}',flush=True)

    print('Loading taxonomy...',flush=True)
    nodes,tnames,merged=load_taxdump(a.taxdump_dir)

    print('Parsing B2ZB02 reference...',flush=True)
    ref=parse_atoms(a.pmtyr)
    ref_cus,ref_ca,cuA_pos,cuB_pos=build_ref_anchor_sites(ref)
    ref_seqs=ordered_seqs(ref); rcs=ref_core_seqs(ref_ca)
    # B2ZB02 CuA/CuB centroids: CuA = Cu nearer cuA_pos imidazoles
    print(f'  CuA His: {cuA_pos}   CuB His: {cuB_pos}',flush=True)
    # decide which physical Cu is CuA
    imid=defaultdict(list)
    for at in ref:
        if at['resn']=='HIS' and at['atom'] in IMIDAZOLE and at['seq'].isdigit():
            imid[int(at['seq'])].append(at)
    cuA_to0=sum(min(d3(x,ref_cus[0]) for x in imid[p]) for p in cuA_pos)
    cuA_to1=sum(min(d3(x,ref_cus[1]) for x in imid[p]) for p in cuA_pos)
    ref_CuA=vec(ref_cus[0] if cuA_to0<cuA_to1 else ref_cus[1])
    ref_CuB=vec(ref_cus[1] if cuA_to0<cuA_to1 else ref_cus[0])

    # foldseek best hit per accession
    fs={}
    with open(a.fs_tsv) as f:
        for line in f:
            p=line.rstrip('\n').split('\t')
            if len(p)<8: continue
            acc=p[0].split('_taxID_')[0]
            qstart=int(p[2]); tstart=int(p[5]); cigar=p[-1]
            qtm=float(p[7]) if p[7] else 0.0
            if acc not in fs or qtm>fs[acc][3]:
                fs[acc]=(qstart,tstart,cigar,qtm)

    POS_LABELS=[('CuA_His1',cuA_pos[0]),('CuA_His2',cuA_pos[1]),('CuA_His3',cuA_pos[2]),
                ('CuB_His1',cuB_pos[0]),('CuB_His2',cuB_pos[1]),('CuB_His3',cuB_pos[2])]

    results=[]; no_fs=0; no_seed=0; no_cu=0
    for i,row in enumerate(pool,1):
        acc=row['accession']
        fn=acc2cif.get(acc)
        if not fn: continue
        if acc not in fs: no_fs+=1; continue
        qa=parse_atoms(os.path.join(a.cif_dir,fn))
        cus=get_cu(qa)
        if len(cus)<2: no_cu+=1; continue

        gplddt=None
        for at in qa:  # global pLDDT fallback if not in atoms; parse header separately below
            pass
        # global pLDDT from header
        gplddt=read_global_plddt(os.path.join(a.cif_dir,fn))

        cu1,cu2=cus[0],cus[1]
        cu1_coord=get_coordination(qa,cu1,a.coord_cutoff)
        cu2_coord=get_coordination(qa,cu2,a.coord_cutoff)

        qca=ca_map(qa); qres3=resn3_map(qa); qseqs=ordered_seqs(qa)
        qstart,tstart,cigar,qtm=fs[acc]
        seed=helix_seed(qca,qseqs,qstart,tstart,cigar,ref_ca,ref_seqs)
        pos_res={}; icp_rmsd=None
        if seed is None:
            no_seed+=1
            for lab,_ in POS_LABELS: pos_res[lab]=('---',999.0,None)
            R=t=None
        else:
            R,t,_=seed
            R,t,icp_rmsd,_=icp_refine(qca,R,t,ref_ca,rcs)
            qca_t={s:R@c+t for s,c in qca.items()}
            for lab,p in POS_LABELS:
                pos_res[lab]=nearest_res(p,ref_ca,qca_t,qres3,a.ca_cutoff)

        # Cu assignment in B2ZB02 frame
        if R is not None:
            cu1v=R@vec(cu1)+t; cu2v=R@vec(cu2)+t
            d1a=float(np.linalg.norm(cu1v-ref_CuA)); d1b=float(np.linalg.norm(cu1v-ref_CuB))
            d2a=float(np.linalg.norm(cu2v-ref_CuA)); d2b=float(np.linalg.norm(cu2v-ref_CuB))
        else:
            d1a=d1b=d2a=d2b=999.0
        cu1_as=assign_cu(d1a,d1b,a.assign_cutoff)
        cu2_as=assign_cu(d2a,d2b,a.assign_cutoff)
        if cu1_as==cu2_as and cu1_as!='neither':
            s1=d1a if cu1_as=='CuA' else d1b
            s2=d2a if cu2_as=='CuA' else d2b
            if s1<=s2:
                other='CuB' if cu1_as=='CuA' else 'CuA'
                od=d2a if other=='CuA' else d2b
                cu2_as=other if od<=a.assign_cutoff else 'neither'
            else:
                other='CuB' if cu2_as=='CuA' else 'CuA'
                od=d1a if other=='CuA' else d1b
                cu1_as=other if od<=a.assign_cutoff else 'neither'

        taxid=extract_taxid(fn); lin=lineage(taxid,nodes,tnames,merged) if taxid else {}
        res={'accession':acc,'failed_step':row.get('failed_step',''),
             'rmsd':f"{icp_rmsd:.3f}" if icp_rmsd is not None else '',
             'cu_cu_distance':row.get('cu_dist',''),
             'cu1_assignment':cu1_as,'cu1_dist_CuA':f"{d1a:.2f}",'cu1_dist_CuB':f"{d1b:.2f}",
             'cu1_plddt':f"{cu1['bfac']:.1f}",'cu1_coordination':fmt_coord(cu1_coord),
             'cu2_assignment':cu2_as,'cu2_dist_CuA':f"{d2a:.2f}",'cu2_dist_CuB':f"{d2b:.2f}",
             'cu2_plddt':f"{cu2['bfac']:.1f}",'cu2_coordination':fmt_coord(cu2_coord),
             'global_plddt':f"{gplddt:.1f}" if gplddt else '',
             'hmm_coverage':row.get('hmm_coverage',''),
             'n_coord_his':row.get('n_his',''),
             'species':lin.get('species',''),'genus':lin.get('genus',''),
             'family':lin.get('family',''),'order':lin.get('order',''),
             'phylum':lin.get('phylum',''),'superkingdom':lin.get('superkingdom','')}
        qaxyz=atom_xyz_map(qa)
        for lab,_ in POS_LABELS:
            r3,dist,rnum=pos_res[lab]
            res[lab]=r3
            res[f'{lab}_ca_dist']=f"{dist:.2f}" if dist<900 else ''
            res[f'{lab}_resnum']=str(rnum) if rnum is not None else ''
            cdist=''
            if R is not None and rnum is not None and r3 in COORD_ATOMS:
                cu_ref=ref_CuA if lab.startswith('CuA') else ref_CuB
                best_cd=None
                for aname in COORD_ATOMS[r3]:
                    if aname in qaxyz[rnum]:
                        coord_t=R@qaxyz[rnum][aname]+t
                        cd=float(np.linalg.norm(coord_t-cu_ref))
                        if best_cd is None or cd<best_cd: best_cd=cd
                if best_cd is not None: cdist=f"{best_cd:.2f}"
            res[f'{lab}_coord_dist']=cdist
        results.append(res)
        if i%200==0: print(f'  [{i}/{len(pool)}]',flush=True)

    print(f'\nProcessed {len(results)}  (no_fs={no_fs} no_cu={no_cu} no_seed={no_seed})',flush=True)
    if not results: return
    fields=list(results[0].keys())
    with open(a.output,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t'); w.writeheader()
        w.writerows(results)
    asg=Counter(f"{r['cu1_assignment']}+{r['cu2_assignment']}" for r in results)
    print('Cu assignment pairs:')
    for k,n in asg.most_common(): print(f'  {k}: {n}')
    print(f'Written {len(results)} rows to {a.output}')


def read_global_plddt(cif_path):
    with open(cif_path) as f:
        for line in f:
            if line.startswith('_ma_qa_metric_global.metric_value'):
                parts=line.strip().split()
                if len(parts)>=2:
                    try: return float(parts[1])
                    except ValueError: pass
                nxt=next(f,'').strip()
                try: return float(nxt)
                except (ValueError,StopIteration): pass
    return None


if __name__=='__main__': main()
