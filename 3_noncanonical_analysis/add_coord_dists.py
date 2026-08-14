#!/usr/bin/env python3
"""Add coord_dist columns to noncanonical_analysis.tsv using helix-anchored Kabsch+ICP.

Reads the existing TSV + CIF files + foldseek cigar data, performs the same
alignment as stage3_extract_align.py, and measures the distance from each
coordinatable substitution's nearest coordinating atom to the canonical Cu
position in the B2ZB02 frame."""

import argparse, csv, os, re, math
from collections import defaultdict
import numpy as np

AA3 = {'ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS',
       'MET','PHE','PRO','SER','THR','TRP','TYR','VAL'}
IMIDAZOLE = {'ND1','NE2','CD2','CE1'}
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

def ca_map(atoms):
    m={}
    for a in atoms:
        if a['atom']=='CA' and a['seq'].isdigit(): m[int(a['seq'])]=np.array([a['x'],a['y'],a['z']])
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

def build_ref_anchor_sites(ref):
    cus=get_cu(ref); ref_ca=ca_map(ref)
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
    g0min=min(groups[0]); g1min=min(groups[1])
    cuA_g,cuB_g=(0,1) if g0min<g1min else (1,0)
    return cus,ref_ca,sorted(groups[cuA_g]),sorted(groups[cuB_g])


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--cif-dir',required=True)
    ap.add_argument('--pmtyr',required=True)
    ap.add_argument('--fs-tsv',required=True)
    ap.add_argument('--input-tsv',required=True)
    a=ap.parse_args()

    ref=parse_atoms(a.pmtyr)
    ref_cus,ref_ca,cuA_pos,cuB_pos=build_ref_anchor_sites(ref)
    ref_seqs=ordered_seqs(ref); rcs=ref_core_seqs(ref_ca)

    imid=defaultdict(list)
    for at in ref:
        if at['resn']=='HIS' and at['atom'] in IMIDAZOLE and at['seq'].isdigit():
            imid[int(at['seq'])].append(at)
    cuA_to0=sum(min(d3(x,ref_cus[0]) for x in imid[p]) for p in cuA_pos)
    cuA_to1=sum(min(d3(x,ref_cus[1]) for x in imid[p]) for p in cuA_pos)
    ref_CuA=vec(ref_cus[0] if cuA_to0<cuA_to1 else ref_cus[1])
    ref_CuB=vec(ref_cus[1] if cuA_to0<cuA_to1 else ref_cus[0])
    print(f"Ref CuA: {ref_CuA}")
    print(f"Ref CuB: {ref_CuB}")

    POS_LABELS=[('CuA_His1',cuA_pos[0]),('CuA_His2',cuA_pos[1]),('CuA_His3',cuA_pos[2]),
                ('CuB_His1',cuB_pos[0]),('CuB_His2',cuB_pos[1]),('CuB_His3',cuB_pos[2])]

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

    acc2cif={}
    for fn in os.listdir(a.cif_dir):
        if fn.endswith('.cif'):
            acc2cif[fn.split('_taxID_')[0]]=fn
    print(f"{len(acc2cif)} CIFs indexed")

    with open(a.input_tsv) as f:
        rows=list(csv.DictReader(f, delimiter='\t'))
    print(f"{len(rows)} rows in input TSV")

    row_by_acc={r['accession']:r for r in rows}
    need=set()
    for r in rows:
        for lab,_ in POS_LABELS:
            if r.get(lab,'') in COORD_ATOMS:
                need.add(r['accession']); break
    print(f"{len(need)} structures have coordinatable substitutions")

    coord_dists={}
    n_done=0
    for acc in sorted(need):
        fn=acc2cif.get(acc)
        if not fn or acc not in fs: continue
        qa=parse_atoms(os.path.join(a.cif_dir,fn))
        qca=ca_map(qa); qseqs=ordered_seqs(qa); qaxyz=atom_xyz_map(qa)
        qres3=resn3_map(qa)
        qstart,tstart,cigar,_=fs[acc]
        seed=helix_seed(qca,qseqs,qstart,tstart,cigar,ref_ca,ref_seqs)
        if seed is None: continue
        R,t,_=seed
        R,t,_,_=icp_refine(qca,R,t,ref_ca,rcs)
        qca_t={s:R@c+t for s,c in qca.items()}

        row=row_by_acc[acc]
        for lab,p in POS_LABELS:
            r3=row.get(lab,'')
            if r3 not in COORD_ATOMS: continue
            rc=ref_ca[p]; bd=999.0; bs=None
            for s,c in qca_t.items():
                dd=float(np.linalg.norm(c-rc))
                if dd<bd: bd=dd; bs=s
            if bs is None or bd>3.0: continue
            if qres3.get(bs,'')!=r3: continue
            cu_ref=ref_CuA if lab.startswith('CuA') else ref_CuB
            best_cd=None
            for aname in COORD_ATOMS[r3]:
                if aname in qaxyz[bs]:
                    coord_t=R@qaxyz[bs][aname]+t
                    cd=float(np.linalg.norm(coord_t-cu_ref))
                    if best_cd is None or cd<best_cd: best_cd=cd
            if best_cd is not None:
                coord_dists[(acc,lab)]=f"{best_cd:.2f}"
        n_done+=1
        if n_done%100==0: print(f"  [{n_done}/{len(need)}]")

    print(f"Computed {len(coord_dists)} coord_dist values from {n_done} structures")

    new_fields=list(rows[0].keys())
    for lab,_ in POS_LABELS:
        col=f"{lab}_coord_dist"
        if col not in new_fields: new_fields.append(col)

    for r in rows:
        acc=r['accession']
        for lab,_ in POS_LABELS:
            r[f"{lab}_coord_dist"]=coord_dists.get((acc,lab),'')

    with open(a.input_tsv,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=new_fields,delimiter='\t')
        w.writeheader(); w.writerows(rows)

    n_filled=sum(1 for v in coord_dists.values() if v)
    n_over4=sum(1 for v in coord_dists.values() if v and float(v)>4.0)
    print(f"\nAugmented {a.input_tsv}: {n_filled} coord_dist values, {n_over4} exceed 4.0 A")

if __name__=='__main__': main()
