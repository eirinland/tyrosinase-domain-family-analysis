#!/usr/bin/env python3
"""Reproduce stage-3 helix-anchored Kabsch+ICP superposition and write PDB output."""
import re, math, sys
import numpy as np

HELIX_RANGES = [(34,46),(65,83),(203,211),(226,244)]
ANCHORS = [42,60,69,204,208,231]
ANCHOR_LABELS = ['CuA_His1','CuA_His2','CuA_His3','CuB_His1','CuB_His2','CuB_His3']
AA3 = {'ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS',
       'MET','PHE','PRO','SER','THR','TRP','TYR','VAL'}

def parse_cif_atoms(path):
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
                    atoms.append({'group':r.get('_atom_site.group_PDB','ATOM'),
                        'serial':int(r.get('_atom_site.id','0')),
                        'elem':r['_atom_site.type_symbol'].upper(),
                        'atom':r['_atom_site.label_atom_id'].upper(),
                        'resn':r['_atom_site.label_comp_id'].upper(),
                        'chain':r.get('_atom_site.label_asym_id','A'),
                        'seq':r['_atom_site.label_seq_id'],
                        'bfac':float(r.get('_atom_site.B_iso_or_equiv','0') or 0),
                        'x':float(r['_atom_site.Cartn_x']),
                        'y':float(r['_atom_site.Cartn_y']),
                        'z':float(r['_atom_site.Cartn_z'])})
                except (KeyError,ValueError): continue
    return atoms

def vec(a): return np.array([a['x'],a['y'],a['z']])
def ca_map(atoms):
    m={}
    for a in atoms:
        if a['atom']=='CA' and a['seq'].isdigit(): m[int(a['seq'])]=vec(a)
    return m
def resn3_map(atoms):
    m={}
    for a in atoms:
        if a['seq'].isdigit() and a['resn'] in AA3:
            s=int(a['seq'])
            if s not in m: m[s]=a['resn']
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
    for _ in range(iters):
        Qt=(R@Qbase.T).T+t
        P=[]; Q=[]
        for rc in ref_core:
            dd=np.linalg.norm(Qt-rc,axis=1); j=int(dd.argmin())
            if dd[j]<=cutoff: P.append(Qbase[j]); Q.append(rc)
        if len(P)<8: break
        R,t=kabsch(np.array(P),np.array(Q))
    return R,t

def write_pdb(atoms, R, t, outpath):
    with open(outpath,'w') as f:
        for a in atoms:
            v = R @ vec(a) + t
            group = 'HETATM' if a['group']=='HETATM' or a['resn'] not in AA3|{'HOH','CU','ZN','FE','MG','CA','NA','CL'} else 'ATOM'
            if a['resn'] not in AA3: group = 'HETATM'
            seq = int(a['seq']) if a['seq'].isdigit() else 0
            f.write(f"{group:<6s}{a['serial']:>5d} {a['atom']:<4s} {a['resn']:>3s} {a['chain']:>1s}{seq:>4d}    "
                    f"{v[0]:>8.3f}{v[1]:>8.3f}{v[2]:>8.3f}{1.0:>6.2f}{a['bfac']:>6.2f}          {a['elem']:>2s}\n")
        f.write('END\n')

def main():
    base = sys.argv[0].rsplit('/',1)[0] if '/' in sys.argv[0] else '.'

    ref_cif = f"{base}/PmTYR_B2ZB02.cif"
    query_cif = f"{base}/Oomycota_3His_H6__A0ABR1FRY0.cif"

    cigar_str = "23M9I6M4I12M1D4M4D42M1I1M8I7M4I16M2D20M2I4M12I5M2I7M3I2M1D2M39I12M23I7M16I9M6I11M5I3M3I3M2I3M4I4M1D20M2I7M2I2M4I2M10I4M1D14M1I2M12I13M5I4M"
    qstart, tstart = 143, 6

    print("Parsing reference B2ZB02...")
    ref = parse_cif_atoms(ref_cif)
    ref_ca = ca_map(ref)
    ref_seqs = ordered_seqs(ref)
    rcs = [s for s in ref_ca if in_helix(s)]
    ref_res3 = resn3_map(ref)

    print("Parsing query A0ABR1FRY0...")
    query = parse_cif_atoms(query_cif)
    qca = ca_map(query)
    qseqs = ordered_seqs(query)
    qres3 = resn3_map(query)

    print(f"Reference: {len(ref_ca)} CA, core helix CA: {len(rcs)}")
    print(f"Query: {len(qca)} CA")

    print("\nHelix-anchored Kabsch seed...")
    seed = helix_seed(qca, qseqs, qstart, tstart, cigar_str, ref_ca, ref_seqs)
    if seed is None:
        print("ERROR: helix seed failed (<8 correspondences)")
        return
    R, t, n_seed = seed
    print(f"  {n_seed} helix CA pairs used for Kabsch seed")

    print("ICP refinement...")
    R, t = icp_refine(qca, R, t, ref_ca, rcs)

    qca_t = {s: R @ c + t for s, c in qca.items()}

    print("\nPosition mapping (3.0 A cutoff):")
    for label, anc in zip(ANCHOR_LABELS, ANCHORS):
        rc = ref_ca[anc]
        bd = 999.0; bs = None
        for s, c in qca_t.items():
            dd = float(np.linalg.norm(c - rc))
            if dd < bd: bd = dd; bs = s
        mapped = qres3.get(bs, 'UNK') if bs and bd <= 3.0 else '---'
        ref_r = ref_res3.get(anc, '?')
        print(f"  {label:12s} ref={ref_r}{anc:<4d} -> query={mapped}{f' ({bs}, {bd:.2f}A)' if bs else ''}")

    out_query = f"{base}/A0ABR1FRY0_stage3_superimposed.pdb"
    print(f"\nWriting superimposed query -> {out_query}")
    write_pdb(query, R, t, out_query)

    out_ref = f"{base}/B2ZB02_reference_frame.pdb"
    print(f"Writing reference (identity) -> {out_ref}")
    I = np.eye(3)
    z = np.zeros(3)
    write_pdb(ref, I, z, out_ref)

    print("\nDone. Load both PDB files in PyMOL to see the stage-3 superposition.")
    print("The query is in the EXACT frame that stage-3 uses for nearest-CA lookup.")

if __name__ == '__main__':
    main()
