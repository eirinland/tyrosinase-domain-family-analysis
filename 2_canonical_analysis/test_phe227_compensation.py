"""
Test whether the Phe227->Val substitution is structurally compensated by
another phenylalanine occupying the vacated aromatic cavity.

Method (alignment-free, Cu-anchored): each aromatic ring is described by its
distances to the two Cu atoms, stored sorted (dmin,dmax) so the descriptor is
invariant to Cu labelling and to overall orientation. We locate the normal
Phe227 ring slot from PmTYR (B2ZB02, Phe227=F) and from the control set, then
ask whether Phe227=Val structures have a Phe ring in that same slot.
"""
import glob, os, math, statistics
import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'electrostatics'))
from extract_features import parse_cif, dist

RING={'PHE':['CG','CD1','CD2','CE1','CE2','CZ'],'TYR':['CG','CD1','CD2','CE1','CE2','CZ'],
      'TRP':['CG','CD1','CD2','NE1','CE2','CE3','CZ2','CZ3','CH2']}

def cu_pair(atoms):
    cus=[(a['x'],a['y'],a['z']) for a in atoms if a['elem']=='CU']
    if len(cus)<2: return None
    best=None
    for i in range(len(cus)):
        for j in range(i+1,len(cus)):
            d=dist(cus[i],cus[j])
            if best is None or d<best[0]: best=(d,cus[i],cus[j])
    return best[1],best[2]

def aromatics(atoms):
    """Return list of (restype, (dmin,dmax) placeholder filled later, centroid)."""
    groups={}
    for a in atoms:
        if a['res'] in RING and a['name'] in RING[a['res']]:
            groups.setdefault((a['chain'],a['seq'],a['res']),[]).append(a)
    out=[]
    for (ch,seq,res),ats in groups.items():
        if len(ats)<4: continue
        c=(sum(x['x'] for x in ats)/len(ats),sum(x['y'] for x in ats)/len(ats),sum(x['z'] for x in ats)/len(ats))
        out.append((res,c))
    return out

def descr(path):
    atoms=parse_cif(path)
    cp=cu_pair(atoms)
    if not cp: return None
    cu1,cu2=cp
    mid=tuple((cu1[k]+cu2[k])/2 for k in range(3))
    res=[]
    for restype,c in aromatics(atoms):
        if dist(mid,c)<=12:
            d1,d2=dist(cu1,c),dist(cu2,c)
            res.append((restype,min(d1,d2),max(d1,d2)))
    return res

def load(listfile):
    accs=[l.strip() for l in open(listfile) if l.strip()]
    return set(accs)

valset=load('/tmp/comp_val.txt'); ctrlset=load('/tmp/comp_ctrl.txt')
cifs=glob.glob('/tmp/comp_cifs/*.cif')
V=[]; F=[]
for p in cifs:
    acc=os.path.basename(p).split('_taxID_')[0]
    d=descr(p)
    if d is None: continue
    (V if acc in valset else F).append((acc,d))
print(f"parsed Val={len(V)}  ctrl(Phe227=F)={len(F)}")

# reference Phe227 slot from B2ZB02
ref=descr(glob.glob('cifs/B2ZB02_taxID_*_model.cif')[0]) if glob.glob('cifs/B2ZB02_taxID_*_model.cif') else None
if ref:
    phes=sorted([r for r in ref if r[0]=='PHE'],key=lambda x:x[1])
    print("B2ZB02 near-site PHE rings (restype,dmin,dmax):",[(round(a[1],1),round(a[2],1)) for a in phes])

def count_near(desclist, slot, tol=1.6, types=('PHE',)):
    """fraction of structures with an aromatic of given types within tol of slot (dmin,dmax)."""
    s_dmin,s_dmax=slot
    hits=0
    for acc,d in desclist:
        ok=any(r[0] in types and abs(r[1]-s_dmin)<=tol and abs(r[2]-s_dmax)<=tol for r in d)
        if ok: hits+=1
    return 100*hits/len(desclist)

# Identify the two conserved near-site Phe slots from the CONTROL set (Phe227=F):
# bin all control PHE (dmin,dmax), find the two densest clusters.
import numpy as np
allF=[(r[1],r[2]) for acc,d in F for r in d if r[0]=='PHE']
arr=np.array(allF)
# crude cluster: round to 1A grid, find top cells
from collections import Counter
grid=Counter((round(x[0]),round(x[1])) for x in allF)
print("\nTop control(Phe227=F) PHE (dmin,dmax) clusters [rounded]:",grid.most_common(5))

print("\nFor each candidate Phe slot: % of structures with a PHE there, and with ANY aromatic there")
for slot,_ in grid.most_common(4):
    fF=count_near(F,slot,types=('PHE',)); fV=count_near(V,slot,types=('PHE',))
    aF=count_near(F,slot,types=('PHE','TYR','TRP')); aV=count_near(V,slot,types=('PHE','TYR','TRP'))
    print(f"  slot {slot}: PHE  F={fF:.0f}% V={fV:.0f}% | ANY-aromatic F={aF:.0f}% V={aV:.0f}%")

# scalar summary: mean # PHE and # aromatics within 9A of midpoint
def near_counts(desclist,rmax=9.0):
    nphe=[];narom=[]
    for acc,d in desclist:
        nphe.append(sum(1 for r in d if r[0]=='PHE' and r[2]<=rmax))
        narom.append(sum(1 for r in d if r[2]<=rmax))
    return statistics.mean(nphe),statistics.mean(narom)
fp,fa=near_counts(F); vp,va=near_counts(V)
print(f"\nMean # PHE rings within 9A of Cu-midpoint:  F(227=Phe)={fp:.2f}  V(227=Val)={vp:.2f}  (diff={fp-vp:+.2f})")
print(f"Mean # aromatic rings within 9A:            F={fa:.2f}  V={va:.2f}  (diff={fa-va:+.2f})")
