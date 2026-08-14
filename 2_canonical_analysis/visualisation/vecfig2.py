"""Vector-space landscape of all unique second-shell signatures (PCoA).
Panel A: characterised activity (PyMOL palette); panel B: taxonomic lineage.
Also writes a hotspot-highlight figure for the divergent novel groups.
Reads position_vectors.csv and characterized_PPOs.xlsx (one level up) and
visualisation/taxonomy_lookup.csv; writes the figures next to this script."""
import csv, math
from pathlib import Path
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.spatial.distance import pdist, squareform
from scipy.sparse.linalg import eigsh
from collections import Counter, defaultdict

BASE=Path(__file__).resolve().parent.parent
OUT=Path(__file__).resolve().parent
POSITIONS=['Gly46','Phe65','Trp68','Glu195','Asn205','Arg209','Val218','Ala221','Phe227','His230','thioether']

import pandas as pd
_cdf=pd.read_excel(f"{BASE}/characterized_PPOs.xlsx")
CHARACTERIZED=dict(zip(_cdf['Accession'].astype(str).str.strip(), _cdf['Activity'].astype(str).str.strip()))

# Activity palette (PyMOL set_color values, matches the structure renders)
ACT_COLOR={'TYR':(0.30,0.55,0.76),'CaOx':(0.53,0.71,0.82),'AUS':(0.62,0.80,0.78),'oAPO':(0.65,0.82,0.65),
           'oMP':(0.42,0.67,0.46),'DHICA ox':(0.94,0.81,0.69),'DCT':(0.85,0.65,0.54),'hemocyanin':(0.71,0.62,0.78)}
ACT_ORDER=['TYR','CaOx','AUS','oAPO','oMP','DHICA ox','DCT','hemocyanin']

def get_parts(row):
    parts=[]
    for pos in POSITIONS:
        val=row.get(pos,'?') or '?'
        if pos=='thioether':
            val='C' if val in ('C','C*') else '-'
        elif pos=='Gly46' and row.get('Gly46_ss')=='c' and val!='?':
            val='@'
        parts.append(val)
    return tuple(parts)

# taxonomy
tax={}
with open(f"{BASE}/visualisation/taxonomy_lookup.csv") as f:
    for r in csv.DictReader(f): tax[r['accession']]=r.get('kingdom','?')

vec_count=Counter(); vec_parts={}; vec_acts=defaultdict(set); vec_king=defaultdict(Counter)
with open(f"{BASE}/position_vectors.csv") as f:
    for row in csv.DictReader(f):
        p=get_parts(row); vec_count[p]+=1; vec_parts[p]=p
        acc=row['accession']
        if acc in CHARACTERIZED: vec_acts[p].add(CHARACTERIZED[acc])
        vec_king[p][tax.get(acc,'?')]+=1

vectors=sorted(vec_count.keys())            # ALL unique vectors, no filter
n=len(vectors)
print("total structures:",sum(vec_count.values()),"| ALL unique vectors:",n)
print("characterised unique vectors:",sum(1 for v in vectors if v in vec_acts))

# integer-encode each position for hamming
Mint=np.zeros((n,len(POSITIONS)),dtype=int)
for j in range(len(POSITIONS)):
    col=[v[j] for v in vectors]; cats={c:i for i,c in enumerate(sorted(set(col)))}
    Mint[:,j]=[cats[c] for c in col]
D=squareform(pdist(Mint,metric='hamming'))   # fraction of 11 positions differing
# classical MDS (PCoA) via double-centering + top-2 eigvecs
D2=D**2
B=-0.5*(D2 - D2.mean(0,keepdims=True) - D2.mean(1,keepdims=True) + D2.mean())
w,V=eigsh(B,k=2,which='LA'); idx=np.argsort(w)[::-1]; w=w[idx]; V=V[:,idx]
XY=V*np.sqrt(np.clip(w,0,None))
print("top-2 eigenvalue share of positive spread: %.1f%%"%(100*np.clip(w,0,None).sum()/np.clip(np.linalg.eigvalsh(B),0,None).sum()))

counts=np.array([vec_count[v] for v in vectors],float)
logc=np.log10(counts)
psize=6+(logc-logc.min())/(logc.max()-logc.min())*70
is_char=np.array([v in vec_acts for v in vectors])
act_of=[sorted(vec_acts[v])[0] if v in vec_acts else None for v in vectors]

KING_COLOR={'Fungi':(0.95,0.89,0.07),'Metazoa':(0.25,0.55,0.83),'Bacteria':(0.34,0.74,0.88),
            'Viridiplantae':(0.13,0.53,0.30),'Other Eukaryota':(0.52,0.24,0.90),'Unknown':(0.62,0.62,0.62)}
KING_MAP={'Fungi':'Fungi','Animals':'Metazoa','Bacteria':'Bacteria','Plants':'Viridiplantae','Oomycota':'Other Eukaryota','Archaea':'Unknown','?':'Unknown'}
KING_ORDER=['Fungi','Metazoa','Bacteria','Viridiplantae','Other Eukaryota','Unknown']
topk=[KING_MAP.get(vec_king[v].most_common(1)[0][0],'Unknown') for v in vectors]

def panel_activity(ax):
    bg=~is_char
    ax.scatter(XY[bg,0],XY[bg,1],s=psize[bg],c='#EAEAEA',alpha=0.5,linewidths=0,rasterized=True)
    present=[a for a in ACT_ORDER if a in set(x for x in act_of if x)]
    for a in present:
        m=np.array([x==a for x in act_of])
        ax.scatter(XY[m,0],XY[m,1],s=psize[m]*2.2+45,color=ACT_COLOR[a],marker='D',
                   edgecolors='black',linewidths=1.0,zorder=5,label=f'{a} ({m.sum()})')
    ax.set_title("A   Characterised activities in the second-shell vector space",fontsize=11)
    ax.set_xticks([]);ax.set_yticks([])
    for s in ax.spines.values(): s.set_visible(False)
    return present

def panel_kingdom(ax):
    for k in KING_ORDER:
        m=np.array([t==k for t in topk])
        if m.sum()==0: continue
        ax.scatter(XY[m,0],XY[m,1],s=psize[m],color=KING_COLOR[k],alpha=0.7,linewidths=0,rasterized=True,label=f'{k}')
    ax.set_title("B   Taxonomy of the same space",fontsize=11)
    ax.set_xticks([]);ax.set_yticks([])
    for s in ax.spines.values(): s.set_visible(False)

fig,(axA,axB)=plt.subplots(1,2,figsize=(14,6.6))
present=panel_activity(axA)
ah=[Line2D([],[],marker='D',ls='',mfc=ACT_COLOR[a],mec='black',mew=0.6,ms=8,label=a) for a in present]
sh=[Line2D([],[],marker='o',ls='',mfc='#b0b0b0',mec='none',ms=np.sqrt(6+(math.log10(s)-logc.min())/(logc.max()-logc.min())*70)*0.7,label=f'{s}') for s in [1,30,300]]
l1=axA.legend(handles=ah,loc='upper left',fontsize=8.5,frameon=False,title='Characterised activity',title_fontsize=9);axA.add_artist(l1)
axA.legend(handles=sh,loc='lower left',fontsize=8,frameon=False,title='structures / signature',title_fontsize=8,labelspacing=1.0,borderpad=0.7)
panel_kingdom(axB)
kh=[Line2D([],[],marker='o',ls='',mfc=KING_COLOR[k],mec='none',ms=8,label=k) for k in KING_ORDER]
axB.legend(handles=kh,loc='upper left',fontsize=8.5,frameon=False,title='Lineage',title_fontsize=9)
fig.suptitle(f"Second-shell vector space of the canonical PPO fold — all {n:,} unique signatures",fontsize=13,y=0.99)
fig.text(0.5,0.015,"Each point = one unique second-shell signature (10 positions + thioether); area ∝ structures carrying it. "
         "Grey = no structurally characterised PPO with that signature. Layout = classical MDS (PCoA) on Hamming distance; both panels share coordinates.",
         ha='center',fontsize=8,style='italic')
plt.tight_layout(rect=[0,0.04,1,0.97])
plt.savefig(f"{OUT}/vector_space_all_main.pdf",bbox_inches='tight')
plt.savefig(f"{OUT}/vector_space_all_main.png",dpi=170,bbox_inches='tight')
print("saved vector_space_all_main")

# ---------- HOTSPOT HIGHLIGHT figure ----------
iH=POSITIONS.index('His230'); iP=POSITIONS.index('Phe227'); iN=POSITIONS.index('Asn205'); iR=POSITIONS.index('Arg209'); iV=POSITIONS.index('Val218')
g119=np.array([v[iH]=='Y' and v[iN]=='D' and v[iR]=='G' and v[iV]=='N' for v in vectors])
gpV=np.array([v[iP]=='V' for v in vectors])
gpW=np.array([v[iP]=='W' for v in vectors])
def nstr(mask): return int(counts[mask].sum())
print("hotspot nodes: 119-core=%d (%d struct); Phe227=Val=%d (%d); Phe227=Trp=%d (%d)"%(
    g119.sum(),nstr(g119),gpV.sum(),nstr(gpV),gpW.sum(),nstr(gpW)))
figH,axH=plt.subplots(figsize=(7.8,7.2))
axH.scatter(XY[:,0],XY[:,1],s=psize,c='#E7E7E7',alpha=0.5,linewidths=0,rasterized=True)
axH.scatter(XY[is_char,0],XY[is_char,1],s=34,facecolors='none',edgecolors='#555555',
            linewidths=0.7,marker='D',zorder=4,label='characterised signatures (53)')
HL=[(g119,'His230=Tyr core (Fusarium / Pseudomonas)','#D7263D','o'),
    (gpV,'Phe227=Val (actinobacteria + ascomycetes)','#F46036','s'),
    (gpW,'Phe227=Trp (black yeasts)','#8E24AA','^')]
for mask,lab,col,mk in HL:
    axH.scatter(XY[mask,0],XY[mask,1],s=np.sqrt(counts[mask])*8+35,color=col,marker=mk,
                edgecolors='black',linewidths=0.6,alpha=0.92,zorder=6,label=f'{lab}  (n={nstr(mask)})')
axH.set_title("Divergent novel-fraction groups in the vector space",fontsize=12)
axH.set_xticks([]); axH.set_yticks([])
for s in axH.spines.values(): s.set_visible(False)
axH.legend(loc='upper left',fontsize=8.5,frameon=False)
figH.text(0.5,0.02,"Same PCoA layout as the main figure. Grey = all 3,897 signatures; open diamonds = the 53 characterised signatures. "
          "Point area is proportional to the number of structures carrying each signature.",ha='center',fontsize=8,style='italic')
plt.tight_layout(rect=[0,0.04,1,1])
plt.savefig(f"{OUT}/vector_space_hotspots.pdf",bbox_inches='tight')
plt.savefig(f"{OUT}/vector_space_hotspots.png",dpi=175,bbox_inches='tight')
print("saved vector_space_hotspots")
