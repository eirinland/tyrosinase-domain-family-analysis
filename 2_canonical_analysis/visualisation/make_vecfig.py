import pandas as pd, numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.spatial.distance import pdist, squareform
from sklearn.manifold import TSNE

BASE="/sessions/intelligent-blissful-knuth/mnt/2_canonical_analysis"
NODES=f"{BASE}/visualisation/network_nodes.csv"
ALLV=f"{BASE}/all_vectors.csv"
OUT="/sessions/intelligent-blissful-knuth/mnt/outputs"

POS=['Gly46','Phe65','Trp68','Glu195','Asn205','Arg209','Val218','Ala221','Phe227','His230','thioether']

df=pd.read_csv(NODES)
for c in POS: df[c]=df[c].astype(str)
df['is_char']=df['characterized'].astype(str).str.lower().eq('yes')
print("nodes:",len(df)," characterized:",df['is_char'].sum())

def hamming_int(frame, cols):
    Mint=np.zeros((len(frame),len(cols)),dtype=int)
    for j,c in enumerate(cols):
        vals=frame[c].astype(str).values
        cats={v:i for i,v in enumerate(pd.unique(vals))}
        Mint[:,j]=[cats[v] for v in vals]
    return squareform(pdist(Mint, metric='hamming'))

def pcoa(D,k=2):
    n=D.shape[0]; D2=D**2
    J=np.eye(n)-np.ones((n,n))/n
    B=-0.5*J@D2@J
    w,V=np.linalg.eigh(B)
    idx=np.argsort(w)[::-1]; w=w[idx]; V=V[:,idx]
    L=np.sqrt(np.clip(w[:k],0,None))
    return V[:,:k]*L

D=hamming_int(df,POS)
XY_mds=pcoa(D,2)
try:
    XY_tsne=TSNE(n_components=2,metric='precomputed',init='random',perplexity=30,random_state=0).fit_transform(D)
except Exception as e:
    print("tsne fail:",e); XY_tsne=XY_mds

ACT=['TYR','CaOx','oMP','oAPO','AUS','DHICA ox','DCT','hemocyanin']
ACT_COLORS={'TYR':'#0072B2','CaOx':'#E69F00','oMP':'#009E73','oAPO':'#CC79A7',
            'AUS':'#C7B200','DHICA ox':'#D55E00','DCT':'#444444','hemocyanin':'#56B4E9'}
KING_COLORS={'Fungi':'#1b9e77','Animals':'#d95f02','Plants':'#7570b3',
             'Bacteria':'#e7298a','Oomycota':'#66a61e','?':'#bdbdbd'}

def plot_activity(ax,XY,title):
    bg=(~df['is_char']).values
    ax.scatter(XY[bg,0],XY[bg,1],s=np.sqrt(df.loc[bg,'count'])*1.3,c='#D5D5D5',
               alpha=0.55,linewidths=0,rasterized=True)
    present=[a for a in ACT if a in set(df.loc[df['is_char'],'activity'])]
    for act in present:
        grp=df[df['is_char'] & (df['activity']==act)]
        ii=grp.index.values
        ax.scatter(XY[ii,0],XY[ii,1],s=np.sqrt(grp['count'])*3+30,
                   c=ACT_COLORS.get(act,'#333'),edgecolors='black',linewidths=0.6,
                   label=f'{act} ({len(grp)})',zorder=5)
    ax.set_title(title,fontsize=11)
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values(): sp.set_visible(False)
    return present

def plot_kingdom(ax,XY,title):
    order=['Fungi','Animals','Plants','Bacteria','Oomycota','?']
    for k in order:
        grp=df[df['top_kingdom']==k]
        if len(grp)==0: continue
        ii=grp.index.values
        ax.scatter(XY[ii,0],XY[ii,1],s=np.sqrt(grp['count'])*1.5,
                   c=KING_COLORS.get(k,'#999'),alpha=0.6,linewidths=0,
                   label=f'{k}',rasterized=True)
    ax.set_title(title,fontsize=11); ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values(): sp.set_visible(False)

# ---------- class-collapsed ----------
CLASS={}
for a in 'AVLIM': CLASS[a]='aliphatic'
for a in 'FWY': CLASS[a]='aromatic'
for a in 'STNQ': CLASS[a]='polar'
for a in 'KRH': CLASS[a]='basic'
for a in 'DE': CLASS[a]='acidic'
CLASS['C']='Cys'; CLASS['G']='Gly'; CLASS['P']='Pro'
def coll(x): return CLASS.get(x,x)
resid_pos=[p for p in POS if p!='thioether']

def collapsed_unique(frame):
    tmp=frame.copy()
    for c in resid_pos: tmp[c]=frame[c].astype(str).map(coll)
    sig=tmp[POS].astype(str).agg('|'.join,axis=1)
    return sig

# reduction stats on count>=3 set
sig_exact=df[POS].astype(str).agg('|'.join,axis=1)
sig_coll=collapsed_unique(df)
print("count>=3 set: residue-exact unique signatures:",sig_exact.nunique(),
      " class-collapsed:",sig_coll.nunique())

# reduction stats on FULL set
try:
    allv=pd.read_csv(ALLV)
    for c in POS:
        if c in allv.columns: allv[c]=allv[c].astype(str)
    sig_exact_all=allv[POS].astype(str).agg('|'.join,axis=1)
    sig_coll_all=collapsed_unique(allv)
    print("FULL set: residue-exact unique signatures:",sig_exact_all.nunique(),
          " class-collapsed:",sig_coll_all.nunique(), " total structures:",allv['count'].sum())
except Exception as e:
    print("all_vectors read issue:",e)

# build collapsed node set from count>=3 nodes (aggregate)
dfx=df.copy(); dfx['csig']=sig_coll
agg=dfx.groupby('csig').agg(count=('count','sum'),
                            is_char=('is_char','max'),
                            activity=('activity', lambda s: next((a for a,c in zip(s, dfx.loc[s.index,'is_char']) if c), 'none')),
                            top_kingdom=('top_kingdom', lambda s: s.value_counts().idxmax())).reset_index()
# representative residue-class columns for distance
rep=dfx.drop_duplicates('csig').set_index('csig')
for c in POS:
    agg[c]=agg['csig'].map(rep[c] if c=='thioether' else rep[c].map(coll))
print("collapsed nodes (count>=3):",len(agg)," characterized:",int(agg['is_char'].sum()))
Dc=hamming_int(agg,POS); XYc=pcoa(Dc,2)
aggc=agg

def plot_activity_frame(ax,frame,XY,title):
    bg=(~frame['is_char'].astype(bool)).values
    ax.scatter(XY[bg,0],XY[bg,1],s=np.sqrt(frame.loc[bg,'count'])*1.3,c='#D5D5D5',
               alpha=0.55,linewidths=0,rasterized=True)
    chv=set(frame.loc[frame['is_char'].astype(bool),'activity'])
    present=[a for a in ACT if a in chv]
    for act in present:
        grp=frame[frame['is_char'].astype(bool) & (frame['activity']==act)]
        ii=grp.index.values
        ax.scatter(XY[ii,0],XY[ii,1],s=np.sqrt(grp['count'])*3+30,
                   c=ACT_COLORS.get(act,'#333'),edgecolors='black',linewidths=0.6,
                   label=f'{act} ({len(grp)})',zorder=5)
    ax.set_title(title,fontsize=11); ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values(): sp.set_visible(False)

# ---------- COMBINED OPTIONS FIGURE ----------
fig,axs=plt.subplots(2,2,figsize=(13,12))
present=plot_activity(axs[0,0],XY_mds,f"A  PCoA / classical MDS — coloured by activity\n{len(df)} signatures (count≥3), {int(df['is_char'].sum())} characterised")
axs[0,0].legend(loc='upper left',fontsize=8,frameon=False,title='Characterised activity',title_fontsize=8)
plot_activity(axs[0,1],XY_tsne,"B  t-SNE — coloured by activity\n(same data, cleaner islands)")
axs[0,1].legend(loc='upper left',fontsize=8,frameon=False)
plot_kingdom(axs[1,0],XY_mds,"C  PCoA / MDS — coloured by kingdom\n(where the unmapped space sits)")
axs[1,0].legend(loc='upper left',fontsize=8,frameon=False,title='Top kingdom',title_fontsize=8)
plot_activity_frame(axs[1,1],aggc,XYc,f"D  Class-collapsed PCoA — coloured by activity\n{len(agg)} physicochemical-class signatures")
axs[1,1].legend(loc='upper left',fontsize=8,frameon=False)
# size legend
for s_ex,lab in [(3,'3'),(50,'50'),(300,'300')]:
    axs[0,0].scatter([],[],s=np.sqrt(s_ex)*1.3,c='#999',label=f'_{lab}')
fig.suptitle("Second-shell vector space of 21,928 canonical PPO-fold structures",fontsize=13,y=0.995)
fig.text(0.5,0.005,"Each point = one unique second-shell signature; point area ∝ number of structures. Grey = no characterised PPO with that signature. "
                   "2-D positions are an approximate projection of categorical Hamming distance.",ha='center',fontsize=8,style='italic')
plt.tight_layout(rect=[0,0.02,1,0.98])
plt.savefig(f"{OUT}/vector_space_options.png",dpi=150,bbox_inches='tight')
plt.savefig(f"{OUT}/vector_space_options.pdf",bbox_inches='tight')
print("saved options figure")

# ---------- polished MAIN 2-panel figure: shared PCoA layout ----------
fig2,(axA,axB)=plt.subplots(1,2,figsize=(14,6.6))
# Panel A: activity
present=plot_activity(axA,XY_mds,"A   Characterised activities within the realised vector space")
# annotate the largest characterised node
cc=df[df['is_char']]
big=cc.loc[cc['count'].idxmax()]
axA.annotate(f"{big['activity']} ({big['count']} structures)",
             (XY_mds[big.name,0],XY_mds[big.name,1]),
             textcoords="offset points",xytext=(8,8),fontsize=8,
             arrowprops=dict(arrowstyle='-',lw=0.5,color='#555'))
act_handles=[Line2D([],[],marker='o',ls='',mfc=ACT_COLORS[a],mec='black',mew=0.5,ms=8,
                    label=f'{a} ({len(df[df.is_char & (df.activity==a)])})') for a in present]
size_handles=[Line2D([],[],marker='o',ls='',mfc='#bbb',mec='none',ms=np.sqrt(s)*0.55,
                     label=f'{s}') for s in [5,50,300]]
leg1=axA.legend(handles=act_handles,loc='upper left',fontsize=8.5,frameon=False,
                title='Characterised activity',title_fontsize=9)
axA.add_artist(leg1)
axA.legend(handles=size_handles,loc='lower left',fontsize=8,frameon=False,
           title='structures / signature',title_fontsize=8,labelspacing=1.1,borderpad=0.8)
# Panel B: kingdom (SAME coordinates)
plot_kingdom(axB,XY_mds,"B   Taxonomy of the same space (where the unmapped diversity sits)")
king_handles=[Line2D([],[],marker='o',ls='',mfc=KING_COLORS[k],mec='none',ms=8,label=k)
              for k in ['Fungi','Animals','Plants','Bacteria','Oomycota','?'] if (df['top_kingdom']==k).any()]
axB.legend(handles=king_handles,loc='upper left',fontsize=8.5,frameon=False,
           title='Top kingdom',title_fontsize=9)
fig2.text(0.5,0.02,"Each point = one unique second-shell signature (10 superposed second-shell positions + thioether); "
                   "point area ∝ number of structures carrying it (count ≥ 3 shown). Grey = no characterised PPO with that exact signature. "
                   "Layout = classical MDS (PCoA) on categorical Hamming distance; both panels share identical coordinates.",
          ha='center',fontsize=8,style='italic',wrap=True)
plt.tight_layout(rect=[0,0.045,1,1])
plt.savefig(f"{OUT}/vector_space_main.pdf",bbox_inches='tight')
plt.savefig(f"{OUT}/vector_space_main.png",dpi=170,bbox_inches='tight')
print("saved main 2-panel figure")
print("DONE")
