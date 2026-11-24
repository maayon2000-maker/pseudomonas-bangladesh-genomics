from Bio import Phylo
import pandas as pd, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
BG='#fcfcfb'; INK='#0b0b0b'; INK2='#52514e'
COL={'Bangladesh':'#2a78d6','India':'#eb6834','Pakistan':'#1baf7a'}
m=pd.read_csv('tree_meta.csv'); m['st']=m.st.astype(str).replace({'-':'n/a'}); meta=m.set_index('asm_acc')
t=Phylo.read('pa_tree.treefile','newick'); t.root_at_midpoint(); t.ladderize()
tips=t.get_terminals(); ypos={c:i for i,c in enumerate(tips)}
CLIP=0.10
def setx(cl,x):
    cl._x=x
    for c in cl.clades:
        bl=c.branch_length or 0; c._clipped=bl>CLIP; setx(c,x+min(bl,CLIP))
setx(t.root,0)
def sety(cl):
    if cl.is_terminal(): cl._y=ypos[cl]
    else:
        for c in cl.clades: sety(c)
        ys=[c._y for c in cl.clades]; cl._y=(min(ys)+max(ys))/2
sety(t.root)
xmax=max(c._x for c in tips)
fig,ax=plt.subplots(figsize=(10.5,20),dpi=120); fig.patch.set_facecolor(BG); ax.set_facecolor(BG)
xc=xmax+0.012
def draw(cl):
    for c in cl.clades:
        ax.plot([cl._x,cl._x],[cl._y,c._y],color='#7a7872',lw=0.8,solid_capstyle='butt')
        ax.plot([cl._x,c._x],[c._y,c._y],color='#7a7872',lw=0.8,solid_capstyle='butt')
        if c._clipped:
            xm=(cl._x+c._x)/2; ax.plot([xm-0.004,xm+0.004],[c._y-0.45,c._y+0.45],color=BG,lw=3,zorder=3)
            for dx in (-0.004,0.004): ax.plot([xm+dx-0.0015,xm+dx+0.0015],[c._y-0.45,c._y+0.45],color=INK2,lw=0.8,zorder=4)
        if not c.is_terminal():
            s=c.confidence
            if s is not None and s<70: ax.scatter(c._x,c._y,s=18,facecolor=BG,edgecolor=INK2,linewidth=0.8,zorder=5)
        draw(c)
draw(t.root)
for c in tips:
    nm=c.name; y=c._y
    if nm=='Reference':
        ax.plot([c._x,xc],[y,y],color='#d6d4cd',lw=0.5,ls=(0,(1,2)),zorder=1)
        ax.scatter(xc,y,s=34,facecolor=BG,edgecolor=INK2,linewidth=1.2,zorder=5)
        ax.text(xc+0.03,y,'PAO1 reference genome',va='center',fontsize=7.8,color=INK2); continue
    r=meta.loc[nm]; ax.plot([c._x,xc],[y,y],color='#d6d4cd',lw=0.5,ls=(0,(1,2)),zorder=1)
    if r.country=='Bangladesh': ax.axhspan(y-0.5,y+0.5,xmin=0,xmax=1,color='#e6effb',zorder=0,lw=0)
    ax.scatter(xc,y,s=34,color=COL[r.country],zorder=5)
    if bool(r.carb): ax.scatter(xc+0.016,y,s=30,marker='s',color=INK,zorder=5)
    ax.text(xc+0.03,y,f'ST{r.st}' if r.st!='n/a' else 'ST n/a',va='center',fontsize=7.6,color=INK)
    ax.text(xc+0.075,y,r['cat'],va='center',fontsize=7.6,color=INK2)
ax.text(xc-0.003,-2.6,'Country',ha='right',fontsize=8,color=INK2,va='center')
ax.text(xc+0.016,-2.6,'Carbapenemase\ngene',ha='center',fontsize=7.2,color=INK2,va='center')
ax.text(xc+0.03,-2.6,'Sequence type',ha='left',fontsize=8,color=INK2,va='center'); ax.text(xc+0.075,-2.6,'Source',ha='left',fontsize=8,color=INK2,va='center')
# scale bar
ax.plot([0.01,0.03],[len(tips)+0.8]*2,color=INK,lw=1.5); ax.text(0.02,len(tips)+1.8,'0.02 substitutions per site',ha='center',fontsize=7.8,color=INK2)
ax.set_ylim(len(tips)+3,-4.5); ax.set_xlim(-0.005,xc+0.14); ax.axis('off')
h=[Line2D([],[],marker='o',ls='',color=v,markersize=7,label=k) for k,v in COL.items()]
h+= [Line2D([],[],marker='s',ls='',color=INK,markersize=6,label='Predicted carbapenemase gene'),
     Line2D([],[],marker='o',ls='',markerfacecolor=BG,markeredgecolor=INK2,markersize=6,label='Branch support below 70%'),
     Line2D([],[],marker='s',ls='',color='#e6effb',markersize=9,label='Bangladeshi genome (row shading)')]
fig.legend(handles=h,loc='lower center',ncol=3,frameon=False,fontsize=8.5,bbox_to_anchor=(0.5,0.005))
fig.suptitle('Core-SNP maximum-likelihood tree of 122 P. aeruginosa genomes (60 from Bangladesh)',x=0.02,ha='left',fontsize=12.5,fontweight='bold',color=INK,y=0.992)
fig.text(0.02,0.972,'IQ-TREE, GTR+F+ASC+G4, 228,715 variable core sites vs PAO1, 1,000 ultrafast bootstraps. Midpoint-rooted for display only. Two long branches shortened (//).',fontsize=8.3,color=INK2)
fig.subplots_adjust(left=0.01,right=0.99,top=0.955,bottom=0.045); fig.savefig('fig15_tree.png',facecolor=BG)
print(xmax)
