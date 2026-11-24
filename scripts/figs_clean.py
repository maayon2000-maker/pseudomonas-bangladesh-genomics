import pandas as pd, numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from math import sqrt
BG='#fcfcfb'; INK='#0b0b0b'; INK2='#52514e'; COL={'Bangladesh':'#2a78d6','India':'#eb6834','Pakistan':'#1baf7a'}
d=pd.read_csv('clusters_main_clean.csv')
def w(k,n,z=1.96):
    p=k/n; den=1+z*z/n; c=(p+z*z/(2*n))/den; h=z*sqrt(p*(1-p)/n+z*z/(4*n*n))/den; return 100*(c-h),100*(c+h)
dr=d.sort_values('n50',ascending=False).drop_duplicates(['cl_99.99','country'])
fig,axs=plt.subplots(1,2,figsize=(11,5),dpi=160,sharey=True); fig.patch.set_facecolor(BG)
for ax,(df,t) in zip(axs,[(d,'All genomes'),(dr,'One genome per country per clone group\n(ANI ≥ 99.99%)')]):
    ax.set_facecolor(BG)
    for i,c in enumerate(['Bangladesh','India','Pakistan']):
        s=df[df.country==c]; k=int(s.carb.sum()); n=len(s); lo,hi=w(k,n); p=100*k/n
        ax.bar(i,p,color=COL[c],width=0.55,zorder=2); ax.plot([i,i],[lo,hi],color=INK,lw=1.4,zorder=3)
        ax.text(i+0.34,p,f'{k}/{n}\n{p:.1f}%',va='center',fontsize=9,color=INK)
    ax.set_xticks(range(3)); ax.set_xticklabels(['Bangladesh','India','Pakistan'],fontsize=10,color=INK)
    ax.set_title(t,fontsize=10.5,color=INK,loc='left'); ax.grid(axis='y',color='#ecebe6',zorder=0)
    for s_ in ['top','right','left']: ax.spines[s_].set_visible(False)
    ax.spines['bottom'].set_color('#cfcdc6'); ax.tick_params(axis='y',length=0,labelcolor=INK2); ax.set_xlim(-0.6,2.9)
axs[0].set_ylabel('Genomes with a predicted carbapenemase gene (%)',fontsize=10,color=INK2); axs[0].set_ylim(0,45)
fig.suptitle('Predicted carbapenemase carriage by country, laboratory strains removed',x=0.01,ha='left',fontsize=12.5,fontweight='bold',color=INK,y=0.985)
fig.text(0.01,0.895,'Bars: proportion with a carbapenemase gene (AMRFinderPlus); lines: 95% Wilson interval. Genes are predicted, phenotype not tested.',fontsize=8.8,color=INK2)
fig.legend(handles=[Patch(color=v,label=k) for k,v in COL.items()],loc='lower center',ncol=3,frameon=False,fontsize=9)
fig.subplots_adjust(top=0.76,bottom=0.17,left=0.07,right=0.98,wspace=0.08); fig.savefig('fig13_carb_clean.png',facecolor=BG)
# fig14 sharing
rows=[]
for col,lab in [('cl_99.9','ANI ≥ 99.9%'),('cl_99.95','ANI ≥ 99.95%'),('cl_99.99','ANI ≥ 99.99%')]:
    g=d.groupby(col).country.agg(set); bd=d[d.country=='Bangladesh']
    I=sum('India' in g[x] for x in bd[col]); P=sum('Pakistan' in g[x] for x in bd[col]); A=sum(len(g[x]-{'Bangladesh'})>0 for x in bd[col]); rows.append((lab,I,P,A,len(bd)))
fig,ax=plt.subplots(figsize=(9,5),dpi=160); fig.patch.set_facecolor(BG); ax.set_facecolor(BG)
x=np.arange(3); wdt=0.26
for j,(nm,col,idx) in enumerate([('Shares a clone group with India','#eb6834',1),('with Pakistan','#1baf7a',2),('with either','#52514e',3)]):
    for i,r in enumerate(rows):
        v=r[idx]; ax.bar(i+(j-1)*wdt,100*v/r[4],width=wdt-0.03,color=col,zorder=2); ax.text(i+(j-1)*wdt,100*v/r[4]+1.2,f'{v}',ha='center',fontsize=9,color=INK)
ax.set_xticks(x); ax.set_xticklabels([r[0] for r in rows],color=INK); ax.set_ylim(0,85); ax.grid(axis='y',color='#ecebe6',zorder=0)
for s_ in ['top','right','left']: ax.spines[s_].set_visible(False)
ax.spines['bottom'].set_color('#cfcdc6'); ax.tick_params(axis='y',length=0,labelcolor=INK2); ax.set_ylabel('Bangladeshi genomes (% of 60)',color=INK2)
ax.legend(handles=[Patch(color='#eb6834',label='Shares a clone group with India'),Patch(color='#1baf7a',label='Shares a clone group with Pakistan'),Patch(color='#52514e',label='Shares with either country')],frameon=False,fontsize=9,loc='upper right')
fig.suptitle('How many Bangladeshi genomes belong to clone groups seen in India or Pakistan?',x=0.01,ha='left',fontsize=12,fontweight='bold',color=INK,y=0.98)
fig.text(0.01,0.905,'Clone group = genomes joined by skani ANI at the stated threshold. Numbers above bars are genome counts. Lab strains removed.',fontsize=8.8,color=INK2)
fig.subplots_adjust(top=0.84,left=0.08,right=0.98,bottom=0.1); fig.savefig('fig14_sharing_clean.png',facecolor=BG)
print(rows)
