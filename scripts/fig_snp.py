import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
BG='#fcfcfb'; INK='#0b0b0b'; INK2='#52514e'; BD='#2a78d6'; IN='#eb6834'
rows=[ # label, nearest India/Pak, nearest Bangladeshi
 ("ST357 · rectal swab, 2023",32,8),
 ("ST357 · tracheal aspirate, 2023",32,8),
 ("ST357 · wastewater, 2022",48,38),
 ("ST5184 · blood (year not recorded)",28,None),
 ("ST308 · blood, 2024",46,None),
 ("ST773 · pus, 2022",31,17),
 ("ST not assigned (ST773 clone) · dental caries, 2019",28,17),
 ("ST274 · sesame plant, 2024",200,None),
 ("ST2556 · blood (year not recorded)",150,None)]
rows=sorted(rows,key=lambda r:r[1])
fig,ax=plt.subplots(figsize=(10.5,5.6),dpi=160); fig.patch.set_facecolor(BG); ax.set_facecolor(BG)
ax.axvspan(20,25,color='#e9e7e1',zorder=0)
ax.text(22.4,len(rows)-0.25,"20–25 SNPs: range reported for\nrecent transmission in hospital\nstudies (different methods)",ha='center',va='bottom',fontsize=8,color=INK2)
for i,(lab,i_d,b_d) in enumerate(rows):
    y=i
    if b_d is not None:
        ax.plot([b_d,i_d],[y,y],color='#cfcdc6',lw=1.6,zorder=1)
        ax.scatter(b_d,y,s=80,color=BD,zorder=3,edgecolor=BG,linewidth=1.5)
        ax.text(b_d/1.12,y,str(b_d),ha='right',va='center',fontsize=9,color=INK)
    ax.scatter(i_d,y,s=80,color=IN,zorder=3,edgecolor=BG,linewidth=1.5)
    ax.text(i_d*1.12,y,str(i_d),ha='left',va='center',fontsize=9,color=INK)
ax.set_yticks(range(len(rows))); ax.set_yticklabels([r[0] for r in rows],fontsize=9,color=INK)
ax.set_xscale('log'); ax.set_xlim(5,330); ax.set_xticks([5,10,20,50,100,200]); ax.set_xticklabels(['5','10','20','50','100','200'],fontsize=9,color=INK2)
ax.set_xlabel("SNP differences (core sites shared by the compared genomes)",fontsize=10,color=INK2)
ax.set_ylim(-0.7,len(rows)+1.0)
for s in ['top','right','left']: ax.spines[s].set_visible(False)
ax.spines['bottom'].set_color('#cfcdc6'); ax.tick_params(axis='y',length=0); ax.tick_params(axis='x',color='#cfcdc6')
ax.grid(axis='x',color='#ecebe6',lw=0.8,zorder=0); ax.set_axisbelow(True)
ax.legend(handles=[Line2D([],[],marker='o',ls='',color=IN,markersize=8,label='Nearest Indian or Pakistani genome'),
                   Line2D([],[],marker='o',ls='',color=BD,markersize=8,label='Nearest other Bangladeshi genome')],
          loc='lower right',frameon=False,fontsize=9,labelcolor=INK)
fig.suptitle("How close are Bangladeshi P. aeruginosa genomes to their nearest relatives?",x=0.01,ha='left',fontsize=12.5,color=INK,fontweight='bold',y=0.985)
fig.text(0.01,0.915,"Nine Bangladeshi genomes in clone groups shared with India or Pakistan (lab strains excluded). Predicted from public assemblies.",fontsize=9,color=INK2)
fig.subplots_adjust(left=0.36,right=0.97,top=0.86,bottom=0.12)
fig.savefig('fig12_snp_nearest_relatives.png',facecolor=BG)
