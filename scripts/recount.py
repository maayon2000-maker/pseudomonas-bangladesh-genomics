import pandas as pd, numpy as np
from scipy.stats import fisher_exact
from math import sqrt
d=pd.read_csv('clusters_main.csv')
excl=['GCA_051063055.1','GCA_051063075.1','GCA_058267545.1','GCA_058267575.1','GCA_058890545.1','GCA_058931605.1','GCA_059878655.1','GCA_059878685.1']
print('main set',len(d),d.country.value_counts().to_dict(), 'excluded present:',d.asm_acc.isin(excl).sum())
def wilson(k,n,z=1.96):
    if n==0: return (np.nan,np.nan)
    p=k/n; den=1+z*z/n; c=(p+z*z/(2*n))/den; h=z*sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return (100*(c-h),100*(c+h))
def block(df,tag):
    print('\n=====',tag,'n=',len(df),df.country.value_counts().to_dict())
    for c in ['Bangladesh','India','Pakistan']:
        s=df[df.country==c]; k=int(s.carb.sum()); lo,hi=wilson(k,len(s))
        print(f'crude carbapenemase {c}: {k}/{len(s)} = {100*k/len(s):.1f}% (Wilson {lo:.1f}-{hi:.1f})')
    for col in ['cl_99.9','cl_99.95','cl_99.99']:
        bd=df[df.country=='Bangladesh']
        # cluster membership by country
        g=df.groupby(col).country.agg(set)
        def sh(row,o): return o in g[row[col]]
        bd_any=sum(any(x!='Bangladesh' for x in g[r[col]]) for _,r in bd.iterrows())
        bd_in=sum('India' in g[r[col]] for _,r in bd.iterrows())
        bd_pk=sum('Pakistan' in g[r[col]] for _,r in bd.iterrows())
        print(f'BD genomes sharing a clone group ({col}): any other country {bd_any}/{len(bd)}, India {bd_in}, Pakistan {bd_pk}')
    # dereplicate: one per country per cluster 99.99 (best N50)
    for col in ['cl_99.99']:
        dr=df.sort_values('n50',ascending=False).drop_duplicates([col,'country'])
        print(f'-- dereplicated one/country/cluster ({col}): ',dr.country.value_counts().to_dict())
        for c in ['Bangladesh','India','Pakistan']:
            s=dr[dr.country==c]; k=int(s.carb.sum()); lo,hi=wilson(k,len(s))
            print(f'   {c}: {k}/{len(s)} = {100*k/len(s):.1f}% (Wilson {lo:.1f}-{hi:.1f})')
        b=dr[dr.country=='Bangladesh']
        for o in ['India','Pakistan']:
            s=dr[dr.country==o]
            t=[[int(b.carb.sum()),len(b)-int(b.carb.sum())],[int(s.carb.sum()),len(s)-int(s.carb.sum())]]
            OR,p=fisher_exact(t); print(f'   Fisher BD vs {o}: p={p:.3f}, OR={OR:.2f}')
    # ST-level dereplication
    dr=df.sort_values('n50',ascending=False).drop_duplicates(['st','country'])
    print('-- ST-level (one per ST per country):',dr.country.value_counts().to_dict())
    b=dr[dr.country=='Bangladesh']
    for c in ['Bangladesh','India','Pakistan']:
        s=dr[dr.country==c]; k=int(s.carb.sum()); lo,hi=wilson(k,len(s)); print(f'   {c}: {k}/{len(s)} = {100*k/len(s):.1f}% ({lo:.1f}-{hi:.1f})')
    for o in ['India','Pakistan']:
        s=dr[dr.country==o]
        t=[[int(b.carb.sum()),len(b)-int(b.carb.sum())],[int(s.carb.sum()),len(s)-int(s.carb.sum())]]
        OR,p=fisher_exact(t); print(f'   Fisher BD vs {o}: p={p:.3f}')
block(d,'BEFORE (as in earlier analysis)')
d2=d[~d.asm_acc.isin(excl)].copy()
block(d2,'AFTER removing 8 genomes')
d2.to_csv('clusters_main_clean.csv',index=False)
# BD carb list
print(d2[(d2.country=='Bangladesh')&d2.carb][['asm_acc','st','cat','cl_99.99']])
