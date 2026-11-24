import pandas as pd, numpy as np
from scipy.stats import fisher_exact, binom

def multipletests(ps,method='holm'):
    p=np.array(ps); m=len(p)
    if method=='bonferroni': return None,np.minimum(1,p*m)
    o=np.argsort(p); adj=np.empty(m); run=0
    for i,idx in enumerate(o):
        run=max(run,(m-i)*p[idx]); adj[idx]=min(1,run)
    return None,adj
c=pd.read_csv('clusters_main_clean.csv',dtype={'st':str}); d=pd.read_csv('derep_main.csv',dtype={'st':str})
keep=set(c.asm_acc); dm=d[d.asm_acc.isin(keep)]
print('derep',dm.country.value_counts().to_dict())
def sets(c,dm):
    S={'All genomes':c,'One per country, BioProject, ST':dm}
    for col,l in [('cl_99.9','ANI>=99.9'),('cl_99.95','ANI>=99.95'),('cl_99.99','ANI>=99.99')]:
        S['One per country per clone group, '+l]=c.sort_values('n50',ascending=False).drop_duplicates([col,'country'])
    return S
def run(flag,label):
    cc=c.copy(); cc['carb']=cc.asm_acc.map(flag) if callable(flag) else cc.carb
    carb=dict(zip(cc.asm_acc,cc.carb)); dd=dm.copy(); dd['carb']=dd.asm_acc.map(carb)
    out=[]; ps=[]
    for name,df in sets(cc,dd).items():
        row=[name]
        b=df[df.country=='Bangladesh']; kb=int(b.carb.sum()); nb=len(b)
        for o in ['India','Pakistan']:
            s=df[df.country==o]; k=int(s.carb.sum()); n=len(s)
            p=fisher_exact([[kb,nb-kb],[k,n-k]])[1]; ps.append(p); row+= [f'{k}/{n}',round(p,4)]
        out.append([row[0],f'{kb}/{nb}']+row[1:])
    print('\n##',label); 
    for r in out: print(r)
    return ps
p0=run(None,'as reported')
# strict: partial hits <90% coverage -> negative
x=['GCA_000506885.1']
flag=lambda a: False if a in x else bool(c.set_index('asm_acc').carb[a])
p1=run(flag,'strict: 2 partial hits (<90% cov) = negative')
for nm,ps in [('reported',p0),('strict',p1)]:
    h=multipletests(ps,method='holm')[1]; b=multipletests(ps,method='bonferroni')[1]
    print(nm,'raw',[round(p,4) for p in ps]); print(' holm',[round(v,3) for v in h]); print(' min holm',round(h.min(),3))
# power
def power(nb,pb,n2,p2,alpha=0.05):
    pw=0
    from scipy.stats import fisher_exact
    kb=np.arange(nb+1); k2=np.arange(n2+1)
    wb=binom.pmf(kb,nb,pb); w2=binom.pmf(k2,n2,p2)
    # exact fisher power; heavy -> use normal approx two-proportion z with continuity as approximation
    return None
from scipy.stats import norm
def pw_z(nb,n2,p1,p2):
    pbar=(nb*p1+n2*p2)/(nb+n2); se0=np.sqrt(pbar*(1-pbar)*(1/nb+1/n2)); se1=np.sqrt(p1*(1-p1)/nb+p2*(1-p2)/n2)
    return norm.cdf((abs(p1-p2)-1.96*se0)/se1)
for lab,nb,n2 in [('BD 48 vs IN 306',48,306),('BD 48 vs PK 84',48,84)]:
    for base in [0.10,0.18,0.25]:
        # smallest absolute increase giving >=80% power (BD higher or lower)
        for sign in [+1,-1]:
            for dlt in np.arange(0.01,0.5,0.005):
                p1=base+sign*dlt
                if 0<p1<1 and pw_z(nb,n2,p1,base)>=0.8: print(lab,'comparator',base,'BD',('higher' if sign>0 else 'lower'),'by',round(dlt*100),'pp (BD',round(p1*100),'%)'); break
