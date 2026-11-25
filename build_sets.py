import sys, os, re; sys.path.insert(0,'.')
from hpo import *
import pandas as pd, json
# Corrections (Nov 12, 2026). Both off reproduces the original framework exactly.
#   DROP_ZERO=1  annotations recorded in 0 of n individuals (0/n, 0%) are not positive matches
#   SOD_SEP=1    Septo-optic dysplasia (HP:0100842) is only its own inclusion flag and no longer
#                counts toward the optic and midline features
DROP_ZERO=os.environ.get('DROP_ZERO','0')=='1'
SOD_SEP=os.environ.get('SOD_SEP','0')=='1'
ONLY=os.environ.get('ONLY')            # run a single term set, e.g. ONLY=S1c_curated
TAG=os.environ.get('TAG','')
def zero_freq(f):
    f=str(f).strip(); m=re.match(r'^(\d+)/(\d+)$',f)
    if m: return int(m.group(1))==0
    m=re.match(r'^([\d.]+)%$',f)
    return bool(m) and float(m.group(1))==0
T,C,ALT=load_obo('data/hp.obo')
g=pd.read_pickle('data/g2p.pkl'); g['hpo_id']=g.hpo_id.map(lambda x: ALT.get(x,x))
g['zero_freq']=g.frequency.map(zero_freq)
print(f"zero-frequency annotation rows: {int(g.zero_freq.sum())} (DROP_ZERO={DROP_ZERO}, SOD_SEP={SOD_SEP})")
if DROP_ZERO: g=g[~g.zero_freq]
D=lambda *r: set().union(*[desc(x,C) for x in r])
MID=D('HP:0001273','HP:0007375')-{'HP:0011750'}  # drop cavum septum pellucidum if present
cav=[t for t in MID if 'avum' in T[t]['name']]; MID-=set(cav)
OPT_S=D('HP:0008058')|{'HP:0000587'}; OPT_B=D('HP:0000587','HP:0000648')
PIT_C=D('HP:0012503','HP:0034977','HP:0011734','HP:0034323','HP:0000863')
PIT_HH=D('HP:0000044')
EXACT={'midline':{"HP:0001331","HP:0007375","HP:0001274","HP:0007370","HP:0033725"},
 'optic':{"HP:0000609","HP:0008058","HP:0000587"},
 'pituitary':{"HP:0000871","HP:0000830","HP:0040075","HP:0011747","HP:0012503","HP:0000864","HP:0011751"}}
from curated import MID as CM, OPT as CO, PIT as CP
if SOD_SEP: CM=CM-{'HP:0100842'}; CO=CO-{'HP:0100842'}
# Term sets. The published analysis uses S1c_curated only (run_all.sh sets ONLY=S1c_curated).
# The other sets were exploratory comparisons made while curating the terms.
S={
 'S1c_curated':{'midline':CM,'optic':CO,'pituitary':CP},
 'S1c_curated_bare864':{'midline':CM,'optic':CO,'pituitary':CP|{'HP:0000864'}},
 'S0_exact_original':EXACT,
 'S1_expanded_core':{'midline':MID,'optic':OPT_S,'pituitary':PIT_C},
 'S2_S1+broad_optic':{'midline':MID,'optic':OPT_B|OPT_S,'pituitary':PIT_C},
 'S3_S1+HH':{'midline':MID,'optic':OPT_S,'pituitary':PIT_C|PIT_HH},
 'S4_S1+bare864':{'midline':MID,'optic':OPT_S,'pituitary':PIT_C|{'HP:0000864'}},
}
def run(dom, same_disease=False, add_sod_term=True):
    rows=[]
    for gid,sub in g.groupby('ncbi_gene_id'):
        hit={d:set(sub.disease_id[sub.hpo_id.isin(v)]) for d,v in dom.items()}
        present=[d for d,v in hit.items() if v]
        # same-disease co-occurrence: any disease with >=2 domains
        dz={}
        for d,v in hit.items():
            for x in v: dz.setdefault(x,set()).add(d)
        same=max([len(v) for v in dz.values()]+[0])
        sodt=(sub.hpo_id=='HP:0100842').any()
        nm=dict(zip(sub.hpo_id,sub.hpo_name))
        rec=dict(ncbi_gene_id=gid,gene_symbol=sub.gene_symbol.iloc[0],
            midline=bool(hit['midline']),optic=bool(hit['optic']),pituitary=bool(hit['pituitary']),
            n_domains=len(present),max_domains_same_disease=same,SOD_term=sodt)
        # traceability: matched terms and supporting disease IDs per feature
        for d,v in dom.items():
            m_=sub[sub.hpo_id.isin(v)]
            rec[d+'_terms']='; '.join(sorted({f'{nm[t]} ({t})' for t in m_.hpo_id}))
            rec[d+'_diseases']='; '.join(sorted(v_ for v_ in hit[d]))
        rec['SOD_term_diseases']='; '.join(sorted(set(sub.disease_id[sub.hpo_id=='HP:0100842'])))
        rec['same_disease_ids']='; '.join(f"{x} ({' + '.join(sorted(v))})" for x,v in sorted(dz.items()) if len(v)>=2)
        rows.append(rec)
    df=pd.DataFrame(rows)
    df['selected']=(df.n_domains>=2)|(add_sod_term & df.SOD_term)
    df['selected_same_disease']=(df.max_domains_same_disease>=2)|(add_sod_term & df.SOD_term)
    return df
rep=pd.read_csv('literature_genes.csv',dtype=str)
R=set(rep[rep.in_original_27=='yes'].gene)   # the original 27 published-report genes
# diagnostic and candidate genes in the cohort (Individual numbers as in the paper); printed for reference only
coh={'NFIA':1,'TUBB2B':2,'DOCK4':3,'HESX1':3,'LSS':4,'COL4A1':5,'SPTB':6,'ASXL1':6,'USP5':6,'SLCO1A2':7,'ARMCX5':8}
out={}
for name,dom in S.items():
    if ONLY and name!=ONLY: continue
    df=run(dom, add_sod_term=(name!='S0_exact_original'))
    df.to_csv(f'outputs/{name}{TAG}_genes.csv',index=False)
    sel=df[df.selected]; sel3=sel[sel.n_domains==3]
    ss=df[df.selected_same_disease]
    n_same_id=int((sel.max_domains_same_disease>=2).sum())   # >=2 features under one disease ID, no SOD-term pass
    print(f"\n{name}: selected={len(sel)} all3={len(sel3)} same-disease (incl. SOD-term pass)={len(ss)} "
          f">=2 features under one disease ID={n_same_id} | reported recovered={len(R&set(sel.gene_symbol))}/27 (same-dis {len(R&set(ss.gene_symbol))})")
    print("  missed reported:", sorted(R-set(sel.gene_symbol)))
    print("  cohort:", {k:(int(df[df.gene_symbol==k].n_domains.iloc[0]) if (df.gene_symbol==k).any() else 'noHPO') for k in coh})
    out[name]=df
