from pathlib import Path
from itertools import product
import ast,json
import numpy as np,pandas as pd
from scipy import stats
OUT=(Path(__import__('os').environ['CRRM_WORKSPACE']) / 'Result5_analysis');R4=(OUT.parent/'Result4') if (OUT.parent/'Result4').exists() else (OUT.parent/'Result4_analysis');plan=json.loads((OUT/'parameters/analysis_plan_v1.json').read_text())
code=(OUT/'scripts/02_score_and_compare.py').read_text();tree=ast.parse(code)
for name in ['infer','bh','hc3']:
    n=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name==name);exec(ast.get_source_segment(code,n))
core=pd.read_csv(OUT/'parameters/frozen_CRRM388.csv');sets=json.loads((OUT/'parameters/score_gene_sets.json').read_text());sign=np.where(core.frozen_direction.eq('up'),1.,-1.);targets=['CRRM388','IFN-list overlap','Outside IFN lists','Neutrophil degranulation'];scores=pd.read_csv(OUT/'results/profile_scores.csv');old=pd.read_csv(R4/'results/rank_score_profiles.csv');ranks=[]
for source in ['GSE38900_GPL6884','GSE42026']:
    x=pd.read_csv(OUT/'data'/f'{source}_log_expression.csv.gz',index_col=0);g=scores[scores.source_label.eq(source)].copy();assert list(x.columns)==g.geo_accession.tolist();rr=(stats.rankdata(x.to_numpy().T,axis=1)-1)/(len(x)-1)-.5;idx={id:i for i,id in enumerate(x.index)}
    for t in targets:
        members=[id for id in core.human_ensembl if id in sets[t] and id in idx];sg=core.set_index('human_ensembl').loc[members].frozen_direction.map({'up':1.,'down':-1.}).to_numpy();g[t]=(rr[:,[idx[id] for id in members]]*sg).mean(1)
    ranks.append(g)
for source in ['Infant_PBMC','CD1c_DC']:
    g=scores[scores.source_label.eq(source)].copy();o=old[old.source_label.eq(source)].set_index('profile_id').loc[g.profile_id]
    for t in targets:g[t]=o[t].to_numpy()
    ranks.append(g)
rank=pd.concat(ranks,ignore_index=True);rank.drop(columns=[c for c in rank if c.startswith('contact_')]).to_csv(OUT/'results/rank_profile_scores.csv',index=False)
rows=[]
for s in plan['specifications']:
    if s['total_family']=='Direct virus total':continue
    g=rank[rank.source_label.eq(s['source'])&rank.analysis_celltype.eq(s['cell'])&rank.view.eq('canonical')]
    if s['source'] in ['Infant_PBMC','CD1c_DC']:g=g[g.cells.ge(30)]
    a=g[g.analysis_group.eq(s['case'])];b=g[g.analysis_group.eq(s['control'])]
    if s['paired']:
        aa=a[a.time.isin(['4h','8h','18h'])].groupby(['analysis_unit','time'])[targets].mean();bb=b[b.time.isin(['4h','8h','18h'])].groupby(['analysis_unit','time'])[targets].mean();common=aa.index.intersection(bb.index);aa=aa.loc[common];bb=bb.loc[common];a=aa.groupby(level=0).mean();b=bb.groupby(level=0).mean()
    for t in targets:rows.append(dict(key=s['key'],family=s['total_family'],score=t,method='Signed transcriptome rank',n_case=len(a),n_control=len(b),**infer(a[t],b[t],s['paired'])))
rs=pd.DataFrame(rows)
for (family,target),ix in rs.groupby(['family','score']).groups.items():rs.loc[ix,'adjusted_P']=bh(rs.loc[ix,'P'])
rs.to_csv(OUT/'results/rank_induction_and_replication_effects.csv',index=False)
comp=[];ref=pd.read_csv(OUT/'results/PBMC_fixed_healthy_reference.csv',index_col=0).iloc[:,0];keep=set(pd.read_csv(OUT/'results/PBMC_composition_values.csv').analysis_unit)
four=rank[rank.source_label.eq('Infant_PBMC')&rank.analysis_celltype.isin(['B_cell','T_cell','NK_cell','Monocyte'])&rank.analysis_unit.isin(keep)&rank.view.eq('canonical')].copy();four['fraction']=four.cells/four.groupby('analysis_unit').cells.transform('sum')
for t in ['CRRM388','IFN-list overlap','Outside IFN lists']:
    m=four.pivot(index='analysis_unit',columns='analysis_celltype',values=t);fractions=four.pivot(index='analysis_unit',columns='analysis_celltype',values='fraction').reindex(columns=m.columns);groups=four.drop_duplicates('analysis_unit').set_index('analysis_unit').loc[m.index,'analysis_group'];w=ref.reindex(m.columns).to_numpy();ia=np.flatnonzero(groups.eq('RSV'));ib=np.flatnonzero(groups.eq('COVID'));ih=np.flatnonzero(groups.eq('Healthy'))
    for method,y in [('Observed composition',(m*fractions).sum(1)),('Fixed healthy composition',m@w)]:
        random=np.random.default_rng(20261113);boot=[]
        for _ in range(2000):
            a=random.choice(ia,len(ia),replace=True);b=random.choice(ib,len(ib),replace=True);h=random.choice(ih,len(ih),replace=True);ww=fractions.iloc[h].mean().to_numpy();aa=(m.iloc[a].to_numpy()*(fractions.iloc[a].to_numpy() if method=='Observed composition' else ww)).sum(1).mean();bb=(m.iloc[b].to_numpy()*(fractions.iloc[b].to_numpy() if method=='Observed composition' else ww)).sum(1).mean();boot.append(aa-bb)
        comp.append(dict(score=t,method=method,scoring='Signed transcriptome rank',estimate=y[groups.eq('RSV')].mean()-y[groups.eq('COVID')].mean(),lower=np.quantile(boot,.025),upper=np.quantile(boot,.975),n_RSV=len(ia),n_COVID=len(ib),n_healthy_reference=len(ih)))
pd.DataFrame(comp).to_csv(OUT/'results/PBMC_rank_composition_effects.csv',index=False)

# Gene HC3 sensitivity uses ALL388 hypotheses in the family, regardless of the5replicated findings.
source='GSE38900_GPL6884';x=pd.read_csv(OUT/'data'/f'{source}_log_expression.csv.gz',index_col=0);m=pd.read_csv(OUT/'data'/f'{source}_sample_metadata.csv');r=pd.read_csv(OUT/'data/gene_reference_parameters.csv');r=r[r.source.eq(source)].set_index('human_ensembl');sel=m.condition.isin(['RSV','Influenza']);gene=[]
for id in core.human_ensembl:
    if id not in x.index or not r.loc[id,'variable']:continue
    sg=1 if core.set_index('human_ensembl').loc[id,'frozen_direction']=='up' else -1;y=(x.loc[id].to_numpy()-r.loc[id,'mean'])/r.loc[id,'SD']*sg
    constant=(np.ptp(y[m.condition.eq('RSV')])<=1e-12 and np.ptp(y[m.condition.eq('Influenza')])<=1e-12)
    rr=hc3(y[sel],m.loc[sel,'condition'].eq('RSV').to_numpy(),m.loc[sel,'age_months'].to_numpy(),m.loc[sel,'sex'].to_numpy());
    if constant:rr.update(lower=np.nan,upper=np.nan,P=np.nan)
    gene.append(dict(human_ensembl=id,**rr))
gene=pd.DataFrame(gene);gene['adjusted_P_BH388']=bh(gene.P,388);gene.to_csv(OUT/'results/Dallas_gene_age_sex_HC3.csv',index=False)
rep=pd.read_csv(OUT/'results/locked_gene_replication.csv');rep=rep[rep.replicated].merge(gene,on='human_ensembl');rep.to_csv(OUT/'results/replicated_genes_age_sex_sensitivity.csv',index=False)
print(rs[rs.score.eq('CRRM388')][['key','estimate','adjusted_P']].to_string(index=False))
print(rep[['display_symbol','estimate','adjusted_P_BH388']].to_string(index=False))
