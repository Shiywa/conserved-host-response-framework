"""Rebuild only GSE246622; retain other studies' raw results and original families."""
from pathlib import Path
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import ast,hashlib,itertools,json,math
import numpy as np,pandas as pd
from scipy import stats

BASE=(Path(__import__('os').environ['CRRM_WORKSPACE']) / 'clinical_otherwisehealthy_262_r1');OUT=BASE/'main'
REV=Path(str(Path(__import__('os').environ['CRRM_WORKSPACE']) / ''))
OLD=REV/'phase20_Result3_clinical_analysis_20261112'
source=(BASE/'audit/original_scripts/main/04_analyze_clinical_CRRM388.py').read_text()
tree=ast.parse(source)
for node in tree.body:
 if isinstance(node,ast.FunctionDef) and node.name in ['truth','padjust','welch','paired','ols_hc3','independent_permutation','signflip']:
  exec(compile(ast.Module(body=[node],type_ignores=[]),'<original-statistical-functions>','exec'))
rng=np.random.default_rng(20261120)
def write(df,name):df.to_csv(OUT/'tables'/name,index=False)
def read(name):return pd.read_csv(OLD/'tables'/name,float_precision='round_trip')
def replace(df,name,key='dataset',labels=None):
 old=read(name);labels=labels or df[key].unique().tolist()
 return pd.concat([df,old[~old[key].isin(labels)]],ignore_index=True,sort=False)
def dump(p,x):Path(p).write_text(json.dumps(x,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
core=pd.read_csv(OUT/'parameters/frozen_CRRM388.csv').set_index('human_ensembl')
br=pd.read_csv(OUT/'parameters/frozen_functional_branches.csv');branches=br.branch.drop_duplicates().tolist()
sets=json.loads((REV/'Result6_analysis/parameters/score_gene_sets.json').read_text())
data={};sample_rows=[];fcov=[];covrows=[];normrows=[];scov=[]
for ds in ['GSE246622','GSE246622_RMA']:
 x=pd.read_csv(OUT/'data/processed'/f'{ds}_log_expression.csv.gz',index_col=0,float_precision='round_trip')
 m=pd.read_csv(OUT/'data/processed'/f'{ds}_metadata.csv');assert x.columns.tolist()==m.sample_id.tolist() and len(m)==262
 mu=x.mean(axis=1);sd=x.std(axis=1,ddof=1);valid=sd.gt(1e-8)&np.isfinite(sd)
 for g,r in core.iterrows():covrows.append(dict(dataset=ds,human_ensembl=g,display_symbol=r.display_symbol,frozen_direction=r.frozen_direction,present=g in x.index,scored=bool(valid.get(g,False)),reference_sd=sd.get(g,np.nan)))
 x=x.loc[valid];z=x.sub(mu.loc[x.index],axis=0).div(sd.loc[x.index],axis=0)
 sg=core.loc[x.index,'frozen_direction'].map({'up':1,'down':-1});dz=z.mul(sg,axis=0);up=sg.eq(1)
 assert len(x)/388>=.7 and up.sum()/280>=.5 and (~up).sum()/108>=.5
 scores=pd.DataFrame({'Total':dz.mean(axis=0),'Up contribution':dz.loc[up].sum(axis=0)/len(x),'Down contribution':dz.loc[~up].sum(axis=0)/len(x)})
 assert np.max(np.abs(scores.Total-scores['Up contribution']-scores['Down contribution']))<1e-12
 for b in branches:
  ids=br.loc[br.branch.eq(b),'human_ensembl'].tolist();hit=[g for g in ids if g in dz.index]
  eligible=len(hit)>=5 and len(hit)/len(ids)>=.5;scores[b]=dz.loc[hit].mean(0) if eligible else np.nan
  fcov.append(dict(dataset=ds,branch=b,n_frozen=len(ids),n_measured=len(hit),fraction_measured=len(hit)/len(ids),eligible=eligible))
 # These inherited IFN partitions are source measurements, without new hypothesis tests.
 for name in ['IFN-list overlap','Outside IFN lists']:
  ids=sets[name];hit=[g for g in ids if g in dz.index];scores[name]=dz.loc[hit].mean(0) if len(hit)>=5 and len(hit)/len(ids)>=.5 else np.nan
  fcov.append(dict(dataset=ds,branch=name,n_frozen=len(ids),n_measured=len(hit),fraction_measured=len(hit)/len(ids),eligible=len(hit)>=5 and len(hit)/len(ids)>=.5))
 scores=scores.reset_index(drop=True)
 sample_rows.append(pd.concat([m.reset_index(drop=True),scores.add_prefix('score_')],axis=1))
 z.to_csv(OUT/'data/processed'/f'{ds}_gene_z.csv.gz',index_label='human_ensembl',compression='gzip')
 pd.DataFrame({'human_ensembl':x.index,'reference_mean':mu.loc[x.index],'reference_sample_SD':sd.loc[x.index],'frozen_sign':sg}).to_csv(OUT/'parameters'/f'{ds}_gene_reference.csv',index=False)
 normrows.append(dict(dataset=ds,study_family='GSE246622',n_profiles=len(m),n_reference=len(m),n_mapped=len(valid),n_scored=len(x),method='Author batch-corrected log expression;median unique probes' if ds=='GSE246622' else 'Author RMA-only log expression;median unique probes;same population'))
 scov.append(dict(dataset=ds,n_present=len(valid),n_scored=len(x),n_up=int(up.sum()),n_down=int((~up).sum()),eligible=True))
 data[ds]=dict(x=x,z=z,m=m,scores=scores)
write(replace(pd.concat(sample_rows,ignore_index=True),'07_sample_scores.csv'),'07_sample_scores.csv')
for name,new in [('02_member_coverage.csv',covrows),('05_normalization_summary.csv',normrows),('06_coverage_summary.csv',scov),('08_functional_branch_coverage.csv',fcov)]:write(replace(pd.DataFrame(new),name),name)
# Retain only this revision's eligibility ledger for the whole-blood source.
old_manifest=read('01_clinical_source_manifest.csv');ledger=pd.read_csv(BASE/'source_data/GSE246622_population_eligibility.tsv',sep='\t')
wb=ledger.rename(columns={'laboratory_id':'sample_id','clinical_group':'group'})
wb['study_family']='GSE246622';wb['source_visit']=wb.source_visit_type
wb['visit']=np.where(wb.source_visit_type.eq('Healthy'),'Healthy',np.where(wb.source_visit_type.eq('RSV visit'),'Acute','Convalescent'))
wb['unit_id']='GSE246622:'+wb.sample_id;wb['tissue']='Whole blood'
write(pd.concat([wb,old_manifest[old_manifest.study_family.ne('GSE246622')]],ignore_index=True,sort=False),'01_clinical_source_manifest.csv')

o=data['GSE246622'];m=o['m'];keep=m.group.isin(['Severe','Mild']);ix=np.flatnonzero(keep);g=m.loc[keep,'group'].eq('Severe').to_numpy();Y=o['scores'].iloc[ix].reset_index(drop=True)
c=dict(contrast_key='GSE246622_Severe_Mild',dataset='GSE246622',positive='Severe',reference='Mild',role='Primary',n_positive=int(g.sum()),n_reference=int((~g).sum()),n_measured=len(o['x']))
p,nperm,method=independent_permutation(Y.Total,g)
effects=replace(pd.DataFrame([dict(**c,**welch(Y.Total,g),resampling_P=p,n_resamples=nperm,resampling_method=method)]),'09_total_clinical_effects.csv')
ii=effects.role.eq('Primary');effects.loc[ii,'primary_Holm_P']=padjust(effects.loc[ii,'p_value'],'holm',2);effects.loc[ii,'resampling_Holm_P']=padjust(effects.loc[ii,'resampling_P'],'holm',2);write(effects,'09_total_clinical_effects.csv')
fun=[]
for b in branches:
 hit=next(r for r in fcov if r['dataset']=='GSE246622' and r['branch']==b)
 fun.append(dict(**{k:v for k,v in c.items() if k!='n_measured'},n_measured_total=c['n_measured'],branch=b,**welch(Y[b],g),n_frozen=hit['n_frozen'],n_measured_branch=hit['n_measured'],adjustment_family='BH14 strict severity branches'))
fun=replace(pd.DataFrame(fun),'10_functional_clinical_effects.csv')
for role,n in [('Primary',14),('Hospitalization sensitivity',7)]:
 sel=fun.role.eq(role);fun.loc[sel,'BH_adjusted_P']=padjust(fun.loc[sel,'p_value'],'BH',n)
write(fun,'10_functional_clinical_effects.csv')
comp=replace(pd.DataFrame([dict(**c,component=b,**welch(Y[b],g)) for b in ['Up contribution','Down contribution']]),'11_directional_contributions.csv')
for role,n in [('Primary',4),('Hospitalization sensitivity',2)]:
 sel=comp.role.eq(role);comp.loc[sel,'BH_adjusted_P']=padjust(comp.loc[sel,'p_value'],'BH',n)
assert abs(comp.loc[comp.dataset.eq('GSE246622'),'effect'].sum()-effects.loc[effects.dataset.eq('GSE246622'),'effect'].iloc[0])<1e-12
write(comp,'11_directional_contributions.csv')
ee=o['x'].iloc[:,ix].to_numpy();eff=ee[:,g].mean(1)-ee[:,~g].mean(1)
gene=pd.DataFrame({'contrast_key':c['contrast_key'],'dataset':'GSE246622','human_ensembl':o['x'].index,'display_symbol':core.loc[o['x'].index,'display_symbol'].to_numpy(),'frozen_direction':core.loc[o['x'].index,'frozen_direction'].to_numpy(),'log_expression_effect':eff})
write(replace(gene,'12_clinical_gene_effects.csv'),'12_clinical_gene_effects.csv')
infl=[]
for i,unit in enumerate(m.unit_id.iloc[ix]):
 kk=np.arange(len(Y))!=i;infl.append(dict(contrast_key=c['contrast_key'],removed_unit=unit,full_effect=effects.loc[effects.dataset.eq('GSE246622'),'effect'].iloc[0],**welch(Y.Total.to_numpy()[kk],g[kk])))
write(replace(pd.DataFrame(infl),'13_patient_deletion_influence.csv',key='contrast_key'),'13_patient_deletion_influence.csv')
groups=[]
for ds,ob in data.items():
 for (visit,group),v in ob['m'].groupby(['visit','group'],sort=False):
  y=ob['scores'].Total.iloc[v.index];se=y.std(ddof=1)/np.sqrt(len(v));q=stats.t.ppf(.975,len(v)-1)
  groups.append(dict(dataset=ds,visit=visit,group=group,n=len(v),mean=y.mean(),SD=y.std(ddof=1),median=y.median(),CI95_low=y.mean()-q*se,CI95_high=y.mean()+q*se))
write(replace(pd.DataFrame(groups),'14_group_score_summary.csv'),'14_group_score_summary.csv')
acute=m.visit.eq('Acute');rr=stats.spearmanr(m.loc[acute,'group'].map({'Mild':1,'Moderate':2,'Severe':3}),o['scores'].Total[acute]);kw=stats.kruskal(*(o['scores'].Total[m.group.eq(v)] for v in ['Mild','Moderate','Severe']))
write(replace(pd.DataFrame([dict(dataset='GSE246622',n=int(acute.sum()),Spearman_rho=rr.statistic,Spearman_P=rr.pvalue,Kruskal_Wallis=kw.statistic,Kruskal_P=kw.pvalue)]),'15_RSV_only_grade_patterns.csv'),'15_RSV_only_grade_patterns.csv')
grade=[]
for pos,ref in [('Moderate','Mild'),('Severe','Mild'),('Severe','Moderate')]:
 mask=m.group.isin([pos,ref]);gg=m.loc[mask,'group'].eq(pos).to_numpy()
 grade.append(dict(contrast_key=f'GSE246622_{pos}_{ref}',dataset='GSE246622',positive=pos,reference=ref,role='Secondary grade comparison',n_positive=int(gg.sum()),n_reference=int((~gg).sum()),n_measured=len(o['x']),**welch(o['scores'].Total[mask],gg)))
grade=pd.DataFrame(grade);grade['BH3_within_cohort']=padjust(grade.p_value,'BH',3);write(replace(grade,'16_secondary_grade_comparisons.csv'),'16_secondary_grade_comparisons.csv')
adj=[];adjfun=[]
for model in ['unadjusted','adjusted']:
 d=pd.read_csv(OUT/'models'/f'GSE246622_{model}_design.csv');X=d.iloc[:,1:].to_numpy();cc=np.zeros(X.shape[1]);cc[d.columns[1:].tolist().index('high')]=1
 lookup=dict(zip(m.sample_id,m.index));loc=[lookup[s] for s in d.sample_id];V=o['scores'].iloc[loc]
 adj.append(dict(dataset='GSE246622',model=model,n=len(d),covariates='None' if model=='unadjusted' else 'Age spline (3 df), sex',SE_method='HC3',**ols_hc3(V.Total,X,cc)))
 if model=='adjusted':
  for b in branches:
   cv=next(r for r in fcov if r['dataset']=='GSE246622' and r['branch']==b)
   adjfun.append(dict(dataset='GSE246622',branch=b,n=len(d),n_measured_branch=cv['n_measured'],n_frozen=cv['n_frozen'],covariates='Age spline3df, sex',SE_method='HC3',role='Additional exploratory covariate sensitivity',**ols_hc3(V[b],X,cc)))
adj=replace(pd.DataFrame(adj),'17_covariate_adjusted_effects.csv');sel=adj.model.eq('adjusted');adj.loc[sel,'adjusted_model_BH2_P']=padjust(adj.loc[sel,'p_value'],'BH',2);write(adj,'17_covariate_adjusted_effects.csv')
af=replace(pd.DataFrame(adjfun),'30_covariate_adjusted_branch_sensitivity.csv');af['BH14_adjusted_branch_sensitivity_P']=padjust(af.p_value,'BH',14);write(af,'30_covariate_adjusted_branch_sensitivity.csv')

# Common-measurement checks: establish the intersection anew; other study parent scales stay fixed.
other_genes={ds:set(pd.read_csv(OUT/'data/processed'/f'{ds}_gene_z.csv.gz',index_col=0).index) for ds in ['GSE69606','GSE105450']}
common_genes=set(o['z'].index).intersection(*other_genes.values())
old_common=set(read('22_common_measured_members.csv').human_ensembl)
assert common_genes==old_common, 'If coverage changed, dependent other-study common-member scores also require updating'
commonrows=[]
for ds in ['GSE246622']:
 hit=[v for v in data[ds]['z'].index if v in common_genes];sg=core.loc[hit,'frozen_direction'].map({'up':1,'down':-1});v=data[ds]['z'].loc[hit].mul(sg,axis=0).mean(0)
 commonrows.append(dict(**c,n_common=len(hit),n_common_up=int(sg.eq(1).sum()),n_common_down=int(sg.eq(-1).sum()),full_score_correlation=np.corrcoef(data[ds]['scores'].Total,v)[0,1],**welch(v.iloc[ix],g)))
write(replace(pd.DataFrame(commonrows),'21_common_member_sensitivity.csv'),'21_common_member_sensitivity.csv')
write(core.loc[sorted(common_genes)].reset_index()[['human_ensembl','display_symbol','frozen_direction']],'22_common_measured_members.csv')
ob=data['GSE246622_RMA'];mask=ob['m'].group.isin(['Severe','Mild']);gg=ob['m'].loc[mask,'group'].eq('Severe').to_numpy()
rma=dict(**{**c,'dataset':'GSE246622_RMA','contrast_key':'GSE246622_RMA_Severe_Mild','role':'RMA-only sensitivity','n_measured':len(ob['x'])},primary_score_correlation=np.corrcoef(o['scores'].Total,ob['scores'].Total)[0,1],**welch(ob['scores'].Total[mask],gg));write(pd.DataFrame([rma]),'23_RMA_source_sensitivity.csv')
covdist=[]
for group,v in m.groupby('group',sort=False):
 age=v.age_days;covdist.append(dict(dataset='GSE246622',group=group,n=len(v),age_unit='days',age_median=age.median(),age_Q1=age.quantile(.25),age_Q3=age.quantile(.75),n_female=int(v.sex.eq('Female').sum()),n_age_missing=int(age.isna().sum()),n_sex_missing=int(v.sex.isna().sum())))
write(replace(pd.DataFrame(covdist),'25_clinical_covariate_summary.csv'),'25_clinical_covariate_summary.csv')
strict=m.iloc[ix];agelevel=np.where(strict.age_days.lt(90),'Younger','Older');sub=[]
for variable,level,kk in [('Age',v,agelevel==v) for v in ['Younger','Older']]+[('Sex',v,strict.sex.eq(v).to_numpy()) for v in ['Female','Male']]:
 sub.append(dict(dataset='GSE246622',variable=variable,level=level,age_cutoff='90 days' if variable=='Age' else '',n_positive=int(g[kk].sum()),n_reference=int((~g[kk]).sum()),**welch(Y.Total.to_numpy()[kk],g[kk])))
sub=pd.DataFrame(sub);sub['BH4_within_cohort']=padjust(sub.p_value,'BH',4);write(replace(sub,'26_age_sex_subgroup_effects.csv'),'26_age_sex_subgroup_effects.csv')

d=pd.read_csv(OUT/'models/GSE246622_gene_comparison_design.csv');X=d.iloc[:,1:].to_numpy();names=d.columns[1:].tolist();lookup=m.reset_index().set_index('sample_id')['index'];loc=lookup.loc[d.sample_id].to_numpy();ZZ=o['x'].iloc[:,loc].to_numpy().T
cI=np.zeros(X.shape[1]);cI[names.index('groupMild')]=1;cS=np.zeros(X.shape[1]);cS[names.index('groupSevere')]=1;cS[names.index('groupMild')]=-1;C=np.vstack([cI,cS]);coef=C@np.linalg.pinv(X)@ZZ;rho=float(stats.spearmanr(coef[0],coef[1]).statistic)
gene_map=pd.DataFrame({'human_ensembl':o['x'].index,'display_symbol':core.loc[o['x'].index,'display_symbol'].to_numpy(),'frozen_direction':core.loc[o['x'].index,'frozen_direction'].to_numpy(),'Mild_minus_Healthy_adjusted_log_effect':coef[0],'Severe_minus_Mild_adjusted_log_effect':coef[1]});write(gene_map,'27_infection_severity_gene_map.csv')
gg=m.group.to_numpy()[loc];indices=[np.flatnonzero(gg==v) for v in ['Healthy','Mild','Moderate','Severe']];boot=[]
for k in range(1000):
 sample=np.concatenate([rng.choice(v,len(v),replace=True) for v in indices]);w=np.bincount(sample,minlength=len(X));xwx=X.T@(X*w[:,None])
 try:W=C@np.linalg.solve(xwx,X.T);bb=(W*w[None,:])@ZZ;rr=float(stats.spearmanr(bb[0],bb[1]).statistic)
 except np.linalg.LinAlgError:rr=np.nan
 boot.append(dict(draw=k+1,Spearman_rho=rr))
boot=pd.DataFrame(boot);write(boot,'28_patient_bootstrap_gene_map.csv');valid=boot.Spearman_rho.dropna();ci=np.quantile(valid,[.025,.975])
mild=indices[1];healthy=indices[0];severe=indices[3];covcols=[0]+[i for i,n in enumerate(names) if n.startswith('ns(') or n.startswith('sex')]
def disjoint_effect(refids,posids):
 ids=np.r_[refids,posids];A=np.column_stack([np.ones(len(ids)),np.r_[np.zeros(len(refids)),np.ones(len(posids))],X[ids][:,covcols[1:]]]);return (np.linalg.pinv(A)@ZZ[ids])[1]
split=[]
for k in range(200):
 a,b=np.array_split(rng.permutation(mild),2);r1=stats.spearmanr(disjoint_effect(healthy,a),disjoint_effect(b,severe)).statistic;r2=stats.spearmanr(disjoint_effect(healthy,b),disjoint_effect(a,severe)).statistic
 split.append(dict(split=k+1,rho_first_half_infection=r1,rho_second_half_infection=r2,mean_swapped_rho=(r1+r2)/2))
split=pd.DataFrame(split);write(split,'29_disjoint_Mild_split_sensitivity.csv')
relation=dict(n_genes=len(o['x']),Spearman_rho=rho,patient_bootstrap_CI95_low=float(ci[0]),patient_bootstrap_CI95_high=float(ci[1]),bootstrap_draws=1000,valid_draws=len(valid),correlated_gene_P_not_computed=True,shared_Mild_samples_resampled_jointly=True,disjoint_splits=200,disjoint_split_median_rho=float(split.mean_swapped_rho.median()),disjoint_split_min_rho=float(split.mean_swapped_rho.min()),disjoint_split_max_rho=float(split.mean_swapped_rho.max()),covariates='Rebuilt natural age spline3df and sex; fixed within resampling',resampling_unit='Profile; real source patient IDs unique where available')
dump(OUT/'gene_relationship_summary.json',relation)
dump(OUT/'analysis_summary.json',dict(population='Author-defined otherwise-healthy',reference_profiles=262,acute_profiles=206,strict_n=115,coverage=scov,common_members=len(common_genes),gene_relationship=relation,other_raw_results_reused=True))
print(effects[['dataset','effect','CI95_low','CI95_high','p_value','primary_Holm_P','resampling_P','resampling_Holm_P']].to_string(index=False));print(fun[fun.branch.isin(['Neutrophil degranulation','TNF / NF-kB'])][['dataset','branch','effect','CI95_low','CI95_high','p_value','BH_adjusted_P']].to_string(index=False));print(json.dumps(relation,indent=2))
