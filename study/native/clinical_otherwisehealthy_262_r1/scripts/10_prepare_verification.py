"""Prepare independent replay and verify reuse of untouched cohorts/host panels."""
from pathlib import Path
import hashlib,json,shutil
import numpy as np,pandas as pd
BASE=(Path(__import__('os').environ['CRRM_WORKSPACE']) / 'clinical_otherwisehealthy_262_r1')
REV=Path(str(Path(__import__('os').environ['CRRM_WORKSPACE']) / ''))
P20=REV/'phase20_Result3_clinical_analysis_20261112';P22=REV/'phase22_Result3_incremental_information_20261112';R6=REV/'Result6_analysis'
A=BASE/'audit';(A/'model_replay').mkdir(exist_ok=True)
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
read=lambda p:pd.read_csv(p,float_precision='round_trip')
def dump(p,o):Path(p).write_text(json.dumps(o,ensure_ascii=False,indent=2,default=str)+'\n')

# No old-population plan is an active plan. Membership/mapping files stay fixed.
for block in ['main','incremental']:
 for name in ['clinical_analysis_plan_v1.json','exploration_plan_v1.json','declaration_receipt.json','functional_covariate_sensitivity_amendment.json','nested_multiplicity_reporting_addendum.json']:
  p=BASE/block/'parameters'/name
  if p.exists():
   dst=A/'prior_plans'/block/name;dst.parent.mkdir(parents=True,exist_ok=True);shutil.move(p,dst)

preserved=[]
for old,block in [(P20,'main'),(P22,'incremental')]:
 for p in (old/'tables').glob('*.csv'):
  new=BASE/block/'tables'/p.name
  if not new.exists() or p.name.startswith('S'):continue
  x,y=read(p),read(new)
  if 'dataset' not in x or 'dataset' not in y:continue
  keep=lambda d:d.loc[~d.dataset.isin(['GSE246622','GSE246622_RMA'])].reset_index(drop=True)
  a,b=keep(x),keep(y);assert len(a)==len(b),(p.name,len(a),len(b))
  cols=[c for c in x if c in y and not any(k in c for k in ['BH','Holm','adjusted_P'])]
  for c in cols:
   if pd.api.types.is_numeric_dtype(a[c]):
    assert np.allclose(a[c].to_numpy(dtype=float),pd.to_numeric(b[c],errors='coerce'),rtol=0,atol=1e-12,equal_nan=True),(p.name,c)
   else:assert a[c].fillna('').astype(str).tolist()==b[c].fillna('').astype(str).tolist(),(p.name,c)
  preserved.append(dict(block=block,file=p.name,other_study_rows=len(a),unchanged_raw_fields=len(cols)))
# Unchanged recovery and original clinical host scores are literal file copies.
for name in ['18_paired_recovery_values.csv','19_paired_recovery_effects.csv','20_recovery_initial_grade_strata.csv','24_hospitalization_study_strata.csv']:
 assert sha(P20/'tables'/name)==sha(BASE/'main/tables'/name)
assert sha(P22/'tables/14_paired_recovery_ablation.csv')==sha(BASE/'incremental/tables/14_paired_recovery_ablation.csv')
for p in (R6/'results').glob('*.csv'):
 if p.name not in ['clinical_gene_alignment.csv','clinical_gene_alignment_points.csv']:
  assert sha(p)==sha(BASE/'host/results'/p.name),p.name
oldlayout=json.loads((R6/'figures/FigureS7_component_focus_r3_layout.json').read_text())
newlayout=json.loads((BASE/'host/figures/FigureS7_component_focus_r3_layout.json').read_text())
unchanged_axes=[k for k in oldlayout['artists'] if not k.startswith('S7E')]
assert all(oldlayout['artists'][k]==newlayout['artists'][k] for k in unchanged_axes)
dump(A/'unchanged_study_and_host_checks.json',dict(status='PASS',other_study_rows=preserved,unchanged_S7_A_to_D_axes=unchanged_axes,Fig6_not_recomputed=True))

# Persist actual 262-profile comparator reference parameters and coverage.
refrows=[]
for ds in ['GSE246622','GSE246622_RMA']:
 p=BASE/'incremental/data/processed'/f'{ds}_target_log_expression.csv.gz';x=pd.read_csv(p,index_col=0,float_precision='round_trip')
 mu=x.mean(1);sd=x.std(1,ddof=1);valid=np.isfinite(sd)&sd.gt(1e-8)
 z=pd.read_csv(BASE/'incremental/data/processed'/f'{ds}_target_gene_z.csv.gz',index_col=0,float_precision='round_trip')
 rebuilt=x.loc[valid].sub(mu[valid],axis=0).div(sd[valid],axis=0)
 assert np.max(np.abs(rebuilt.to_numpy()-z.to_numpy()))<1e-12
 par=pd.DataFrame(dict(human_ensembl=x.index,reference_mean=mu,reference_sample_SD=sd,scorable=valid))
 par.to_csv(BASE/'incremental/parameters'/f'{ds}_comparator_gene_reference.csv',index=False)
 refrows.append(dict(view=ds,reference_profiles=x.shape[1],measured_target_genes=x.shape[0],variable_target_genes=int(valid.sum())))
dump(A/'reference_replay_checks.json',dict(status='PASS',reference_views=refrows,zero_imputation=False))

families=[];pairs=[]
def family(d,p,q,n,method,label):
 fid=f'F{len(families)+1:03d}';assert 0<len(d)<=n
 families.append(dict(family_id=fid,label=label,n=n,rows=len(d),method=method))
 for j,(pp,qq) in enumerate(zip(d[p],d[q])):pairs.append(dict(family_id=fid,position=j,n=n,method=method,P=pp,stored_adjusted_P=qq))
def table(block,name):return read(BASE/block/'tables'/name)
d=table('main','09_total_clinical_effects.csv');d=d[d.role.eq('Primary')]
for p,q in [('p_value','primary_Holm_P'),('resampling_P','resampling_Holm_P')]:family(d,p,q,2,'holm','Strict clinical totals:'+p)
for name,q,ns in [('10_functional_clinical_effects.csv','BH_adjusted_P',{'Primary':14,'Hospitalization sensitivity':7}),('11_directional_contributions.csv','BH_adjusted_P',{'Primary':4,'Hospitalization sensitivity':2})]:
 d=table('main',name)
 for role,n in ns.items():family(d[d.role.eq(role)],'p_value',q,n,'BH',name+':'+role)
for name,q,n in [('16_secondary_grade_comparisons.csv','BH3_within_cohort',3),('26_age_sex_subgroup_effects.csv','BH4_within_cohort',4)]:
 d=table('main',name);family(d[d.dataset.eq('GSE246622')],'p_value',q,n,'BH',name+':wholeblood')
d=table('main','17_covariate_adjusted_effects.csv');family(d[d.model.eq('adjusted')],'p_value','adjusted_model_BH2_P',2,'BH','Adjusted total scores')
d=table('main','30_covariate_adjusted_branch_sensitivity.csv');family(d,'p_value','BH14_adjusted_branch_sensitivity_P',14,'BH','Adjusted annotations')
d=table('incremental','08_conditional_clinical_models.csv')
for (target,model),a in d[d.strict_Holm2_P.notna()].groupby(['target','model']):family(a,'p_value','strict_Holm2_P',2,'holm','Conditional:'+target+model)
for name,q,n in [('07_LCN2_ablation_effects.csv','ablation_Holm2_P',2),('09_direct_branch_difference.csv','strict_Holm2_P',2),('17_measured_composition_models.csv','composition_Holm2_P',2),('18_disjoint_expression_proxy_models.csv','strict_Holm2_P',2),('21_all_nested_models_Holm8_summary.csv','nested_Holm8_P',8),('12_all25_member_clinical_effects.csv','BH50_P',50)]:
 d=table('incremental',name);a=d if name.startswith('12_') else d[d[q].notna()];family(a,'p_value',q,n,'BH' if n==50 else 'holm',name)
d=read(BASE/'host/results/clinical_gene_alignment.csv');family(d,'permutation_P','adjusted_P_BH3',3,'BH','Three original clinical correspondence tests')
pd.DataFrame(families).to_csv(A/'affected_family_registry.tsv',sep='\t',index=False)
pd.DataFrame(pairs).to_csv(A/'affected_P_value_pairs.tsv',sep='\t',index=False)

modelrows=[]
scores=table('main','07_sample_scores.csv').query("dataset=='GSE246622'").set_index('sample_id')
for r in table('main','17_covariate_adjusted_effects.csv').query("dataset=='GSE246622'").itertuples():
 d=read(BASE/'main/models'/f'GSE246622_{r.model}_design.csv');d['outcome_score']=scores.loc[d.sample_id,'score_Total'].to_numpy();p=A/'model_replay'/f'main_total_{r.model}.csv';d.to_csv(p,index=False);modelrows.append(dict(file=str(p),effect=r.effect,CI95_low=r.CI95_low,CI95_high=r.CI95_high,p_value=r.p_value,label='main_total_'+r.model))
for r in table('main','30_covariate_adjusted_branch_sensitivity.csv').query("dataset=='GSE246622'").itertuples():
 d=read(BASE/'main/models/GSE246622_adjusted_design.csv');d['outcome_score']=scores.loc[d.sample_id,'score_'+r.branch].to_numpy();p=A/'model_replay'/f'main_branch_{len(modelrows)}.csv';d.to_csv(p,index=False);modelrows.append(dict(file=str(p),effect=r.effect,CI95_low=r.CI95_low,CI95_high=r.CI95_high,p_value=r.p_value,label='main_branch_'+r.branch))
for name in ['08_conditional_clinical_models.csv','09_direct_branch_difference.csv','17_measured_composition_models.csv','18_disjoint_expression_proxy_models.csv']:
 d=table('incremental',name);d=d[d.dataset.isin(['GSE246622','GSE246622_RMA'])]
 for r in d.itertuples():modelrows.append(dict(file=str(BASE/'incremental'/r.design_file),effect=r.effect,CI95_low=r.CI95_low,CI95_high=r.CI95_high,p_value=r.p_value,label=r.dataset+'_'+r.target+'_'+r.model))
pd.DataFrame(modelrows).to_csv(A/'HC3_models_for_R.tsv',sep='\t',index=False)
welchrows=[]
s=table('main','07_sample_scores.csv')
for name,targetcol,scoreprefix in [('09_total_clinical_effects.csv',None,'score_Total'),('10_functional_clinical_effects.csv','branch',None),('11_directional_contributions.csv','component',None),('23_RMA_source_sensitivity.csv',None,'score_Total')]:
 d=table('main',name);d=d[d.dataset.isin(['GSE246622','GSE246622_RMA'])]
 for r in d.itertuples():
  v=s[s.dataset.eq(r.dataset)&s.group.isin(['Severe','Mild'])];target=scoreprefix if targetcol is None else 'score_'+getattr(r,targetcol)
  p=A/'model_replay'/f'welch_{len(welchrows)}.csv';pd.DataFrame(dict(high=v.group.eq('Severe').astype(int),score=v[target])).to_csv(p,index=False)
  welchrows.append(dict(file=str(p),effect=r.effect,CI95_low=r.CI95_low,CI95_high=r.CI95_high,p_value=r.p_value,label=name+':'+target))
pd.DataFrame(welchrows).to_csv(A/'Welch_models_for_R.tsv',sep='\t',index=False)
dump(A/'verification_scope.json',dict(affected_families=len(families),P_positions=len(pairs),HC3_models=len(modelrows),Welch_models=len(welchrows),unrelated_main_analyses_rerun=False))
print('Prepared independent R verification:',len(families),'families;',len(modelrows),'HC3 models;',len(welchrows),'Welch contrasts. Untouched raw effects and S7A–D verified.')
