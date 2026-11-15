"""Panel-level result ledger and long-form source data for all four revised figures."""
from pathlib import Path
import json
import numpy as np,pandas as pd
BASE=(Path(__import__('os').environ['CRRM_WORKSPACE']) / 'clinical_otherwisehealthy_262_r1')
def main(n):return pd.read_csv(BASE/'main/tables'/n,float_precision='round_trip')
def inc(n):return pd.read_csv(BASE/'incremental/tables'/n,float_precision='round_trip')
records=[];sources={k:[] for k in ['Fig3','FigS3','FigS4']}
def source(fig,panel,kind,d):
 for r in pd.DataFrame(d).to_dict('records'):sources[fig].append(dict(Figure=fig,Panel=panel,Record=kind,**r))
def result(fig,panel,comparison,target,r,*,study='',np_=np.nan,nr=np.nan,nt=np.nan,nm=np.nan,den=np.nan,p='p_value',q=None,family='Unadjusted/descriptive',role='Reported',interval='Pointwise95%tCI',effect='effect',lo='CI95_low',hi='CI95_high',**extra):
 r=dict(r);records.append(dict(Figure=fig,Panel=panel,Study=study,Comparison=comparison,Target=target,
  N_positive=r.get('n_positive',np_),N_reference=r.get('n_reference',nr),N_total=r.get('n',nt),
  Measured_members=nm,Scoring_denominator=den,Estimate=r.get(effect,np.nan),CI95_low=r.get(lo,np.nan),CI95_high=r.get(hi,np.nan),
  Interval_meaning=interval,P_raw=r.get(p,np.nan) if p else np.nan,P_adjusted=r.get(q,np.nan) if q else np.nan,
  Correction_family=family,Reporting_role=role,**extra))
cov=main('06_coverage_summary.csv').set_index('dataset');scores=main('07_sample_scores.csv');fx=main('09_total_clinical_effects.csv')
group=main('14_group_score_summary.csv').query("dataset=='GSE246622'")
source('Fig3','A','Individual score',scores[scores.dataset.eq('GSE246622')]);source('Fig3','A','Group summary',group)
for r in group.to_dict('records'):result('Fig3','A',r['group'],'Total score distribution',r,study='GSE246622',nm=int(cov.loc['GSE246622','n_scored']),effect='mean',interval='Group mean pointwise95%tCI;boxes display median/IQR',p=None)
source('Fig3','B','Total contrast',fx)
for r in fx.to_dict('records'):
 primary=r['role']=='Primary';result('Fig3','B',r['positive']+' minus '+r['reference'],'Total',r,study=r['dataset'],nm=r['n_measured'],q='primary_Holm_P' if primary else None,family='Holm2 model P' if primary else 'Nominal hospitalization context',role=r['role'],P_resampling=r['resampling_P'],P_resampling_adjusted=r['resampling_Holm_P'],Resampling_family='Separate Holm2' if primary else 'Nominal')
d=main('11_directional_contributions.csv');source('Fig3','C','Weighted directional effect',d)
for r in d.to_dict('records'):
 n=int(cov.loc[r['dataset'],'n_up' if r['component']=='Up contribution' else 'n_down']);result('Fig3','C',r['positive']+' minus '+r['reference'],r['component'],r,study=r['dataset'],nm=n,den=int(cov.loc[r['dataset'],'n_scored']),q='BH_adjusted_P',family='BH4 strict contributions' if r['role']=='Primary' else 'BH2 hospitalization contributions',role=r['role'])
d=main('10_functional_clinical_effects.csv');source('Fig3','D','Annotation contrast',d)
for r in d.to_dict('records'):result('Fig3','D',r['positive']+' minus '+r['reference'],r['branch'],r,study=r['dataset'],nm=r['n_measured_branch'],q='BH_adjusted_P',family=r['adjustment_family'],role=r['role'])
g=main('27_infection_severity_gene_map.csv');rel=json.loads((BASE/'main/gene_relationship_summary.json').read_text());source('Fig3','E','Adjusted gene-effect point',g)
result('Fig3','E','Adjusted Mild–Healthy vs Severe–Mild gene effects','Spearman rho',rel,study='GSE246622',nt=262,nm=rel['n_genes'],effect='Spearman_rho',lo='patient_bootstrap_CI95_low',hi='patient_bootstrap_CI95_high',p=None,interval='1000group-stratified profile-bootstrap percentile95%CI',role='Gene-effect correspondence;genes are correlated features')
d=main('17_covariate_adjusted_effects.csv');source('Fig3','F','Same-case HC3 contrast',d)
for r in d.to_dict('records'):result('Fig3','F','Severe–Mild' if r['dataset']=='GSE246622' else 'Inpatient–Outpatient',r['model']+' total',r,study=r['dataset'],nm=int(cov.loc[r['dataset'],'n_scored']),q='adjusted_model_BH2_P' if r['model']=='adjusted' else None,family='BH2 adjusted total models' if r['model']=='adjusted' else 'Nominal same-case HC3 comparator',role='Secondary covariate sensitivity',interval='Pointwise95%HC3 residual-df tCI')
d=main('19_paired_recovery_effects.csv');rr=d[d.metric.eq('Total')].iloc[0].to_dict();source('Fig3','G','Paired scores',main('18_paired_recovery_values.csv'));source('Fig3','G','Paired effect',pd.DataFrame([rr]))
result('Fig3','G','Recovery minus acute','Total',rr,study='GSE69606',nt=15,nm=int(cov.loc['GSE69606','n_scored']),P_resampling=rr['exact_signflip_P'],role='Secondary paired recovery',interval='Pointwise95%paired-tCI')
source('FigS3','A','Age and sex',scores[scores.dataset.isin(['GSE246622','GSE105450'])][['dataset','sample_id','group','age_days','age_months','sex']])
for r in main('25_clinical_covariate_summary.csv').to_dict('records'):result('FigS3','A',r['group'],'Age median',r,study=r['dataset'],effect='age_median',lo='age_Q1',hi='age_Q3',p=None,interval='Age IQR;not a CI',role='Descriptive covariate distribution',Female_count=r['n_female'])
source('FigS3','B','Member coverage',main('02_member_coverage.csv'));source('FigS3','B','Annotation coverage',main('08_functional_branch_coverage.csv'))
for ds,r in cov.iterrows():result('FigS3','B','Source measurement','Total variable members',r,study=ds,nm=int(r.n_scored),effect='n_scored',p=None,interval='No interval',role='Measurement coverage')
d=main('21_common_member_sensitivity.csv');source('FigS3','C','Common-member contrast',d);source('FigS3','C','Full-member contrast',fx)
for r in d.to_dict('records'):result('FigS3','C',r['positive']+' minus '+r['reference'],'Common measured members',r,study=r['dataset'],nm=r['n_common'],role='Measurement sensitivity')
d=main('24_hospitalization_study_strata.csv');source('FigS3','D','Study-stratum contrast',d)
for r in d.to_dict('records'):result('FigS3','D','Inpatient minus outpatient',r['stratum'],r,study='GSE105450',nm=int(cov.loc['GSE105450','n_scored']),role='Secondary study-stratum context')
loo=main('13_patient_deletion_influence.csv');source('FigS3','E','Unit-deletion effect',loo);source('FigS3','E','Profile-bootstrap draw',main('28_patient_bootstrap_gene_map.csv'));source('FigS3','E','Disjoint-Mild split',main('29_disjoint_Mild_split_sensitivity.csv'))
for k,v in loo.groupby('contrast_key'):
 r=fx[fx.contrast_key.eq(k)].iloc[0].to_dict();r['CI95_low']=v.effect.min();r['CI95_high']=v.effect.max();result('FigS3','E',k,'Delete-one-profile influence',r,study=r['dataset'],nm=r['n_measured'],p=None,interval='Min–max influence range;not a CI',role='Descriptive influence',Omissions=len(v))
result('FigS3','E','Disjoint Mild halves;39/40 profiles','Median swapped-split rho',rel,study='GSE246622',nt=262,nm=rel['n_genes'],effect='disjoint_split_median_rho',lo='disjoint_split_min_rho',hi='disjoint_split_max_rho',p=None,interval='200split influence range;not a CI',role='Descriptive shared-reference sensitivity')
result('FigS3','E','Jointly resampled Mild group','Spearman rho',rel,study='GSE246622',nt=262,nm=rel['n_genes'],effect='Spearman_rho',lo='patient_bootstrap_CI95_low',hi='patient_bootstrap_CI95_high',p=None,interval='1000group-stratified profile-bootstrap percentile95%CI',role='Gene-effect correspondence sensitivity')
d=main('19_paired_recovery_effects.csv').query("metric!='Total'");source('FigS3','F','Paired annotation effect',d)
for r in d.to_dict('records'):result('FigS3','F','Recovery minus acute',r['metric'],r,study='GSE69606',nt=15,q='branch_BH7_P',family='BH7 paired annotations',P_resampling=r['exact_signflip_P'],P_resampling_adjusted=r['branch_signflip_BH7_P'],Resampling_family='Separate BH7 exact sign-flip P',role='Secondary paired recovery',interval='Pointwise95%paired-tCI')
d=inc('07_LCN2_ablation_effects.csv');source('FigS4','A','Annotation-omission contrast',d)
for r in d[d.dataset.isin(['GSE246622','GSE69606'])].to_dict('records'):result('FigS4','A','Severe minus Mild',r['target'],r,study=r['dataset'],nm=r['n_measured'],q='ablation_Holm2_P' if r['target']=='CRRM24' else None,family='Holm2 LCN2-omitted contrast' if r['target']=='CRRM24' else 'Nominal contextual subset',role='Exploratory annotation omission')
d=inc('08_conditional_clinical_models.csv');source('FigS4','B','Conditional model',d);nested=inc('21_all_nested_models_Holm8_summary.csv');source('FigS4','B','Additional nested-model correction',nested)
for r in nested.to_dict('records'):result('FigS4','B','Severe minus Mild',r['model'],r,study=r['dataset'],q='nested_Holm8_P',family='Additional Holm8 nested models',role='Exploratory conditional association',interval='Pointwise95%HC3 residual-df tCI',Focal_Holm2_P=r['strict_Holm2_P'])
d=inc('10_every_gene_omission_influence.csv');source('FigS4','C','Additional-member omission',d)
models=inc('08_conditional_clinical_models.csv')
for ds in ['GSE246622','GSE69606']:
 v=d[d.dataset.eq(ds)&d.model.eq('Clinical')];r=models[models.dataset.eq(ds)&models.target.eq('CRRM24')&models.model.eq('Clinical')].iloc[0].to_dict();r['CI95_low']=v.effect.min();r['CI95_high']=v.effect.max();result('FigS4','C','Severe minus Mild','LCN2-omitted;each additional member omitted',r,study=ds,p=None,interval='Min–max influence range;not a CI',role='Exploratory member influence',Omissions=len(v))
d=inc('17_measured_composition_models.csv');source('FigS4','D','Same-case measured-composition model',d);source('FigS4','D','Eligible measured-composition covariates',inc('16_measured_composition_covariates.csv'))
for r in d.to_dict('records'):result('FigS4','D','Severe minus Mild',r['model'],r,study='GSE246622',q='composition_Holm2_P',family='Holm2 after-composition tests' if 'after_' in r['model'] else 'Nominal same-case baseline',role='Exploratory composition adjustment',interval='Pointwise95%HC3 residual-df tCI')
d=inc('20_descriptive_out_of_fold_comparison.csv');source('FigS4','E','All descriptive performance metrics',d);source('FigS4','E','Held-out probabilities',inc('19_leave_one_unit_out_predictions.csv'))
for r in d.to_dict('records'):
 for metric in ['logloss','Brier','AUC']:
  k='AUC_difference' if metric=='AUC' else metric+'_improvement';result('FigS4','E','Augmented versus baseline',metric,r,study=r['dataset'],effect=k,p=None,interval='No interval or performance test',role='Exploratory internal leave-one-profile evaluation',Baseline_metric=r['Baseline_'+metric],Augmented_metric=r['Augmented_'+metric])
d=inc('14_paired_recovery_ablation.csv');source('FigS4','F','Paired annotation/comparator change',d)
for r in d.to_dict('records'):result('FigS4','F','Recovery minus acute',r['metric'],r,study='GSE69606',nt=15,role=r['role'],interval='Pointwise95%paired-tCI',P_resampling=r.get('exact_signflip_P',np.nan))
s7=pd.read_csv(BASE/'source_data/Fig6_FigS7_SourceData.csv',float_precision='round_trip')
hostfx=pd.read_csv(BASE/'host/results/score_effects.csv');threshold=pd.read_csv(BASE/'host/results/threshold_sensitivities.csv');coverage=pd.read_csv(BASE/'host/results/measurement_coverage.csv')
payload=json.loads((BASE/'host/audit/figure_payload_otherwisehealthy_262_r1.json').read_text())
matrix=payload['supplement']['S7A']
for key in matrix['keys']:
 for method in matrix['methods']:
  if key=='BRD4_Airway' and method.startswith('Minimum'):
   result('FigS7','A',key,method,{},study='BRD4_hSAEC',p=None,family='Original method-specific BH12;inapplicable slot retained',role='Unchanged inapplicable cell-threshold check',interval='Unavailable/inapplicable');continue
  if method.startswith('Minimum'):d=threshold[threshold.key.eq(key)&threshold.score.eq('CRRM388')&threshold.method.eq(method)]
  else:d=hostfx[hostfx.key.eq(key)&hostfx.score.eq('CRRM388')&hostfx.method.eq('Primary' if method=='Exact test' else method)&hostfx.scope.eq('Intervention')]
  assert len(d)==1
  r=d.iloc[0].to_dict();cv=coverage[coverage.source.eq(r['source'])&coverage.method.eq('Primary' if method.startswith('Minimum') or method=='Exact test' else method)&coverage.score.eq('CRRM388')]
  result('FigS7','A',key,method,r,study=r['source'],np_=r['n_case'],nr=r['n_control'],nm=cv.scored_genes.iloc[0] if len(cv)==1 else np.nan,effect='estimate',lo='lower',hi='upper',p='exact_P' if method=='Exact test' else 'P',q='exact_adjusted_P' if method=='Exact test' else 'adjusted_P',family='Original method-specific BH12',role='Unchanged host sensitivity;'+r['status'],interval='Original pointwise95%Welch tCI;exact P separately defined')
for panel in ['S7B','S7C','S7D']:
 for r in s7[s7.Panel.eq(panel)&s7.estimate.notna()].to_dict('records'):
  cv=coverage[coverage.source.eq(r.get('source',''))&coverage.method.eq(r.get('method',''))&coverage.score.eq(r.get('score',''))]
  result('FigS7',panel[2:],str(r.get('key','')),str(r.get('score',r.get('method',''))),r,study=str(r.get('source','')),np_=r.get('n_case',np.nan),nr=r.get('n_control',np.nan),nm=cv.scored_genes.iloc[0] if len(cv)==1 else np.nan,effect='estimate',lo='lower',hi='upper',p='P',q='adjusted_P',family=str(r.get('family','Original family;see source data')),role='Unchanged host result',interval='Original pointwise95%Welch tCI')
for r in s7[s7.Panel.eq('S7D')&s7.full.notna()].to_dict('records'):result('FigS7','D',r['key'],r['score'],r,study=r['key'],nm=r['gene_omissions'],effect='full',lo='minimum',hi='maximum',p=None,family='No inferential family',role='Unchanged single-member influence',interval='Min–max influence range;not a CI',Omissions=r['gene_omissions'])
al=pd.read_csv(BASE/'host/results/clinical_gene_alignment.csv')
for r in al.to_dict('records'):result('FigS7','E',r['clinical_contrast'],'Spearman rho',r,study=r['clinical_contrast'],np_=36 if r['clinical_contrast'].startswith('GSE246622') else 8 if 'Severe' in r['clinical_contrast'] else 15,nr=79 if r['clinical_contrast'].startswith('GSE246622') else 9 if 'Severe' in r['clinical_contrast'] else 15,nm=r['members'],effect='rho',p='permutation_P',q='adjusted_P_BH3',family='BH3 clinical correspondence',role='Exploratory clinical correspondence',interval='No correlation CI;70whole-culture allocations',Culture_units=8,Label_allocations=70)
for fig,rows in sources.items():pd.DataFrame(rows).to_csv(BASE/'source_data'/f'{fig}_SourceData.csv',index=False)
table=pd.DataFrame(records)
for col in ['Comparison','Target','Interval_meaning','Correction_family','Reporting_role']:
 table[col]=table[col].str.replace('95%t','95% t',regex=False).str.replace('95%HC3','95% HC3',regex=False).str.replace('1000group','1,000 group',regex=False).str.replace('200split','200 split',regex=False).str.replace('70whole','70 whole',regex=False)
table.to_csv(BASE/'source_data/CRRM_clinical_panel_results.tsv',sep='\t',index=False)
panels=sorted(set(zip(table.Figure,table.Panel)))
assert len(panels)==24,(len(panels),panels)
(BASE/'audit/panel_source_inventory.json').write_text(json.dumps(dict(panels=panels,result_rows=len(table),source_rows={k:len(v) for k,v in sources.items()}),indent=2)+'\n')
print('Exported',len(table),'result rows covering24panels;', {k:len(v) for k,v in sources.items()})
