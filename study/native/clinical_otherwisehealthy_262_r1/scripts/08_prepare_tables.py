"""Reuse established table structures; rebuild their rows from the new outputs."""
from pathlib import Path
import csv,json,re,runpy
import numpy as np,pandas as pd
BASE=(Path(__import__('os').environ['CRRM_WORKSPACE']) / 'clinical_otherwisehealthy_262_r1')
REV=Path(str(Path(__import__('os').environ['CRRM_WORKSPACE']) / ''))
P20=REV/'phase20_Result3_clinical_analysis_20261112';P22=REV/'phase22_Result3_incremental_information_20261112'
s=(P20/'scripts/08_prepare_clinical_key_tables.py').read_text()
s=s.replace("'54 / 87 acute samples',369,r", "f'{int(r.n_positive)} / {int(r.n_reference)} acute profiles',int(r.n_measured),r")
s=s.replace("f'351 members:255 up,96 down; full-score correlation {r.full_score_correlation:.6f}; descriptive sensitivity'", "f'{int(r.n_common)} members:{int(r.n_common_up)} up,{int(r.n_common_down)} down; full-score correlation {r.full_score_correlation:.6f}; descriptive sensitivity'")
s=s.replace("'313 clinical profiles;369 genes'", "str(len(pd.read_csv(OUT/'data/processed/GSE246622_metadata.csv')))+' clinical profiles;'+str(r['n_genes'])+' genes'")
s=s.replace("source Mild group n=87", "source Mild group n=79")
s=s.replace('patient bootstrapping','profile bootstrapping').replace('patient-bootstrap','profile-bootstrap').replace('patient resamples','profile resamples').replace('shared Mild patients','shared Mild profiles')
s=s.replace("(OUT/'parameters/workbook_spec.json').write_text", "for item in spec:item['notes']=[]\n(OUT/'parameters/workbook_spec.json').write_text")
p=BASE/'main/scripts/08_prepare_clinical_key_tables.py';p.write_text(s);runpy.run_path(str(p))
s=(P22/'scripts/08_prepare_key_workbook.py').read_text().replace('Each patient has one held-out probability per model.','Each profile has one held-out probability per model; complete individual patient linkage is unavailable.')
p=BASE/'incremental/scripts/08_prepare_key_workbook.py';p.write_text(s);runpy.run_path(str(p))
src=(REV/'phase24_streamlined_tables_and_citations_20261112/scripts/01_build_streamlined_tables.py').read_text()
head=src[src.index('HEADERS='):src.index('primary=read_table')]
body=src[src.index('s3=[]'):src.index('summary_specs=')]
env={'Path':Path,'csv':csv,'json':json,'re':re,'inputs':[],'provenance':[], 'specs':[]}
def read_table(table,sheet):
 folder='main' if table=='S3' else 'incremental'
 with (BASE/folder/'tables'/(sheet+'.csv')).open() as f:return list(csv.DictReader(f))
env['read_table']=read_table
exec(compile(head+'\n'+body,'<established-streamlined-table-rows>','exec'),env)
s3,s4=env['s3'],env['s4'];headers=env['HEADERS']
for row in s3+s4:
 for k,v in row.items():
  if isinstance(v,str):
   row[k]=v.replace('patient-bootstrap','profile-bootstrap').replace('patient resamples','profile resamples').replace('donor pairs','participant pairs').replace('clinical subjects','clinical participants')
# Keep every original metric: add four compact AUC rows to the existing single S4 sheet.
perf=pd.read_csv(BASE/'incremental/tables/20_descriptive_out_of_fold_comparison.csv')
for r in perf.itertuples():
 row={k:'' for k in headers};row.update({'Analysis':'Held-out comparison | '+r.role,'Study':r.dataset,
 'Target or contrast':'Reference model versus model adding the LCN2-omitted subset','Model or metric':'AUC difference (augmented minus baseline)',
 'Analysis units':str(r.n)+' profiles','Estimate':r.AUC_difference,'Interval meaning':'Descriptive AUC difference; no interval',
 'Correction family':'No performance P value or multiplicity test','Notes':f'Baseline AUC = {r.Baseline_AUC:.17g}; augmented AUC = {r.Augmented_AUC:.17g}; exploratory internal evaluation'})
 s4.append(row)
def write_csv(name,rows,cols):
 with (BASE/'submission_tables'/name).open('w',newline='') as f:
  wr=csv.DictWriter(f,fieldnames=cols);wr.writeheader();wr.writerows(rows)
write_csv('S3A_Summary.csv',s3,headers);write_csv('S4_Key_results.csv',s4,headers)
gene=pd.read_csv(BASE/'main/tables/S3C_Gene_effects.csv');gene.to_csv(BASE/'submission_tables/S3B_Gene_effects.csv',index=False)
spec=[dict(workbook='TableS3_临床关联与恢复.xlsx',sheet='S3A_Summary',file='S3A_Summary.csv',title='Table S3A. Clinical CRRM results: score and annotation associations, paired recovery, and sensitivities'),
 dict(workbook='TableS3_临床关联与恢复.xlsx',sheet='S3B_Gene_effects',file='S3B_Gene_effects.csv',title='Table S3B. Clinical gene effects for all 388 frozen CRRM members'),
 dict(workbook='TableS4_稳健性与增量.xlsx',sheet='S4_Key_results',file='S4_Key_results.csv',title='Table S4. Exploratory member-omission, conditional association, composition, and internal evaluation results')]
(BASE/'parameters/submission_workbooks.json').write_text(json.dumps(spec,ensure_ascii=False,indent=2)+'\n')

# Update the statistical-design table, without changing any narrative Methods.
with (REV/'Methods_audit/CRRM_Methods_statistical_design_r3.tsv').open() as f:
 title=next(f).strip();reader=csv.DictReader(f,delimiter='\t');cols=reader.fieldnames;design=list(reader)
design=[r for r in design if r['Analysis question']!='Comorbidity population restriction']
population={k:'' for k in cols};population.update({'Analysis question':'Whole-blood analysis population',
 'Comparison':'GSE246622 author-defined otherwise-healthy acute/healthy study set',
 'Inference unit':'Expression profile; real participant/visit identifiers retained when available',
 'Test or model and interval':'Author cohort and respiratory-support selection,followed by deposited ReSVinet-grade eligibility',
 'Hypothesis family and size':'No additional test family', 'Adjustment or unadjusted status':'Not a hypothesis test',
 'Corresponding figures and tables':'Fig. 3; Figs. S3–S4,S7E; Tables S3–S4,S7', 'Reporting role':'Analysis population definition',
 'Necessary notes':'56 Healthy,79 Mild,91 Moderate,36 Severe:262 profiles,206 acute. Six acute records in the author268-record set lacked the required grade. Other clinical studies are not relabeled otherwise-healthy. Full individual patient deduplication is not asserted.'})
reference={k:'' for k in cols};reference.update({'Analysis question':'Whole-blood score reference',
 'Comparison':'Author batch-corrected expression; separate RMA-only view with identical participants/profiles',
 'Inference unit':'All262 retained expression profiles within each source view',
 'Test or model and interval':'Gene mean and sample SD,SD>10^−8,fixed discovery signs and actual measured-member denominators',
 'Hypothesis family and size':'No additional test family','Adjustment or unadjusted status':'Not a hypothesis test',
 'Corresponding figures and tables':'Fig. 3; Figs. S3–S4; Tables S3–S4','Reporting role':'Scoring reference definition',
 'Necessary notes':'369/388 variable core members in each view;351 common members across clinical studies. Omissions and composition analyses inherit new parent parameters; comparator scores use the same population. Rank-score testing was not part of the executed whole-blood clinical workflow and is not newly introduced.'})
design[10:10]=[population,reference]
spline=json.loads((BASE/'parameters/age_spline_parameters.json').read_text())
for r in design:
 q=r['Analysis question']
 if q=='Strict clinical total-score association':r['Necessary notes']='GSE246622:36 Severe/79 Mild in the author-defined otherwise-healthy set. GSE69606:8 Severe/9 Mild. Individual source IDs are checked where available; hospitalization remains a separate phenotype.'
 if q=='Covariate-adjusted total scores':r['Necessary notes']='Whole blood uses the115-case three-df natural age spline (knots55/177 days,boundaries10/356)+sex; hospitalization retains month-age+sex+study stratum. PBMC individual covariates are unavailable.'
 if q=='Infection/severity gene-effect correspondence':r['Necessary notes']='One262-profile four-group age/sex model; age knots60/186 days,boundaries10/363. Shared79-profile Mild group is bootstrapped jointly. The200 Mild splits have39/40-profile halves and swapped roles; their ranges are not CIs.'
 if q=='Internal incremental prediction evaluation':r['Necessary notes']='Whole-blood115-fold evaluation:114 training profiles;261-profile gene reference (56Healthy+91Moderate+114training). The rebuilt115-case spline basis and author batch correction stay fixed across folds; predictor scaling uses training cases only. Full patient-disjointness is not asserted.'
 if q=='Measured and expression-proxy clinical composition':r['Necessary notes']='Exact intersection with the new source set retains60 linked participants; targeted models use33 cases and complete models24. Three PCs are refitted without outcome labels on42 eligible complete records. Before/after models use identical cases.'
 if q=='Perturbation–clinical gene correspondence':r['Necessary notes']='Whole-blood vector:unadjusted log-expression Severe–Mild,36/79 profiles; other clinical vectors unchanged. Measured genes24/21/21 are features,not participants. Eight BRD4 cultures and all70 label allocations are retained.'
 for k in cols:r[k]=re.sub(r'(?<=[A-Za-z]),(?=\S)|,(?=[A-Za-z])',', ',r[k])
 for k in cols:
  for old,new in [('All262','All 262'),('author268-record','author 268-record'),('115-case','115-case'),('the115-case','the 115-case'),('One262-profile','One 262-profile'),('knots55/177','knots 55/177'),('boundaries10/356','boundaries 10/356'),('knots60/186','knots 60/186'),('boundaries10/363','boundaries 10/363'),('Shared79-profile','Shared 79-profile'),('The200 Mild','The 200 Mild'),('have39/40-profile','have 39/40-profile'),('Whole-blood115-fold','Whole-blood 115-fold'),('evaluation:114','evaluation: 114'),(';261-profile','; 261-profile'),('56Healthy+91Moderate+114training','56 Healthy + 91 Moderate + 114 training'),('retains60','retains 60'),('use33','use 33'),('models24','models 24'),('on42','on 42'),('Eight BRD4 cultures and all70','Eight BRD4 cultures and all 70'),('24/21/21 are features,not participants','24/21/21 are features, not participants'),('UK388','UK 388')]:r[k]=r[k].replace(old,new)
with (BASE/'submission_tables/TableS8_统计设计.tsv').open('w',newline='') as f:
 wr=csv.DictWriter(f,fieldnames=cols,delimiter='\t');wr.writeheader();wr.writerows(design)
(BASE/'submission_tables/TableS8_title.txt').write_text(title+'\n')

# Source data for S7: change only its three-statistic family and the whole-blood points.
source=pd.read_csv(REV/'Result6_analysis/reports/Fig6_FigS7_SourceData_component_focus_r2.csv',float_precision='round_trip')
al=pd.read_csv(BASE/'host/results/clinical_gene_alignment.csv');pts=pd.read_csv(BASE/'host/results/clinical_gene_alignment_points.csv')
for rowid,r in source[source.Panel.eq('S7E')].iterrows():
 if r.Record=='Clinical vector statistic':
  v=al[al.clinical_contrast.eq(r.clinical_contrast)].iloc[0]
  for k in al.columns:source.loc[rowid,k]=v[k]
 elif r.Record=='Gene-effect point' and r.clinical_contrast=='GSE246622_Severe_Mild':
  v=pts[pts.clinical_contrast.eq(r.clinical_contrast)&pts.human_ensembl.eq(r.human_ensembl)].iloc[0]
  for k in pts.columns:source.loc[rowid,k]=v[k]
source.to_csv(BASE/'source_data/Fig6_FigS7_SourceData.csv',index=False)
print('Submission tables: S3',len(s3),'summary rows +388genes; S4',len(s4),'rows; S8',len(design),'design rules with header in line1.')
