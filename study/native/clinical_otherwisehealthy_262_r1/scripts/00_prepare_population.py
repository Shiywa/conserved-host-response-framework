"""Freeze the authorized otherwise-healthy population and copy only reusable inputs."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib, json, shutil
import numpy as np
import pandas as pd

OUT=(Path(__import__('os').environ['CRRM_WORKSPACE']) / 'clinical_otherwisehealthy_262_r1')
ROOT=Path(str(Path(__import__('os').environ['CRRM_DATA_ROOT']) / ''))
REV=Path(__import__('os').environ['CRRM_WORKSPACE'])
P20=REV/'phase20_Result3_clinical_analysis_20261112'
P22=REV/'phase22_Result3_incremental_information_20261112'
R6=REV/'Result6_analysis'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(p,o):Path(p).write_text(json.dumps(o,ensure_ascii=False,indent=2,default=str)+'\n')
assert not (OUT/'parameters/population_revision_plan.json').exists(), 'Do not overwrite a frozen revision'
for d in ['parameters','audit','figures','submission_tables','source_data','reports']:
 (OUT/d).mkdir(exist_ok=True)
for d in ['main','incremental']:
 for sub in ['data/processed','tables','models','parameters','scripts','logs','figures']:
  (OUT/d/sub).mkdir(parents=True,exist_ok=True)
for sub in ['results','data/processed','figures','audit','scripts']:
 (OUT/'host'/sub).mkdir(parents=True,exist_ok=True)
inputs=[]
def cp(src,dst):
 dst=Path(dst);dst.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(src,dst)
 inputs.append(dict(path=str(src),sha256=sha(src),working_copy=str(dst.relative_to(OUT))))
for original,new in [(P20,'main'),(P22,'incremental')]:
 for folder in ['parameters','tables']:
  for p in (original/folder).glob('*'):
   if p.is_file():cp(p,OUT/new/folder/p.name)
 # Reuse processed expression and known annotations, never author preprocessing anew.
 for p in (original/'data/processed').glob('*.csv*'):
  cp(p,OUT/new/'data/processed'/p.name)
 for p in (original/'models').glob('*.csv'):
  if not p.name.startswith('GSE246622'):cp(p,OUT/new/'models'/p.name)
 # Original executed scripts are archived as provenance, not executed wholesale.
 for p in (original/'scripts').glob('*'):
  if p.is_file():cp(p,OUT/'audit/original_scripts'/new/p.name)
for p in (R6/'results').glob('*.csv'):cp(p,OUT/'host/results'/p.name)
for p in (R6/'figures').glob('*_layout.json'):cp(p,OUT/'host/figures'/p.name)
for p in (R6/'data/processed').glob('*'):
 if p.is_file() and 'BRD4' in p.name:cp(p,OUT/'host/data/processed'/p.name)

source=REV/'Methods_audit/GSE246622_profile_specimen_visit_patient_r3.tsv'
m=pd.read_csv(source,sep='\t',keep_default_na=False)
author=pd.read_csv(REV/'Methods_audit/RESCEU_author_expression_metadata_r2.tsv',sep='\t',keep_default_na=False)
assert len(m)==533 and author.Name.is_unique
raw=m.merge(author[['Name','Visit','DS_Invasiveventilation']],left_on='laboratory_id',right_on='Name',suffixes=('','_author'),validate='one_to_one')
acute=raw.Visit.eq('RSV_visit')&raw.DS_Invasiveventilation_author.isin(['Non_hospitalized','Hospitalized_wo_IV','Invasive_ventilation'])
healthy=raw.Visit.eq('Healthy')&raw.DS_Invasiveventilation_author.eq('Healthy')
assert acute.sum()==212 and healthy.sum()==56
admissible=acute|healthy
gradeok=raw.analysis_grade.isin(['Healthy','Mild','Moderate','Severe'])
included=admissible&gradeok
raw['included']=included
raw['author_population_member']=admissible
raw['population_basis']=np.where(admissible,'Author respiratory-support selection and published otherwise-healthy cohort designation','Outside the author-defined otherwise-healthy acute/healthy analysis set')
raw['exclusion_reason']=np.select([included,raw.source_visit_type.eq('RSV convalescence visit'),~admissible,~gradeok],['','Convalescent visit','Outside author-defined otherwise-healthy acute/healthy set','Missing deposited ReSVinet acute grade'],default='Unresolved')
raw['analysis_unit']='Expression profile; source patient ID retained where available'
raw['available_batch_field']='Not supplied in author phenotype metadata; author batch correction retained'
raw['clinical_group']=raw.analysis_grade
raw['age_days']=pd.to_numeric(raw.baseline_age_at_visit,errors='coerce')
raw['sex']=raw.Sex
cols=['geo_accession','laboratory_id','source_specimen_label','source_visit_id','source_visit_type','source_patient_id','patient_link_evidence','population_basis','author_population_member','clinical_group','age_days','sex','available_batch_field','included','exclusion_reason','analysis_unit','title','DS_Invasiveventilation_author','Resvinet_severity_groups']
raw[cols].to_csv(OUT/'source_data/GSE246622_population_eligibility.tsv',sep='\t',index=False)
raw.loc[included,cols].to_csv(OUT/'source_data/GSE246622_included_262.tsv',sep='\t',index=False)
raw.loc[~included,cols].to_csv(OUT/'audit/GSE246622_exclusions.tsv',sep='\t',index=False)
counts=raw.loc[included,'analysis_grade'].value_counts().to_dict()
assert counts=={'Moderate':91,'Mild':79,'Healthy':56,'Severe':36},counts
known=raw.loc[included&raw.source_patient_id.ne(''),'source_patient_id']
assert not known.duplicated().any(), 'Known repeated patients require grouped withholding'
ids=raw.loc[included,'laboratory_id'].tolist()
for block,exprsuffix in [('main','log_expression'),('incremental','target_log_expression')]:
 for ds in ['GSE246622','GSE246622_RMA']:
  p=OUT/block/'data/processed'/f'{ds}_metadata.csv'
  meta=pd.read_csv(p);meta=meta[meta.sample_id.isin(ids)].copy()
  meta['included']=True;meta['exclusion_reason']='';meta['reference_for_scaling']=True
  add=raw.loc[included].set_index('laboratory_id')
  for name in ['source_patient_id','source_visit_id','population_basis']:
   meta[name]=add.loc[meta.sample_id,name].to_numpy()
  meta.to_csv(p,index=False)
  xp=OUT/block/'data/processed'/f'{ds}_{exprsuffix}.csv.gz'
  x=pd.read_csv(xp,index_col=0,float_precision='round_trip');x=x.loc[:,meta.sample_id]
  x.to_csv(xp,compression='gzip',index_label='human_ensembl')
  # Old z values cannot enter any new analysis.
  for old in (OUT/block/'data/processed').glob(f'{ds}*gene_z.csv.gz'):old.unlink()
# Do not carry obsolete probe-level exports or old assembled result sheets as current results.
for p in (OUT/'main/data/processed').glob('*native*'):p.unlink()
for folder in ['main','incremental']:
 for p in (OUT/folder/'tables').glob('S[34]*.csv'):p.unlink()
for src in [source,REV/'Methods_audit/RESCEU_author_expression_metadata_r2.tsv',
 P22/'data/public/RESCEU_downstreamAnalysis.R',P22/'data/public/RESCEU_datTraits_counts.txt']:
 inputs.append(dict(path=str(src),sha256=sha(src)))
plan=dict(version='otherwisehealthy_262_r1',frozen_at_utc=datetime.now(timezone.utc).isoformat(),
 authorized_population='Author-defined otherwise-healthy GSE246622 acute/healthy set with deposited ReSVinet grade',
 counts=counts,total_profiles=262,acute_profiles=206,strict_Severe=36,strict_Mild=79,
 all_patient_identity_verified=False,linked_source_patients=int(len(known)),known_patient_duplicates=0,
 keep_frozen_members=388,up=280,down=108,frozen_annotation_subsets=7,
 scoring_reference='All262profiles, separately for author batch-corrected and RMA expression',
 old313_population='Internal history only; not a displayed primary, expanded or sensitivity population',
 spline='R splines::ns(df=3), refit on115severity profiles for score models and on262for four-group gene model; inherit the115-case basis within exploratory folds',
 resampling=dict(seed=20261120,permutations=19999,gene_bootstrap=1000,disjoint_Mild_splits=200),
 multiplicity='Retain all original families; update shared adjusted P using new whole-blood raw P; no r3 Holm3 population-restriction family',
 other_studies='Reuse original scores,effects,raw P and resampling results; recalculate shared adjusted P only',
 narratives='No edits to Word,Results,Abstract,Highlights,Discussion,formal legends or Methods body',
 inputs=inputs)
dump(OUT/'parameters/population_revision_plan.json',plan)
dump(OUT/'parameters/population_revision_receipt.json',dict(plan_sha256=sha(OUT/'parameters/population_revision_plan.json'),new_analysis_started=False))
print(json.dumps({'population_counts':counts,'linked_real_patient_IDs':len(known),'strict_comparison':'36 Severe vs79Mild','new_root':str(OUT)},ensure_ascii=False))
