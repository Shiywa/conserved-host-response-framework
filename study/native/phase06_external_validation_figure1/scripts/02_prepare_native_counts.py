#!/usr/bin/env python3
import csv,gzip,hashlib,json,shutil,os
from pathlib import Path
from datetime import datetime,timezone
import numpy as np,pandas as pd
R=Path(str(Path(__import__('os').environ['CRRM_WORKSPACE']) / 'phase06_external_validation_figure1'))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,x):p.write_text(json.dumps(x,ensure_ascii=False,indent=2)+'\n')
manifest=pd.read_csv(R/'tables/01_external_sample_manifest.csv',keep_default_na=False)
manifest=manifest[manifest.included_primary.astype(str).str.lower().isin(['true','1'])]
records=[];summaries=[];outputs=[]
cohorts=os.environ.get('CRRM_EXTERNAL_COHORTS','GSE270463,GSE309353').split(',')
for g in cohorts:
 m=manifest[manifest.dataset.eq(g)].sort_values(['donor_id','condition','culture_replicate']).copy();matrices=[]
 for row in m.to_dict('records'):
  p=Path(row['local_file']);receipt=json.loads(p.with_name(p.name+'.receipt.json').read_text());assert sha(p)==receipt['sha256']
  if g=='GSE270463':
   d=pd.read_csv(p,sep='\t',dtype={'Feature':str});n_source=len(d);duplicates=d[d.duplicated(keep=False)].copy();d=d.drop_duplicates();n_exact_copies=n_source-len(d)
   if len(duplicates):duplicates.to_csv(R/'data/processed'/f'{row["gsm"]}_exact_duplicate_rows.csv',index=False)
   features=d.Feature;counts=d.Count;column='Count';symbols=None
  else:
   d=pd.read_csv(p,sep='\t',header=None,names=['feature_id','display_symbol','raw_count'],dtype={'feature_id':str});n_source=len(d);n_exact_copies=0;features=d.feature_id;counts=d.raw_count;column='third column, raw gene count';symbols=d.display_symbol
  ids=features.str.replace(r'\.[0-9]+$','',regex=True)
  human=ids.str.fullmatch(r'ENSG\d+')
  assert not ids[human].duplicated().any(),(g,row['gsm'],'ambiguous native human gene key')
  c=pd.to_numeric(counts,errors='raise').to_numpy();assert np.isfinite(c).all() and (c>=0).all()
  s=pd.Series(c[human],index=ids[human],name=row['gsm']);matrices.append(s)
  records.append(dict(dataset=g,gsm=row['gsm'],file=str(p),sha256=sha(p),count_column=column,n_source_features=n_source,n_identical_redundant_rows_removed=n_exact_copies,n_unique_human_features=len(s),n_nonhuman_excluded=int((~human).sum()),n_fractional=int(np.count_nonzero(c[human]!=np.round(c[human]))),library_sum=float(c[human].sum())))
 first=matrices[0].index
 assert all(s.index.equals(first) for s in matrices),'All selected libraries must have identical feature ordering and keys'
 cm=pd.concat(matrices,axis=1);path=R/'data/processed'/f'{g}_native_counts.csv.gz';cm.index.name='feature_id';cm.to_csv(path,compression={'method':'gzip','mtime':0})
 m.to_csv(R/'data/processed'/f'{g}_selected_manifest.csv',index=False)
 outputs.append(dict(dataset=g,path=str(path),sha256=sha(path),manifest=str(R/'data/processed'/f'{g}_selected_manifest.csv'),manifest_sha256=sha(R/'data/processed'/f'{g}_selected_manifest.csv')))
 groups=m.groupby(['donor_id','condition']).size();assert m.donor_id.nunique()==3 and set(m.condition)=={'mock','RSV'} and groups.nunique()==1
 summaries.append(dict(dataset=g,n_selected_profiles=len(m),n_donors=3,profiles_per_donor_condition=int(groups.iloc[0]),n_native_human_features=len(cm),all_counts_integer=bool((cm.to_numpy()==np.round(cm.to_numpy())).all()),count_scale='raw counts',sample_time_h=int(m.time_h.iloc[0])))
pd.DataFrame(records).to_csv(R/'tables/02_count_file_audit.csv',index=False)
pd.DataFrame(summaries).to_csv(R/'tables/03_external_design_summary.csv',index=False)
save(R/'parameters/count_implementation_receipt.json',dict(declared_at=datetime.now(timezone.utc).isoformat(),analysis_plan_sha256=sha(R/'parameters/analysis_plan_v1.json'),cohorts_prepared=cohorts,other_cohort_results_may_exist=(R/'external_validation_summary.json').exists(),count_scale='raw',selected_columns_rule='Source-declared raw count column; CPM/RPKM excluded',exact_duplicate_policy='Remove byte-equivalent table rows once after parsing, never sum duplicated copies; all feature metadata and count/CPM/RPKM fields must agree. Remaining ambiguous human IDs cause an explicit error rather than selection by expression.',source_files=records,processed_matrices=outputs,cohort_designs=summaries))
shutil.copy2(Path(__file__),R/'scripts/02_prepare_native_counts.py')
print(json.dumps(summaries,indent=2))
