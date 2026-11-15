from pathlib import Path
from datetime import datetime,timezone
import ast,hashlib,json
from itertools import combinations
from math import comb
import h5py,numpy as np,pandas as pd
from scipy import sparse,stats

OUT=(Path(__import__('os').environ['CRRM_WORKSPACE']) / 'Result4_analysis');RESOURCE=Path(str(Path(__import__('os').environ.get('CRRM_RESOURCE_ROOT','external_resource')) / ''));META=RESOURCE/'DATA/curated_metadata'
plan=json.loads((OUT/'parameters/analysis_plan_v1.json').read_text());core=pd.read_csv(OUT/'parameters/frozen_CRRM388.csv');sets=json.loads((OUT/'parameters/score_gene_sets.json').read_text());scores=pd.read_csv(OUT/'results/profile_scores.csv');membership=pd.read_parquet(OUT/'data/profile_row_membership.parquet');refs=pd.read_csv(OUT/'data/gene_reference_parameters.csv');permutation_cache={}
src=(OUT/'scripts/04_score_and_test_profiles.py').read_text();tree=ast.parse(src)
for name in ['unpaired_allocations','infer','contrast_rows','bh']:
 node=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name==name);exec(ast.get_source_segment(src,node))
amendment=OUT/'parameters/measurement_robustness_amendment.json'
if not amendment.exists():amendment.write_text(json.dumps(dict(declared_at=datetime.now(timezone.utc).isoformat(),status='EXPLORATORY_AFTER_PRIMARY_CELLULAR_RESULTS',trigger='Contribution audit identified a substantial PIGR contribution despite detection in fewer than1%of PBMC B cells; evaluate measurement sensitivity uniformly, not a gene-set reconstruction',rank_method='Equal-member signed percentile rank of frozen genes within all commonly measured host features of each unit/cell-type count-pseudobulk; no gene-wise variance scaling',expression_supported='Gene detected in at least10%of cells in at least3 biological units, with conditions/time pooled without selecting contrast sign; same gene z-reference as primary',scope='All eight declared comparisons; original primary estimates and BH8/BH56 families unchanged',outcome_selection=False),indent=2)+'\n')
clarification=OUT/'parameters/measurement_robustness_clarification.json'
if not clarification.exists():
 clarification.write_text(json.dumps(dict(declared_at=datetime.now(timezone.utc).isoformat(),before_inspection_of_measurement_sensitivity_results=True,primary_cellular_results_already_known=True,rank_definition='(average transcriptome rank - 1)/(number of common host features - 1), centered at0.5',source_shared_mask='Union of expression-supported genes across the canonical cell identities already declared for inference in that source; use the SAME mask for all its cell identities, permitting composition sensitivity without cell-specific score definitions',gene_support='>=10% detected cells in >=3 units within any declared canonical identity; same primary gene-Z reference; no new core membership',primary_results_changed=False),indent=2)+'\n')
sign=np.where(core.frozen_direction=='up',1.,-1.);rankframes=[];supportframes=[];supportaudit=[]
for source in ['Berg_ALI','CD1c_DC','Infant_PBMC']:
 dest=OUT/'data'/('pseudobulk_host_counts_'+source+'.npz');gp=scores[scores.source_label==source].sort_values('profile_id').drop(columns=['Up contribution','Down contribution']).copy();gids=gp.profile_id.to_numpy();local={g:i for i,g in enumerate(gids)};mp=membership[membership.source_label==source]
 var=pd.read_parquet(META/'features'/(source+'_var.parquet'));npz=META/'features/gene_availability'/(source+'.npz')
 with np.load(npz,allow_pickle=False) as z:host=var.host_gene_flag.to_numpy(bool)&z['availability'].all(axis=0)
 if dest.exists():
  with np.load(dest,allow_pickle=False) as z:full=z['sum_host_counts'];assert np.array_equal(z['profile_ids'],gids)
 else:
  with h5py.File(RESOURCE/'DATA/expression'/(source+'.h5ad'),'r') as h:
   n,p=map(int,h['X'].attrs['shape']);indptr=h['X/indptr'][:];maps=[]
   for _,mm in mp.groupby('view'):
    ids=np.full(n,-1,dtype=int);ids[mm.expression_row_index.to_numpy()]=[local[g] for g in mm.profile_id];maps.append(ids)
   full=np.zeros((len(gp),int(host.sum())))
   for a in range(0,n,10000):
    b=min(n,a+10000);lo,hi=int(indptr[a]),int(indptr[b]);x=sparse.csr_matrix((h['X/data'][lo:hi],h['X/indices'][lo:hi],indptr[a:b+1]-lo),shape=(b-a,p))[:,host]
    for ids in maps:
     idx=ids[a:b];r=np.flatnonzero(idx>=0)
     if len(r):agg=sparse.csr_matrix((np.ones(len(r)),(idx[r],r)),shape=(len(gp),b-a));full+= (agg@x).toarray()
  np.savez_compressed(dest,sum_host_counts=full,profile_ids=gids)
 hostgenes=var.loc[host,'native_gene_id'].astype(str).str.replace(r'\.\d+$','',regex=True).to_numpy();jmap={g:i for i,g in enumerate(hostgenes)};encoded=np.array([g in jmap for g in core.human_ensembl]);pos=np.flatnonzero(encoded);ix=np.array([jmap[g] for g in core.human_ensembl[encoded]])
 rank=(stats.rankdata(full,axis=1,method='average')-1)/max(full.shape[1]-1,1)-.5
 for target,genes in sets.items():
  selected=encoded&core.human_ensembl.isin(genes).to_numpy();columns=[jmap[g] for g in core.human_ensembl[selected]];gp[target]=np.mean(rank[:,columns]*sign[selected],axis=1) if len(columns)>=5 else np.nan
 rankframes.append(gp)
 with np.load(OUT/'data'/('native_profiles_'+source+'.npz'),allow_pickle=False) as z:mean=z['mean_log_expression'];detect=z['expressing_cells']/z['cells'][:,None]
 original=scores.set_index('profile_id').loc[gids].reset_index().drop(columns=['Up contribution','Down contribution']);r=refs[refs.source_label==source].set_index('human_ensembl').reindex(core.human_ensembl);variable=r.reference_variable.to_numpy();zexpr=np.full_like(mean,np.nan);zexpr[:,variable]=(mean[:,variable]-r.reference_mean.to_numpy()[variable])/r.reference_SD.to_numpy()[variable]*sign[variable]
 source_supported=np.zeros(388,dtype=bool)
 primary_types={s['celltype'] for s in plan['primary_comparisons'] if s['source']==source and s['view']=='canonical'}
 for cell in primary_types:
  use=(original.analysis_celltype==cell)&(original.view=='canonical')&original.unit_resolved&(original.cells>=30)
  if source=='Infant_PBMC':use&=original.infection_status.isin(['Healthy','Mild RSV','Moderate RSV','Severe RSV'])
  unitvectors=[np.mean(detect[np.array(list(rows))],axis=0) for _,rows in original.loc[use].groupby('analysis_unit').groups.items()]
  if unitvectors:source_supported|=np.sum(np.array(unitvectors)>=.1,axis=0)>=3
 source_supported&=variable
 for cell,view in original[['analysis_celltype','view']].drop_duplicates().itertuples(index=False,name=None):
  selected=(original.analysis_celltype==cell)&(original.view==view)&original.unit_resolved&(original.cells>=30)
  if source=='Infant_PBMC':selected&=original.infection_status.isin(['Healthy','Mild RSV','Moderate RSV','Severe RSV'])
  unitvectors=[]
  for unit,rows in original.loc[selected].groupby('analysis_unit').groups.items():unitvectors.append(np.mean(detect[np.array(list(rows))],axis=0))
  supported=source_supported
  gg=original[(original.analysis_celltype==cell)&(original.view==view)].copy();rows=gg.index.to_numpy()
  for target,genes in sets.items():
   mask=supported&core.human_ensembl.isin(genes).to_numpy();gg[target]=np.mean(zexpr[rows][:,mask],axis=1) if mask.sum()>=5 else np.nan;supportaudit.append(dict(source_label=source,cell_identity=cell,view=view,score=target,supported_genes=int(mask.sum()),nominal_members=len(genes)))
  supportframes.append(gg)
 print(source,'rank and detection-supported profiles complete',flush=True)
allrank=pd.concat(rankframes,ignore_index=True);allsupport=pd.concat(supportframes,ignore_index=True);allrank.to_csv(OUT/'results/rank_score_profiles.csv',index=False);allsupport.to_csv(OUT/'results/detection_supported_score_profiles.csv',index=False);pd.DataFrame(supportaudit).to_csv(OUT/'results/detection_supported_gene_counts.csv',index=False)
out=[]
for method,data in [('Signed transcriptome rank',allrank),('Common expression-supported genes',allsupport)]:
 scores=data
 for spec in plan['primary_comparisons']:
  for target in ['CRRM388','IFN-list overlap','Outside IFN lists','Neutrophil degranulation','Degranulation without LCN2']:
   row,_=contrast_rows(spec,target);row['measurement_method']=method;row['status_of_analysis']='Exploratory measurement sensitivity; original primary analyses retained';out.append(row)
result=pd.DataFrame(out)
for (_,target),ix in result.groupby(['measurement_method','score']).groups.items():result.loc[ix,'descriptive_BH8_P']=bh(result.loc[ix,'P'])
result.to_csv(OUT/'results/measurement_method_sensitivities.csv',index=False)
print(result[(result.score=='CRRM388')][['comparison_id','measurement_method','estimate','P','descriptive_BH8_P','status']].to_string(index=False),flush=True)
