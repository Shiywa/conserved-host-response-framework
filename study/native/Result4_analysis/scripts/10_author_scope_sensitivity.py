from pathlib import Path
from datetime import datetime,timezone
from itertools import combinations
from math import comb
import ast,json
import h5py,numpy as np,pandas as pd
from scipy import sparse,stats
OUT=(Path(__import__('os').environ['CRRM_WORKSPACE']) / 'Result4_analysis');RESOURCE=Path(str(Path(__import__('os').environ.get('CRRM_RESOURCE_ROOT','external_resource')) / ''))
plan=json.loads((OUT/'parameters/analysis_plan_v1.json').read_text());core=pd.read_csv(OUT/'parameters/frozen_CRRM388.csv');sets=json.loads((OUT/'parameters/score_gene_sets.json').read_text())
meta=pd.read_csv(OUT/'data/profile_metadata.csv');original=pd.read_csv(OUT/'results/profile_scores.csv');gp=meta[meta.source_label.eq('Infant_PBMC')&meta.view.eq('canonical')&meta.unit_resolved].sort_values('profile_id').copy()
cells=pd.read_parquet(OUT/'data/atlas_cell_metadata.parquet',columns=['cell_uid','source_label','expression_row_index','author_intent_sensitivity_included'])
cells=cells[cells.source_label.eq('Infant_PBMC')]
mp=pd.read_parquet(OUT/'data/profile_row_membership.parquet');mp=mp[mp.source_label.eq('Infant_PBMC')&mp.view.eq('canonical')&mp.profile_id.isin(gp.profile_id)]
mp=mp.merge(cells[['cell_uid','author_intent_sensitivity_included']],on='cell_uid',validate='one_to_one');mp=mp[mp.author_intent_sensitivity_included.astype(str).str.lower().eq('true')]
assert len(mp)==222725,('Unexpected author-intent scope',len(mp))
local={v:i for i,v in enumerate(gp.profile_id)}
var=pd.read_parquet(RESOURCE/'DATA/curated_metadata/features/Infant_PBMC_var.parquet')
with np.load(RESOURCE/'DATA/curated_metadata/features/gene_availability/Infant_PBMC.npz') as z:av=z['availability'].all(axis=0)
idx={g:i for i,g in enumerate(var.native_gene_id.astype(str).str.replace(r'\.\d+$','',regex=True))}
j=np.array([idx.get(g,-1) for g in core.human_ensembl]);measured=np.array([k>=0 and av[k] for k in j]);pos=np.flatnonzero(measured);host=var.host_gene_flag.to_numpy(bool)
with h5py.File(RESOURCE/'DATA/expression/Infant_PBMC.h5ad','r') as h:
    n,p=map(int,h['X'].attrs['shape']);ptr=h['X/indptr'][:];m=np.full(n,-1,dtype=int);m[mp.expression_row_index.to_numpy()]=mp.profile_id.map(local).to_numpy()
    sums=np.zeros((len(gp),388));counts=np.zeros(len(gp),dtype=int)
    for a in range(0,n,10000):
        b=min(a+10000,n);lo,hi=int(ptr[a]),int(ptr[b]);x=sparse.csr_matrix((h['X/data'][lo:hi],h['X/indices'][lo:hi],ptr[a:b+1]-lo),shape=(b-a,p));total=np.asarray(x[:,host].sum(axis=1)).ravel();r=np.flatnonzero((m[a:b]>=0)&(total>0))
        if not len(r):continue
        selected=x[:,j[pos]].toarray().astype(float);np.divide(selected*10000,total[:,None],out=selected,where=total[:,None]>0);np.log1p(selected,out=selected)
        agg=sparse.csr_matrix((np.ones(len(r)),(m[a:b][r],r)),shape=(len(gp),b-a));sums[:,pos]+=agg@selected;counts+=np.bincount(m[a:b][r],minlength=len(gp))
        if b==n:print('Author-scope extraction complete:',int(counts.sum()),'cells',flush=True)
means=np.full_like(sums,np.nan);np.divide(sums,counts[:,None],out=means,where=counts[:,None]>0)
assert counts.sum()==222725 and np.any(counts>=30)
refs=pd.read_csv(OUT/'data/gene_reference_parameters.csv');rr=refs[refs.source_label.eq('Infant_PBMC')].set_index('human_ensembl').reindex(core.human_ensembl);variable=rr.reference_variable.to_numpy();sign=np.where(core.frozen_direction.eq('up'),1.,-1.);zz=np.full_like(means,np.nan);zz[:,variable]=(means[:,variable]-rr.reference_mean.to_numpy()[variable])/rr.reference_SD.to_numpy()[variable]*sign[variable]
gp['cells']=counts
for target,genes in sets.items():
    mask=variable&core.human_ensembl.isin(genes).to_numpy();gp[target]=np.mean(zz[:,mask],axis=1) if mask.sum()>=5 else np.nan
gp.to_csv(OUT/'results/PBMC_author_scope_profile_scores.csv',index=False)
scores=gp;permutation_cache={};src=(OUT/'scripts/04_score_and_test_profiles.py').read_text();tree=ast.parse(src)
for name in ['unpaired_allocations','infer','contrast_rows','bh']:
    node=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name==name);exec(ast.get_source_segment(src,node))
rows=[]
for spec in plan['primary_comparisons']:
    if spec['source']!='Infant_PBMC':continue
    for target in ['CRRM388','Outside IFN lists','Neutrophil degranulation']:
        row,_=contrast_rows(spec,target);row['scope']='Predeclared author-intent cell subset; primary gene scaling fixed';rows.append(row)
r=pd.DataFrame(rows)
assert r.status.eq('OK').all() and r.case_units.ge(3).all() and r.control_units.ge(3).all()
for _,ix in r.groupby('score').groups.items():r.loc[ix,'adjusted_P']=bh(r.loc[ix,'P'])
r.to_csv(OUT/'results/PBMC_author_scope_effects.csv',index=False)
(OUT/'audit/author_scope_sensitivity_receipt.json').write_text(json.dumps(dict(status='PASS',completed_at=datetime.now(timezone.utc).isoformat(),retained_known_HC_RSV_cells=int(counts.sum()),original_known_HC_RSV_cells=229699,cell_flag='author_intent_sensitivity_included',original_scaling_preserved=True,author_published_DE_population_independently_confirmed=False,interpretation='Uses the descriptor author-intent sensitivity flag; it does not assert exact reconstruction of the published DE cell set.'),indent=2)+'\n')
print(r[['comparison_id','score','case_units','control_units','estimate','adjusted_P']].to_string(index=False))
