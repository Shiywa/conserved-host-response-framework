"""Diagnose inherited epithelial mixture without relabeling the descriptor."""
from pathlib import Path
from datetime import datetime,timezone
import importlib.util,json,hashlib
import h5py,numpy as np,pandas as pd
from scipy import sparse,stats

OUT=(Path(__import__('os').environ['CRRM_WORKSPACE']) / 'Result6_analysis');RESOURCE=Path(str(Path(__import__('os').environ.get('CRRM_RESOURCE_ROOT','external_resource')) / ''))
def dump(p,x):p.write_text(json.dumps(x,indent=2,ensure_ascii=False)+'\n')
assert (OUT/'audit/scoring_receipt.json').exists()
amendment=dict(declared_at=datetime.now(timezone.utc).isoformat(),parent_plan_sha256=hashlib.sha256((OUT/'parameters/analysis_plan_v1.json').read_bytes()).hexdigest(),status='EXPLORATORY_DIAGNOSTIC_AFTER_PRIMARY_RESULTS',reason='Inherited epithelial lineage includes1379cells with explicitly flagged Mixed_Macrophage_cluster_review subtype; late epithelial contrast requires sensitivity to flagged mixture and subtype composition.',new_annotations=False,exclusion='Exclude exactly the inherited Mixed_Macrophage_cluster_review flag; do not infer new identities.',composition='AT2/Club/Multiciliated inherited labels have>=30cells in each of11source units based on metadata, before these subset scores; observed fractions versus fixed equal-animal mean uninfected fractions within these3types.',reference='Reuse original mouse gene mean/SD and original signs; no new reference fitting.',tests='9exploratory Welch tests:3methods x total/IFN overlap/outside IFN; BH9 including all9; source-level exact tests also retained.',scope='Inherited animal annotations have not undergone human-specific Gate2B review; this is a within-resource sensitivity, not new subtype validation.')
dest=OUT/'parameters/mouse_identity_diagnostic_amendment.json'
if not dest.exists():dump(dest,amendment)
c=pd.read_parquet(RESOURCE/'CRRM_revision_work/Result4/data/atlas_cell_metadata.parquet',filters=[('source_label','==','RelA_mouse')]);c=c[c.identity_level2.eq('Epithelial_cell')].copy()
c['subtype']=c.gate2b_fine_subtype.astype(str);c['unit']=c.verified_unit_id.astype(str)
c['group']=c.infection_status.astype(str).map({'corn oil':'Uninfected vehicle','corn oil-RSV':'RSV vehicle','TMX-RSV':'RSV RelA CKO'})
codes,unique=pd.factorize(pd.MultiIndex.from_frame(c[['unit','group','subtype']]),sort=True)
meta=pd.DataFrame(list(unique),columns=['unit','group','subtype']);meta['cells']=np.bincount(codes,minlength=len(meta));c['subprofile']=codes
var=pd.read_parquet(RESOURCE/'DATA/curated_metadata/features/RelA_mouse_var.parquet');projection=pd.read_csv(OUT/'parameters/frozen_mouse_member_projection.csv');index={g:i for i,g in enumerate(var.native_gene_id)}
indices=np.array([index.get(g,-1) for g in projection.ensembl]);pos=np.flatnonzero(indices>=0);cols=indices[pos]
rowlookup=np.full(84967,-1,dtype=int);rowlookup[c.expression_row_index.to_numpy()]=codes
sums=np.zeros((len(meta),388));counts=np.zeros(len(meta),int)
with h5py.File(RESOURCE/'DATA/expression/RelA_mouse.h5ad','r') as h:
    ptr=h['X/indptr'][:];n,p=map(int,h['X'].attrs['shape']);host=var.host_gene_flag.to_numpy(dtype=bool)
    for a in range(0,n,8000):
        b=min(a+8000,n);lo,hi=int(ptr[a]),int(ptr[b]);ids=rowlookup[a:b];ok=ids>=0
        if not ok.any():continue
        x=sparse.csr_matrix((h['X/data'][lo:hi],h['X/indices'][lo:hi],ptr[a:b+1]-lo),shape=(b-a,p));xx=x[ok];depth=np.asarray(xx[:,host].sum(1)).ravel().astype(float)
        expr=np.log1p(xx[:,cols].toarray()/depth[:,None]*10000)
        agg=sparse.csr_matrix((np.ones(ok.sum()),(ids[ok],np.arange(ok.sum()))),shape=(len(meta),ok.sum()));sums[:,pos]+=agg@expr;counts+=np.bincount(ids[ok],minlength=len(meta))
assert np.array_equal(counts,meta.cells)
means=sums/counts[:,None]
ref=pd.read_csv(OUT/'data/gene_reference_parameters.csv');ref=ref[ref.source.eq('RelA_mouse')].set_index('human_ensembl').loc[projection.human_ensembl]
sign=np.where(projection.frozen_direction.eq('up'),1.,-1.);valid=ref.variable.to_numpy(dtype=bool)
z=np.full_like(means,np.nan);z[:,valid]=(means[:,valid]-ref.reference_mean.to_numpy()[valid])/ref.reference_SD.to_numpy()[valid]*sign[valid]
sets=json.loads((OUT/'parameters/score_gene_sets.json').read_text());targets=['CRRM388','IFN-list overlap','Outside IFN lists']
for t in targets:meta[t]=z[:,valid&projection.human_ensembl.isin(sets[t])].mean(1)
meta.to_csv(OUT/'results/mouse_epithelial_subtype_scores.csv',index=False)
types=['Alveolar_AT2','Club_cell','Multiciliated_cell'];sub=meta[meta.subtype.isin(types)]
assert (sub.cells>=30).all() and len(sub)==33
w=sub.pivot(index='unit',columns='subtype',values='cells');w=w.div(w.sum(1),axis=0)
healthy=meta[meta.group.eq('Uninfected vehicle')].unit.unique();fixed=w.loc[healthy].mean()
rows=[]
for unit,g in meta.groupby('unit',sort=False):
    nonmixed=g[~g.subtype.eq('Mixed_Macrophage_cluster_review')]
    gg=g[g.subtype.isin(types)].set_index('subtype').loc[types]
    vals={'Exclude flagged mixed cluster':np.average(nonmixed[targets].to_numpy(),weights=nonmixed.cells,axis=0),'Three-type observed composition':np.average(gg[targets].to_numpy(),weights=gg.cells,axis=0),'Three-type fixed healthy composition':fixed.loc[types].to_numpy()@gg[targets].to_numpy()}
    for method,value in vals.items():
        rows.append(dict(unit=unit,group=g.group.iloc[0],method=method,**dict(zip(targets,value))))
scores=pd.DataFrame(rows);scores.to_csv(OUT/'results/mouse_epithelial_composition_scores.csv',index=False)
def test(a,b):
    n,m=len(a),len(b);va=np.var(a,ddof=1)/n;vb=np.var(b,ddof=1)/m;df=(va+vb)**2/(va*va/(n-1)+vb*vb/(m-1));se=np.sqrt(va+vb);d=np.mean(a)-np.mean(b);q=stats.t.ppf(.975,df);return dict(estimate=d,lower=d-q*se,upper=d+q*se,P=2*stats.t.sf(abs(d/se),df),n_case=n,n_control=m)
effects=[]
for method,g in scores.groupby('method',sort=False):
    for t in targets:effects.append(dict(method=method,score=t,**test(g.loc[g.group.eq('RSV RelA CKO'),t].to_numpy(),g.loc[g.group.eq('RSV vehicle'),t].to_numpy())))
ef=pd.DataFrame(effects);p=ef.P.to_numpy();ix=np.argsort(p);v=np.minimum.accumulate((p[ix]*9/np.arange(1,10))[::-1])[::-1];q=np.empty(9);q[ix]=np.minimum(v,1);ef['adjusted_P_BH9']=q
ef.to_csv(OUT/'results/mouse_identity_composition_effects.csv',index=False)
# Reaggregate all subtype means and check against the original annotated-lineage score.
original=pd.read_csv(OUT/'results/profile_scores.csv');original=original[original.source.eq('RelA_mouse')&original.method.eq('Primary')&original.celltype.eq('Epithelial_cell')].set_index('unit')
errors=[]
for unit,g in meta.groupby('unit'):
    errors.extend(np.average(g[targets],weights=g.cells,axis=0)-original.loc[unit,targets].to_numpy(dtype=float))
assert np.max(np.abs(errors))<1e-12
dump(OUT/'audit/mouse_identity_diagnostic_receipt.json',dict(status='PASS',excluded_flagged_cells=int(c.subtype.eq('Mixed_Macrophage_cluster_review').sum()),all_epithelial_cells=len(c),fixed_composition=fixed.to_dict(),retained_types=types,animals=11,primary_reaggregation_max_error=float(np.max(np.abs(errors))),diagnostic_tests=9))
print(ef.to_string(index=False))
