from pathlib import Path
from datetime import datetime,timezone
import json,time
import h5py,numpy as np,pandas as pd
from scipy import sparse

OUT=(Path(__import__('os').environ['CRRM_WORKSPACE']) / 'Result4_analysis');RESOURCE=Path(str(Path(__import__('os').environ.get('CRRM_RESOURCE_ROOT','external_resource')) / ''));META=RESOURCE/'DATA/curated_metadata'
plan=json.loads((OUT/'parameters/analysis_plan_v1.json').read_text());assert plan['status']=='FROZEN_BEFORE_NEW_CRRM388_SCORING'
core=pd.read_csv(OUT/'parameters/frozen_CRRM388.csv');groups=pd.read_csv(OUT/'data/profile_metadata.csv');members=pd.read_parquet(OUT/'data/profile_row_membership.parquet')
receipts=[]
for label in sorted(members.source_label.unique()):
    start=time.monotonic();dest=OUT/'data'/('native_profiles_'+label+'.npz');gp=groups[groups.source_label==label].sort_values('profile_id');groupids=gp.profile_id.to_numpy();global_to_local={g:i for i,g in enumerate(groupids)}
    mp=members[members.source_label==label];var=pd.read_parquet(META/'features'/(label+'_var.parquet'));index={x:i for i,x in enumerate(var.native_gene_id.astype(str).str.replace(r'\.\d+$','',regex=True))}
    feature_indices=np.array([index.get(g,-1) for g in core.human_ensembl],dtype=np.int64);mask_path=META/'features/gene_availability'/(label+'.npz')
    if not mask_path.exists():mask_path=RESOURCE/'DATA/features/gene_availability'/(label+'.npz')
    with np.load(mask_path,allow_pickle=False) as z:
        av=z['availability'].all(axis=0);measured=np.array([i>=0 and av[i] for i in feature_indices])
    positions=np.flatnonzero(measured);target_cols=feature_indices[positions];host=var.host_gene_flag.to_numpy(dtype=bool)
    path=RESOURCE/'DATA/expression'/(label+'.h5ad')
    if dest.exists():
        old=json.loads(dest.with_suffix('.json').read_text());assert old['status']=='PASS' and old['plan_declared_at']==plan['declared_at'];receipts.append(old);print(label,'existing verified extraction',flush=True);continue
    with h5py.File(path,'r') as h:
        n,p=map(int,h['X'].attrs['shape']);indptr=h['X/indptr'][:]
        memberships=[]
        for view,vs in mp.groupby('view',sort=False):
            ids=np.full(n,-1,dtype=np.int32);ids[vs.expression_row_index.to_numpy()]=[global_to_local[g] for g in vs.profile_id];memberships.append(ids)
        sums=np.zeros((len(gp),388));raw_sums=np.zeros_like(sums);detected=np.zeros_like(sums);counts=np.zeros(len(gp),dtype=np.int64);host_sums=np.zeros(len(gp));log_host_sums=np.zeros(len(gp));zero_library=0
        for a in range(0,n,10000):
            b=min(a+10000,n);lo,hi=int(indptr[a]),int(indptr[b]);x=sparse.csr_matrix((h['X/data'][lo:hi],h['X/indices'][lo:hi],indptr[a:b+1]-lo),shape=(b-a,p));total=np.asarray(x.sum(axis=1) if host.all() else x[:,host].sum(axis=1)).ravel().astype(float);valid=total>0;zero_library+=int((~valid).sum())
            selected=x[:,target_cols].toarray().astype(float);log=np.zeros_like(selected);np.divide(10000*selected,total[:,None],out=log,where=valid[:,None]);np.log1p(log,out=log)
            for mem in memberships:
                ids=mem[a:b];ok=(ids>=0)&valid;rows=np.flatnonzero(ok)
                if not len(rows):continue
                agg=sparse.csr_matrix((np.ones(len(rows)),(ids[rows],rows)),shape=(len(gp),b-a))
                sums[:,positions]+=agg@log;raw_sums[:,positions]+=agg@selected;detected[:,positions]+=agg@(selected>0).astype(float)
                counts+=np.bincount(ids[rows],minlength=len(gp));host_sums+=np.bincount(ids[rows],weights=total[rows],minlength=len(gp));log_host_sums+=np.bincount(ids[rows],weights=np.log1p(total[rows]),minlength=len(gp))
            if a==0 or b==n or (b%100000==0):print(label,b,'/',n,'cells;',round(time.monotonic()-start,1),'seconds',flush=True)
        assert zero_library==0,(label,zero_library)
        assert np.array_equal(counts,gp.cells.to_numpy()),label
        means=sums/counts[:,None];means[:,~measured]=np.nan
        np.savez_compressed(dest,profile_ids=groupids,mean_log_expression=means,sum_native_counts=raw_sums,expressing_cells=detected,cells=counts,host_counts=host_sums,mean_log_host_depth=log_host_sums/counts,common_measured=measured,human_ensembl=core.human_ensembl.to_numpy(dtype=str))
        receipt=dict(status='PASS',source=label,physical_cells=n,profile_rows=len(gp),mean_expression_dimensions=list(means.shape),common_measured_core_genes=int(measured.sum()),zero_host_library_cells=zero_library,plan_declared_at=plan['declared_at'],completed_at=datetime.now(timezone.utc).isoformat(),wall_seconds=time.monotonic()-start)
        dest.with_suffix('.json').write_text(json.dumps(receipt,indent=2)+'\n');receipts.append(receipt)
(OUT/'audit/expression_extraction_receipt.json').write_text(json.dumps(dict(status='PASS',human_sources=len(receipts),physical_human_cells=sum(r['physical_cells'] for r in receipts),sources=receipts,original_counts_rewritten=False),indent=2)+'\n')
print('COMPLETE:',sum(r['physical_cells'] for r in receipts),'physical human cells',flush=True)
