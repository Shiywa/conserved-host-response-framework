"""Reuse descriptor counts/annotations; aggregate once at source animal x lineage."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib, json, time
import h5py
import numpy as np
import pandas as pd
from scipy import sparse

OUT=(Path(__import__('os').environ['CRRM_WORKSPACE']) / 'Result6_analysis')
RESOURCE=Path(str(Path(__import__('os').environ.get('CRRM_RESOURCE_ROOT','external_resource')) / ''))
def dump(path,value): path.write_text(json.dumps(value,indent=2)+'\n')
plan=json.loads((OUT/'parameters/analysis_plan_v1.json').read_text())
assert hashlib.sha256((OUT/'parameters/analysis_plan_v1.json').read_bytes()).hexdigest()==json.loads((OUT/'parameters/declaration_receipt.json').read_text())['plan_sha256']
meta=pd.read_csv(OUT/'data/mouse_profile_metadata.csv')
members=pd.read_parquet(OUT/'data/mouse_row_membership.parquet').sort_values('expression_row_index')
var=pd.read_parquet(RESOURCE/'DATA/curated_metadata/features/RelA_mouse_var.parquet')
mapping=pd.read_csv(OUT/'parameters/frozen_mouse_member_projection.csv')
core_ids=mapping.human_ensembl.to_numpy(dtype=str)
index={g:i for i,g in enumerate(var.native_gene_id)}
feature_indices=np.array([index.get(g,-1) for g in mapping.ensembl],dtype=int)
with np.load(RESOURCE/'DATA/curated_metadata/features/gene_availability/RelA_mouse.npz',allow_pickle=False) as z:
    assert np.array_equal(z['genes'].astype(str),var.feature_key.astype(str).to_numpy())
    av=z['availability'].all(axis=0)
measured=np.array([j>=0 and av[j] for j in feature_indices])
pos=np.flatnonzero(measured);selected=feature_indices[pos]
def hcol(h,name):
    g=h['obs'][name]
    if isinstance(g,h5py.Dataset): return g.asstr()[:] if h5py.check_string_dtype(g.dtype) else g[:]
    cats=g['categories'].asstr()[:];codes=g['codes'][:];return cats[codes]
start=time.monotonic()
with h5py.File(RESOURCE/'DATA/expression/RelA_mouse.h5ad','r') as h:
    n,p=map(int,h['X'].attrs['shape'])
    assert (n,p)==(len(members),len(var))
    assert np.array_equal(hcol(h,'cell_uid'),members.cell_uid.to_numpy())
    indptr=h['X/indptr'][:]
    sums=np.zeros((len(meta),p));raw=np.zeros((len(meta),388));det=np.zeros_like(raw)
    counts=np.zeros(len(meta),int);host_counts=np.zeros(len(meta));log_depth=np.zeros(len(meta))
    memberships=members.profile_id.to_numpy();host=var.host_gene_flag.to_numpy(dtype=bool)
    for a in range(0,n,8000):
        b=min(a+8000,n);lo,hi=int(indptr[a]),int(indptr[b])
        x=sparse.csr_matrix((h['X/data'][lo:hi].astype(float),h['X/indices'][lo:hi],indptr[a:b+1]-lo),shape=(b-a,p))
        total=np.asarray(x[:,host].sum(1)).ravel()
        assert (total>0).all() and (x.data>=0).all() and np.isfinite(x.data).all()
        ids=memberships[a:b]
        agg=sparse.csr_matrix((np.ones(b-a),(ids,np.arange(b-a))),shape=(len(meta),b-a))
        xx=x.multiply((10000/total)[:,None]).tocsr();np.log1p(xx.data,out=xx.data)
        sums+=(agg@xx).toarray()
        sel=x[:,selected];raw[:,pos]+=(agg@sel).toarray()
        detected=sel.copy();detected.data[:]=1;det[:,pos]+=(agg@detected).toarray()
        counts+=np.bincount(ids,minlength=len(meta));host_counts+=np.bincount(ids,weights=total,minlength=len(meta));log_depth+=np.bincount(ids,weights=np.log1p(total),minlength=len(meta))
        if a==0 or b==n: print('Mouse profiles:',b,'/',n,'cells;',round(time.monotonic()-start,1),'s',flush=True)
    assert np.array_equal(counts,meta.cells.to_numpy())
    means=sums/counts[:,None];cx=np.full((len(meta),388),np.nan);cx[:,pos]=means[:,selected]
    np.savez_compressed(OUT/'data/mouse_expression_profiles.npz',profile_ids=meta.profile_id.to_numpy(),mean_log_expression=cx,transcriptome_mean_log_expression=means,native_gene_ids=var.native_gene_id.to_numpy(dtype=str),human_ensembl=core_ids,common_measured=measured,sum_native_counts=raw,expressing_cells=det,cells=counts,host_counts=host_counts,mean_log_host_depth=log_depth/counts)
receipt=dict(status='PASS',completed_at=datetime.now(timezone.utc).isoformat(),cells=n,profile_rows=len(meta),native_features=p,mapped_common_CRRM=int(measured.sum()),wall_seconds=time.monotonic()-start,raw_counts_modified=False,annotations_recomputed=False)
dump(OUT/'audit/mouse_extraction_receipt.json',receipt)
print(json.dumps(receipt,indent=2))
