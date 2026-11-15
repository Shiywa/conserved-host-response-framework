from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,shutil
import h5py,numpy as np,pandas as pd

OUT=(Path(__import__('os').environ['CRRM_WORKSPACE']) / 'Result4_analysis')
RESOURCE=Path(str(Path(__import__('os').environ.get('CRRM_RESOURCE_ROOT','external_resource')) / ''))
META=RESOURCE/'DATA/curated_metadata'
REV=Path(str(Path(__import__('os').environ['CRRM_WORKSPACE']) / ''))
SHA=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def dump(p,x):p.write_text(json.dumps(x,ensure_ascii=False,indent=2,default=lambda x:x.item() if isinstance(x,np.generic) else str(x))+'\n')
def hcol(h,key):
    x=h['obs'][key]
    if isinstance(x,h5py.Dataset):return x.asstr()[:] if h5py.check_string_dtype(x.dtype) else x[:]
    cats=x['categories'].asstr()[:] if h5py.check_string_dtype(x['categories'].dtype) else x['categories'][:]
    codes=x['codes'][:];a=np.empty(len(codes),dtype=object);a[:]=None;a[codes>=0]=cats[codes[codes>=0]];return a

core_path=REV/'phase18_Result2_external_response_20261112/parameters/frozen_CRRM388.csv'
branch_path=core_path.parent/'frozen_functional_branches.csv'
core=pd.read_csv(core_path);branches=pd.read_csv(branch_path)
shutil.copy2(core_path,OUT/'parameters/frozen_CRRM388.csv');shutil.copy2(branch_path,OUT/'parameters/frozen_functional_branches.csv')
assert len(core)==388 and core.human_ensembl.is_unique
print('Reading frozen cell index and annotations',flush=True)
cols=['cell_uid','original_cell_id','source_id','study_id','source_label','source_study_alias','sample_id','library_id','donor_id','verified_unit_id','unit_mapping_status','unit_kind','species_namespace','species','tissue','experimental_context','infection_status','time','qc_group_id','assay_record_type','author_overview_scope_included','author_intent_sensitivity_included','author_scope_status','expression_row_index']
cells=pd.read_parquet(META/'metadata/cell_index.parquet',columns=cols)
ann_cols=['cell_uid','current_identity','identity_level1','identity_level2','identity_level3','original_cell_type','original_ann_level_1','gate2b_canonical_identity','gate2b_fine_subtype','gate2b_identity_status','gate2b_review_status','gate2b_evidence_grade']
annotations=pd.read_parquet(META/'annotation/cell_annotation.parquet',columns=ann_cols)
assert cells.cell_uid.is_unique and annotations.cell_uid.is_unique and len(cells)==len(annotations)==2474742
cells=cells.merge(annotations,on='cell_uid',how='left',validate='one_to_one',sort=False)
unknown={'','NULL','Unknown','unknown','UNRESOLVED'}
cells['unit_resolved']=cells.verified_unit_id.notna()&~cells.verified_unit_id.astype(str).isin(unknown)&~cells.unit_mapping_status.astype(str).isin(['UNRESOLVED','UNKNOWN','NULL',''])
for c in cells:
    if c not in ['cell_uid','original_cell_id','expression_row_index'] and cells[c].dtype=='object':cells[c]=cells[c].astype('category')
cells.to_parquet(OUT/'data/atlas_cell_metadata.parquet',index=False)
summaries=[];coverage=[];checks=[]
for label,g in cells.groupby('source_label',observed=True,sort=True):
    p=RESOURCE/'DATA/expression'/(str(label)+'.h5ad');v=pd.read_parquet(META/'features'/(str(label)+'_var.parquet'))
    gg=g.sort_values('expression_row_index');row_indices=gg.expression_row_index.to_numpy()
    assert np.array_equal(row_indices,np.arange(len(g))),label
    with h5py.File(p,'r') as h:
        shape=tuple(int(x) for x in h['X'].attrs['shape']);assert shape==(len(g),len(v)),(label,shape)
        uid=hcol(h,'cell_uid');assert np.array_equal(uid,gg.cell_uid.to_numpy()),label
        data=h['X/data'][:10000];assert np.isfinite(data).all() and (data>=0).all(),label
        summary=dict(source_label=str(label),species=str(g.species_namespace.iloc[0]),cells=len(g),cells_with_resolved_units=int(g.unit_resolved.sum()),recorded_resolved_units=int(g.loc[g.unit_resolved,'verified_unit_id'].nunique()),sample_labels=int(g.sample_id.nunique()),library_labels=int(g.library_id.nunique()),native_features=len(v),matrix_shape=str(shape),matrix_encoding=str(h['X'].attrs['encoding-type']),matrix_dtype=str(h['X/data'].dtype),local_expression_bytes=p.stat().st_size,expression_path=str(p),rna_record_types=';'.join(sorted(g.assay_record_type.astype(str).unique())),contexts=';'.join(sorted(g.experimental_context.astype(str).unique())),conditions=';'.join(sorted(g.infection_status.astype(str).unique())))
    masks=META/'features/gene_availability'/(str(label)+'.npz')
    if not masks.exists():masks=RESOURCE/'DATA/features/gene_availability'/(str(label)+'.npz')
    summary['measurement_mask_path']=str(masks)
    with np.load(masks,allow_pickle=False) as m:
        assert np.array_equal(m['genes'].astype(str),v.feature_key.astype(str).to_numpy()),label
        av=m['availability'];assert av.shape[1]==len(v)
        if summary['species']=='human':
            index={x:i for i,x in enumerate(v.native_gene_id.astype(str).str.replace(r'\.\d+$','',regex=True))}
            count=[]
            for r in core.itertuples():
                j=index.get(r.human_ensembl);vals=av[:,j] if j is not None else np.zeros(av.shape[0],dtype=bool)
                coverage.append(dict(source_label=str(label),human_ensembl=r.human_ensembl,display_symbol=r.display_symbol,frozen_direction=r.frozen_direction,encoded_in_native_features=j is not None,measured_in_every_library=bool(vals.all()),measured_in_any_library=bool(vals.any()),libraries_measured=int(vals.sum()),libraries_total=len(vals)))
                count.append(bool(vals.all()))
            present=core[np.array(count)]
            summary.update(CRRM388_encoded=int(core.human_ensembl.isin(index).sum()),CRRM388_common_measured=int(sum(count)),up_common_measured=int((present.frozen_direction=='up').sum()),down_common_measured=int((present.frozen_direction=='down').sum()),CRRM_projection='direct stable human Ensembl IDs')
        else:summary.update(CRRM_projection='native species retained; frozen cross-species mapping required before scoring')
    summaries.append(summary);checks.append(dict(source_label=str(label),cell_uid_order='PASS',native_dimensions='PASS',nonnegative_native_X_sample='PASS'))
    print(str(label),len(g),'cells; known units',summary['recorded_resolved_units'],'CRRM common',summary.get('CRRM388_common_measured','native species'),flush=True)

pd.DataFrame(summaries).to_csv(OUT/'audit/atlas_source_summary.csv',index=False)
pd.DataFrame(coverage).to_csv(OUT/'audit/CRRM388_source_measurement.csv',index=False)
for names,filename in [(['source_label','infection_status','unit_mapping_status'],'atlas_condition_unit_counts.csv'),(['source_label','identity_level1','identity_level2'],'atlas_identity_counts.csv'),(['source_label','verified_unit_id','unit_mapping_status','identity_level2','infection_status','time'],'atlas_unit_identity_condition_counts.csv')]:
    cells.groupby(names,dropna=False,observed=True).size().rename('cells').reset_index().to_csv(OUT/'audit'/filename,index=False)
for label in ['Berg_ALI','Infant_PBMC','CD1c_DC','iAT2_macrophage']:
    g=cells[cells.source_label==label]
    g.groupby(['original_cell_type','identity_level2','gate2b_canonical_identity'],dropna=False,observed=True).size().rename('cells').reset_index().to_csv(OUT/'audit'/(label+'_identity_crosswalk.csv'),index=False)
comparisons=pd.read_csv(META/'comparison/comparison_eligibility.tsv',sep='\t',dtype=str,keep_default_na=False)
keep=['comparison_id','study_label','comparison_kind','cell_identity','case_condition','control_condition','independent_units_case','independent_units_control','qualified_pair_count','case_cells','control_cells','pairing_status','sampling_time','eligibility','effective_decision','confounding_flags','age_sex_full_rank','condition_identifiable_from_recorded_batch']
comparisons[keep].to_csv(OUT/'audit/inherited_comparison_catalog.csv',index=False)
comparisons.groupby(['study_label','comparison_kind','eligibility'],dropna=False).size().rename('comparisons').reset_index().to_csv(OUT/'audit/source_comparison_scope.csv',index=False)
inputs=[core_path,branch_path,META/'metadata/cell_index.parquet',META/'annotation/cell_annotation.parquet',META/'comparison/comparison_eligibility.tsv',META/'metadata/sample_metadata.tsv']
receipt=dict(status='PASS',completed_at=datetime.now(timezone.utc).isoformat(),all_frozen_index_records=len(cells),sources=len(summaries),species_counts=cells.species_namespace.value_counts().to_dict(),identity_counts=cells.identity_level1.value_counts().to_dict(),local_H5AD_records=sum(s['cells'] for s in summaries),source_order_checks=checks,CRRM388_genes_unchanged=True,CRRM388_scores_computed=False,annotation_recomputed=False,descriptor_project_edited=False,inputs=[dict(path=str(p),sha256=SHA(p)) for p in inputs])
dump(OUT/'audit/atlas_audit_receipt.json',receipt)
print(json.dumps({k:v for k,v in receipt.items() if k not in ['source_order_checks','inputs']},ensure_ascii=False,indent=2),flush=True)
