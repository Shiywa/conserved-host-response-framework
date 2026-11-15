from pathlib import Path
from datetime import datetime,timezone
import hashlib,json
import numpy as np,pandas as pd

OUT=(Path(__import__('os').environ['CRRM_WORKSPACE']) / 'Result4_analysis');RESOURCE=Path(str(Path(__import__('os').environ.get('CRRM_RESOURCE_ROOT','external_resource')) / ''));META=RESOURCE/'DATA/curated_metadata';REV=Path(str(Path(__import__('os').environ['CRRM_WORKSPACE']) / ''))
SHA=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def dump(p,x):p.write_text(json.dumps(x,ensure_ascii=False,indent=2,default=lambda x:x.item() if isinstance(x,np.generic) else str(x))+'\n')
assert json.loads((OUT/'audit/atlas_audit_receipt.json').read_text())['status']=='PASS'
cells=pd.read_parquet(OUT/'data/atlas_cell_metadata.parquet');human=cells[cells.species_namespace=='human'].copy()
assert len(human)==2119831
for c in human:
    if isinstance(human[c].dtype,pd.CategoricalDtype):human[c]=human[c].astype(object)
human['analysis_unit']=np.where(human.unit_resolved,human.verified_unit_id.astype(str),'recorded_library::'+human.source_label.astype(str)+'::'+human.library_id.astype(str))
human['analysis_celltype']=human.identity_level2.fillna('Unknown').astype(str);human['view']='canonical'
for c in ['infection_status','time','study_id','analysis_unit']:human[c]=human[c].fillna('Unknown').astype(str)
b=human[human.source_label=='Berg_ALI'].copy();author=pd.read_csv(OUT/'data/Berg_author_cell_metadata.csv',dtype={'donor':str});assert author.original_cell_id.is_unique
b['author_full_barcode']=b.original_cell_id.astype(str).str.rsplit(':',n=1).str[-1]
b=b.merge(author[['original_cell_id','donor','Infection','time','res.0.5','nCount_RNA','nCount_VIRUS','percVirus','Infection_load','percent.mt']],left_on='author_full_barcode',right_on='original_cell_id',how='left',validate='many_to_one',suffixes=('','_author'))
matched=b['donor'].notna()&b.unit_resolved
assert matched.sum()==b.unit_resolved.sum()==10759
assert (b.loc[matched,'verified_unit_id']=='Berg_ALI::'+b.loc[matched,'donor']).all()
assert (b.loc[matched,'infection_status']==b.loc[matched,'Infection']).all()
assert (b.loc[matched,'time']==b.loc[matched,'time_author']).all()
payload=b.loc[matched,['cell_uid','author_full_barcode','donor','Infection','time','res.0.5','nCount_RNA','nCount_VIRUS','percVirus','Infection_load','percent.mt']].copy();payload.to_parquet(OUT/'data/Berg_author_metadata_exact_crosswalk.parquet',index=False)
fine=human[human.cell_uid.isin(payload.cell_uid)].merge(payload[['cell_uid','res.0.5']],on='cell_uid',validate='one_to_one');fine['analysis_celltype']=fine['res.0.5'];fine['view']='author_epithelial';fine=fine.drop(columns=['res.0.5'])
profiles=pd.concat([human,fine],ignore_index=True)
keys=['source_label','study_id','analysis_unit','unit_resolved','analysis_celltype','infection_status','time','view']
keyindex=pd.MultiIndex.from_frame(profiles[keys]);codes,uniques=pd.factorize(keyindex,sort=True)
profiles['profile_id']=codes
groups=pd.DataFrame(list(uniques),columns=keys);groups['profile_id']=np.arange(len(groups));counts=np.bincount(codes,minlength=len(groups));groups['cells']=counts
groups.to_csv(OUT/'data/profile_metadata.csv',index=False)
profiles[['cell_uid','source_label','expression_row_index','profile_id','view']].to_parquet(OUT/'data/profile_row_membership.parquet',index=False)

geo=pd.read_csv(OUT/'data/Infant_PBMC_public_GEO_covariates.csv');geo['participant_id']=geo.title.str.extract(r'(p(?:RSV|HC|CoV)\d+)')[0]
geo['age_months']=pd.to_numeric(geo['age (months):ch1'],errors='coerce');geo['sex']=geo['Sex:ch1'].str.lower().map({'female':'F','male':'M'})
assert geo.participant_id.dropna().is_unique
clinical=geo[geo.participant_id.str.startswith(('pRSV','pHC'),na=False)].copy();assert len(clinical)==36 and clinical[['age_months','sex']].notna().all().all()
clinical['analysis_unit']='Infant_PBMC::'+clinical.participant_id
assert set(clinical.analysis_unit)==set(groups.loc[(groups.source_label=='Infant_PBMC')&groups.unit_resolved,'analysis_unit'])
clinical[['analysis_unit','geo_accession','participant_id','age_months','sex']].to_csv(OUT/'data/patient_covariates_internal.csv',index=False)

core=pd.read_csv(OUT/'parameters/frozen_CRRM388.csv');branches=pd.read_csv(OUT/'parameters/frozen_functional_branches.csv')
s1=pd.read_csv(REV/'phase17_key_supplement_20261112/tables/S1A_CRRM_genes.csv');outside=set(s1.loc[s1['IFN reference overlap']=='No','Ensembl gene ID']);overlap=set(s1.loc[s1['IFN reference overlap']=='Yes','Ensembl gene ID']);assert len(outside)==309 and len(overlap)==79
gene_sets={'CRRM388':core.human_ensembl.tolist(),'Up':core.loc[core.frozen_direction=='up','human_ensembl'].tolist(),'Down':core.loc[core.frozen_direction=='down','human_ensembl'].tolist(),'IFN-list overlap':sorted(overlap),'Outside IFN lists':sorted(outside)}
for name,g in branches.groupby('branch',sort=False):gene_sets[name]=g.human_ensembl.tolist()
gene_sets['Degranulation without LCN2']=[g for g in gene_sets['Neutrophil degranulation'] if g not in set(core.loc[core.display_symbol=='LCN2','human_ensembl'])]
assert len(gene_sets['Degranulation without LCN2'])==24
dump(OUT/'parameters/score_gene_sets.json',gene_sets)
comparisons=[]
for cell in ['B_cell','T_cell','NK_cell','Monocyte']:comparisons.append(dict(id='PBMC_'+cell,source='Infant_PBMC',view='canonical',celltype=cell,case=['Mild RSV','Moderate RSV','Severe RSV'],control=['Healthy'],design='unpaired_participants',primary=True))
for cell,view in [('Epithelial_cell','canonical'),('Ciliated','author_epithelial')]:comparisons.append(dict(id='Berg_'+cell,source='Berg_ALI',view=view,celltype=cell,case=['inf'],control=['ctrl'],design='within_donor_time_balanced',primary=True))
for cell in ['Dendritic_cell','Monocyte']:comparisons.append(dict(id='DC_'+cell,source='CD1c_DC',view='canonical',celltype=cell,case=['RSV'],control=['Untreated'],design='paired_donors',primary=True))
assert len(comparisons)==8
planned=[]
for x in comparisons:
    g=groups[(groups.source_label==x['source'])&(groups.view==x['view'])&(groups.analysis_celltype==x['celltype'])&groups.unit_resolved&(groups.cells>=30)]
    cases=g[g.infection_status.isin(x['case'])];controls=g[g.infection_status.isin(x['control'])];common=set(cases.analysis_unit)&set(controls.analysis_unit)
    planned.append(dict(**x,case_units_metadata=int(cases.analysis_unit.nunique()),control_units_metadata=int(controls.analysis_unit.nunique()),paired_units_metadata=len(common)))
pd.DataFrame(planned).to_csv(OUT/'parameters/primary_comparison_design.csv',index=False)

plan=dict(status='FROZEN_BEFORE_NEW_CRRM388_SCORING',declared_at=datetime.now(timezone.utc).isoformat(),scope='Use the Data Descriptor as the source backbone; all 15 sources audited; source-relative human expression localization across 13 source labels, followed by donor/participant-level within-source inference',exploratory_status='Retrospective cellular analysis after bulk results and earlier historical 70-gene analyses were known; no new score outcomes used to select this protocol',all_resource_cells=2474742,human_expression_cells=2119831,human_source_labels=13,animal_context='Native namespaces and cells retained in audit; species projection required before cross-species CRRM scoring; animal effects not pooled with human inference',annotation='Reuse frozen project canonical lineages; Berg epithelial subtypes imported only through exact full barcode and concordant donor/condition/time; candidate identities retain their evidence status',scoring=dict(input='native source X counts; each cell log1p(10000 * host_gene_count / total_measured_host_counts)',aggregation='mean normalized log expression per source study, biological unit, cell identity, condition and time; unresolved-unit observations retained for descriptive profiles only',features='same frozen388 IDs and signs; common measured features within source; no missing-feature imputation',reference='source study-relative gene means/SD across eligible canonical profiles, with equal total weight for each biological unit; exclude unresolved units in three inferential sources and exclude COVID/Other from PBMC inference scaling',numerical_SD_guard=1e-8,total_eligibility='>=70% of388 and >=50% of each frozen direction; branches >=5 genes and >=50% nominal membership',direction='equal per-gene weights with the frozen up/down sign; positive signed change indicates alignment with the experimental response, not activation of every annotated pathway'),minimum_cells_per_profile=30,minimum_inferential_units=3,primary_comparisons=planned,error_control=dict(total='BH8 across eight declared total-score comparisons; unavailable comparisons retained as P=1 in the family',branches='BH56 across seven fixed overlapping branches and eight declared comparisons',additional_components='separate BH16 for IFN-list-overlap and outside-list scores across the eight comparisons',secondary='Berg bystander/control and infected/bystander; PBMC Severe/Mild; covariate-adjusted PBMC estimates and LCN2 omission each retain separate declared reporting families; no choosing minimum P across methods'),inference=dict(primary='paired t for biological-unit contrasts and Welch t for independent participants; pointwise95% t intervals',Berg='first subtract case/control within donor and the same time point; equally average eligible matched-time contrasts within each donor; one inferential observation per donor',sensitivities='exact donor sign flips; patient-label resampling; minimum-cell thresholds15/50; complete common-time comparisons; author-intent population where mapped; no cell-level biological replication'),composition=dict(population='36 Healthy/RSV participants; four common B/T/NK/Monocyte lineages with adequate observed cells in every patient',observed='sample mean score using observed fractions within the four-lineage compartment',standardized='same cell-lineage means weighted by the equal-patient mean healthy lineage composition',uncertainty='stratified patient bootstrap; no causal mediation or exact whole-blood attribution'),source_overlap='Berg/Infant_PBMC legacy copies are aliases of the same descriptor source, never additional validation studies; HLCA components retained by original study_id; healthy/history/memory/post-exposure sources keep their original roles',CRRM388_scores_computed=False,metadata_audit_receipt_sha256=SHA(OUT/'audit/atlas_audit_receipt.json'),frozen_genes_sha256=SHA(OUT/'parameters/frozen_CRRM388.csv'),functional_branches_sha256=SHA(OUT/'parameters/frozen_functional_branches.csv'))
dump(OUT/'parameters/analysis_plan_v1.json',plan)
dump(OUT/'parameters/declaration_receipt.json',dict(declared_at=plan['declared_at'],analysis_plan_sha256=SHA(OUT/'parameters/analysis_plan_v1.json'),scores_computed=False,raw_author_barcode_matches=len(payload),physical_human_cells=len(human),profiles=len(groups)))
print(json.dumps({'status':plan['status'],'human_cells':len(human),'source_labels':13,'profiles':len(groups),'exact_Berg_author_matches':len(payload),'primary_comparisons':len(comparisons)},ensure_ascii=False,indent=2))
