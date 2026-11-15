"""Fix sources, contrasts and correction families before new CRRM scoring."""
from pathlib import Path
from datetime import datetime, timezone
from html.parser import HTMLParser
import gzip, hashlib, json, re, shutil
import h5py
import numpy as np
import pandas as pd

OUT = (Path(__import__('os').environ['CRRM_WORKSPACE']) / 'Result6_analysis')
ROOT = Path(str(Path(__import__('os').environ['CRRM_DATA_ROOT']) / ''))
REV = ROOT / 'revisions/crrm_revision_20261111'
RESOURCE = Path(str(Path(__import__('os').environ.get('CRRM_RESOURCE_ROOT','external_resource')) / ''))
R4 = OUT.parent / 'Result4_analysis'
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()

def dump(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, default=lambda a: a.item() if isinstance(a, np.generic) else str(a))+'\n')

def soft_samples(path):
    rows = []
    for line in gzip.open(path, 'rt'):
        if line.startswith('^SAMPLE = '):
            rows.append({'geo_accession': line.strip().split(' = ', 1)[1]})
        elif rows and line.startswith('!Sample_') and ' = ' in line:
            key, value = line.strip().split(' = ', 1)
            key = key.removeprefix('!Sample_')
            if rows[-1].get(key) != value:
                rows[-1][key] = rows[-1].get(key, '') + ('; ' if key in rows[-1] else '') + value
    return pd.DataFrame(rows)

assert not (OUT/'audit/scoring_receipt.json').exists(), 'Do not redeclare a protocol after scoring.'
for name in ['frozen_CRRM388.csv', 'frozen_functional_branches.csv', 'score_gene_sets.json']:
    shutil.copy2(R4/'parameters'/name, OUT/'parameters'/name)
core = pd.read_csv(OUT/'parameters/frozen_CRRM388.csv')
sets = json.loads((OUT/'parameters/score_gene_sets.json').read_text())
assert len(core)==388 and core.human_ensembl.is_unique
assert core.frozen_direction.value_counts().to_dict()=={'up':280, 'down':108}
assert len(sets['Neutrophil degranulation'])==25 and len(sets['IFN-list overlap'])==79

# The deposited file is a normalized abundance matrix, despite its file name.
human_path = ROOT/'data/raw/geo_processed/GSE179353/GSE179353_Master_count_matrix_geo.csv.gz'
human = pd.read_csv(human_path, index_col=0)
assert human.gene_id.is_unique and human.columns.tolist()==['gene_id']+['AB'+str(i) for i in range(1,13)]
assert np.isfinite(human.iloc[:,1:].to_numpy()).all() and (human.iloc[:,1:].to_numpy()>=0).all()
hm = soft_samples(OUT/'data/GSE179353_family.soft.gz')
assert hm.geo_accession.tolist()==['GSM'+str(i) for i in range(5416085,5416097)]
hm['column'] = ['AB'+str(i) for i in range(1,13)]
hm['group'] = ['Uninfected baseline']*4+['RSV']*4+['RSV + BRD4i']*4
hm['source'] = 'BRD4_hSAEC'
hm['unit'] = hm.geo_accession
hm['celltype'] = 'Airway epithelium'
hm['replicate'] = hm.title.str.extract(r'replicate (\d+)')[0].astype(int)
hm['mapping_basis'] = 'Deposited AB1..AB12 column order matched to GEO ordered 12-sample design; evaluated against published non-CRRM marker patterns, not optimized using CRRM outcomes'
hm['pairing'] = 'Independent experimental replicates; matching replicate numbers do not establish pairing'
hm[['column','geo_accession','title','group','source','unit','celltype','replicate','mapping_basis','pairing','characteristics_ch1']].to_csv(OUT/'data/human_sample_manifest.csv', index=False)
human.to_csv(OUT/'data/human_author_TPM.csv.gz', index=False, compression='gzip')

# Independent cross-check of the ordered reconstruction against the source paper.
class Plain(HTMLParser):
    def __init__(self): super().__init__(); self.parts=[]
    def handle_data(self, s): self.parts.append(s)
parser=Plain(); parser.feed((OUT/'data/BRD4_article.html').read_text())
article=' '.join(parser.parts)
assert 'transcripts per million' in article
markers={'MX1':'ENSG00000157601','IFNL3':'ENSG00000197110','BRD4':'ENSG00000141867'}
checks=[]
for symbol, gene in markers.items():
    if gene in set(human.gene_id):
        values=human.set_index('gene_id').loc[gene]
        checks.append({'gene':symbol, 'gene_id':gene, **{g:float(values.loc[z.column].mean()) for g,z in hm.groupby('group',sort=False)}})
pd.DataFrame(checks).to_csv(OUT/'audit/ordered_column_mapping_source_marker_check.csv', index=False)

cells=pd.read_parquet(R4/'data/atlas_cell_metadata.parquet', filters=[('source_label','==','RelA_mouse')]).copy()
for col in cells:
    if isinstance(cells[col].dtype, pd.CategoricalDtype): cells[col]=cells[col].astype(object)
cells=cells.sort_values('expression_row_index').reset_index(drop=True)
assert len(cells)==84967 and cells.unit_resolved.all()
assert np.array_equal(cells.expression_row_index.to_numpy(), np.arange(len(cells)))
cells['unit']=cells.verified_unit_id
cells['group']=cells.infection_status.map({'corn oil':'Uninfected vehicle','corn oil-RSV':'RSV vehicle','TMX-RSV':'RSV RelA CKO'})
assert cells.group.notna().all()
cells['celltype']=cells.identity_level2
cells['source']='RelA_mouse'
keys=['source','unit','group','celltype']
codes, unique=pd.factorize(pd.MultiIndex.from_frame(cells[keys]),sort=True)
meta=pd.DataFrame(list(unique), columns=keys); meta['profile_id']=np.arange(len(meta)); meta['cells']=np.bincount(codes,minlength=len(meta))
meta.to_csv(OUT/'data/mouse_profile_metadata.csv',index=False)
cells['profile_id']=codes
cells[['cell_uid','expression_row_index','profile_id']].to_parquet(OUT/'data/mouse_row_membership.parquet', index=False)
mm=soft_samples(OUT/'data/GSE291797_family.soft.gz')
mm['retained_in_descriptor']=mm.geo_accession.isin(cells.unit.str.rsplit('::',n=1).str[-1].unique())
mm['group']=mm.characteristics_ch1.str.extract(r'treatment: ([^;]+)')[0].map({'corn oil':'Uninfected vehicle','corn oil-RSV':'RSV vehicle','TMX-RSV':'RSV RelA CKO'})
mm['analysis_unit']=mm.geo_accession
mm['status']=np.where(mm.retained_in_descriptor,'Included frozen source unit','Not present in frozen descriptor expression object; not restored or counted')
mm.to_csv(OUT/'data/mouse_sample_manifest.csv',index=False)
assert mm.loc[mm.retained_in_descriptor].group.value_counts().to_dict()=={'Uninfected vehicle':4,'RSV RelA CKO':4,'RSV vehicle':3}
assert mm.loc[~mm.retained_in_descriptor,'geo_accession'].tolist()==['GSM8842521']

counts=meta.pivot_table(index='celltype',columns='group',values='cells',aggfunc=lambda x:int((x>=30).sum()),fill_value=0)
primary_cells=sorted(counts.index[(counts>=3).all(axis=1)].tolist())
assert len(primary_cells)==11
counts.to_csv(OUT/'audit/mouse_lineage_unit_eligibility.csv')
mapping_path=REV/'phase05_mapping_corrected_reconstruction/mapping/babelgene22.9_mouse_human_one_to_one.csv'
mapping=pd.read_csv(mapping_path)
assert mapping.ensembl.is_unique and mapping.human_ensembl.is_unique
v=pd.read_parquet(RESOURCE/'DATA/curated_metadata/features/RelA_mouse_var.parquet')
mapping=core.merge(mapping[['ensembl','human_ensembl']],on='human_ensembl',how='left',validate='one_to_one')
mapping['native_feature_present']=mapping.ensembl.isin(v.native_gene_id)
mapping.to_csv(OUT/'parameters/frozen_mouse_member_projection.csv',index=False)

comparisons=[dict(key='BRD4_Airway', source='BRD4_hSAEC',celltype='Airway epithelium',case='RSV + BRD4i',control='RSV',context_case='RSV',context_control='Uninfected baseline',n_case=4,n_control=4)]
for ct in primary_cells:
    comparisons.append(dict(key='RelA_'+ct,source='RelA_mouse',celltype=ct,case='RSV RelA CKO',control='RSV vehicle',context_case='RSV vehicle',context_control='Uninfected vehicle',n_case=4,n_control=3))

plan=dict(version='host_perturbation_r1', declared_at=datetime.now(timezone.utc).isoformat(),status='FROZEN_BEFORE_NEW_CRRM388_SCORING',
    objective='Test perturbation responsiveness and direct component differences within frozen CRRM388; use clinical associations as a separate cross-context comparison',
    design={'human':'12 independent experimental replicates, not 12 donors; unpaired 4 vs 4 intervention comparison. Baseline deposited as 0h; no uninfected drug arm. No intervention-by-infection interaction or restoration claim.',
            'mouse':'Frozen 84967-cell source, 11 lung source units (4 uninfected vehicle, 3 RSV vehicle, 4 RSV TMX/RelA CKO), day21; no uninfected TMX arm in these transcriptomes. Not an acute replication of the human drug contrast.',
            'missing_mouse_sample':'GSM8842521 / library225 is in GEO but absent from the frozen descriptor; retain existing resource subset without adding a fourth infected control.',
            'human_columns':'Order reconstruction AB1..4 baseline, AB5..8 RSV, AB9..12 RSV+ZL0454; retained as an explicit metadata assumption and corroborated against published marker patterns.'},
    normalization={'human':'log2(1+author-deposited normalized abundance/TPM); never DESeq2 on this decimal matrix; all12 profiles form one study-specific gene mean/sample-SD reference.',
                   'mouse':'native probe counts -> per-cell log1p(10000*host_counts/total_host_counts) -> mean per animal and canonical lineage; all11 prespecified complete lineages form one shared study reference with equal total weight per animal; population weighted SD, guard1e-8.',
                   'rank':'Per-profile signed within-transcriptome average percentile rank; distinct score units; no gene-z/rank effect pooling.',
                   'eligibility':'total >=70%388 and >=50% each original direction; branches >=5 and >=50% frozen members; primary mouse >=30 cells/profile and >=3 animals/group; thresholds15/50 are separate checks.'},
    frozen_gene_count=388, directions={'up':280,'down':108},annotations='Seven original overlapping subsets and original 79/309 IFN-list partition; no new signature or outcome-selected members',
    primary_lineages=primary_cells,comparisons=comparisons,
    correction_families={'intervention_total':'BH12 across human plus all11 eligible mouse lineages',
       'intervention_branches':'BH84 across7 fixed branches x12 contrasts',
       'intervention_IFN_components':'BH24 across2 components x12 contrasts',
       'intervention_contributions':'BH24 across2 signed contributions x12 contrasts',
       'infection_context':'Separate corresponding BH12/BH84/BH24 families; human baseline time differs as deposited',
       'direct_component_comparisons':'BH24 across (TNF/NFkB minus IFN-list overlap) and (degranulation annotation minus IFN-list overlap), in each of12 intervention contrasts; compare per-profile component differences, not separate significance labels',
       'clinical_gene_alignment':'BH3 across the two frozen Severe-minus-Mild gene vectors and paired recovery-minus-acute vector; exact human treatment-label permutation of 8 profiles (70 allocations), retaining correlated genes as a vector; no gene-as-independent-replicate P values',
       'gene_tests':'BH4656 for388 members x12 intervention contrasts, retaining unavailable tests in denominator',
       'sensitivities':'Each method repeats its corresponding full prespecified family; no best-method selection or reduced displayed-row family'},
    inference='Unpaired Welch t and pointwise95%t intervals at independent culture/source-animal level; exact whole-unit studentized permutation as corroboration; both-group near-constant signed values <=1e-12 yield P=1, never numerical significance',
    sensitivities=['Signed transcriptome rank','Human total-abundance rescaling before log','Minimum15/50 mouse cells','Degranulation without LCN2','TNF/degranulation outside IFN lists: direct comparisons to IFN overlap','Every frozen-gene omission influence','Detection-supported member sensitivity for mouse'],
    clinical_bridge='Fixed25-member annotation subset only; actual Result3 frozen severe/mild log expression effects; recovery effects reconstructed from the15 previously verified initially moderate/severe participant pairs. Correlation across gene-effect vectors evaluated by intervention-unit permutation. Not transportable clinical intervention efficacy.',
    interpretation='Possible broad suppression, selective attenuation, persistence or opposite shifts are all reportable. Outside IFN lists is annotation membership, not IFN independence; degranulation annotation is not measured degranulation activity. Mouse TMX contrast is the intervention-associated result under this design; nonrecombined lineage changes are not cell-autonomous proof.')
dump(OUT/'parameters/analysis_plan_v1.json',plan)
inputs=[human_path,mapping_path,R4/'parameters/frozen_CRRM388.csv',R4/'parameters/frozen_functional_branches.csv',R4/'parameters/score_gene_sets.json',RESOURCE/'DATA/curated_metadata/features/RelA_mouse_var.parquet',OUT/'data/GSE179353_family.soft.gz',OUT/'data/GSE291797_family.soft.gz']
dump(OUT/'audit/input_audit.json',dict(status='PASS_WITH_DECLARED_SCOPE',human_profile_count=12,human_library_sums=human.iloc[:,1:].sum().to_dict(),human_matrix_type='Normalized abundance/TPM',mouse_cells=len(cells),mouse_units=11,mouse_GEO_samples=12,missing_frozen_mouse_GSM='GSM8842521',mouse_primary_lineages=primary_cells,mapped_mouse_core_present=int(mapping.native_feature_present.sum()),inputs=[dict(path=str(p),sha256=sha(p)) for p in inputs],source_markers=checks))
dump(OUT/'parameters/declaration_receipt.json',dict(declared_at=plan['declared_at'],plan_sha256=sha(OUT/'parameters/analysis_plan_v1.json'),new_scores_computed=False))
print(json.dumps(dict(primary_contrasts=len(comparisons),mouse_lineages=primary_cells,mouse_mapped_core=int(mapping.native_feature_present.sum()),source_marker_check=checks),indent=2))
