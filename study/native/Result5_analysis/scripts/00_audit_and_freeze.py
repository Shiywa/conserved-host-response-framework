from pathlib import Path
from datetime import datetime, timezone
import gzip,csv,json,hashlib,shutil,re
import numpy as np,pandas as pd,h5py

OUT=(Path(__import__('os').environ['CRRM_WORKSPACE']) / 'Result5_analysis')
ROOT=Path(str(Path(__import__('os').environ['CRRM_DATA_ROOT']) / ''))
RESOURCE=Path(str(Path(__import__('os').environ.get('CRRM_RESOURCE_ROOT','external_resource')) / ''))
R4=(OUT.parent/'Result4_analysis') if (OUT.parent/'Result4_analysis').exists() else RESOURCE/'CRRM_revision_work/Result4'
def dump(p,x):p.write_text(json.dumps(x,indent=2,ensure_ascii=False,default=str)+'\n')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def header(path):
    rows=[]
    with gzip.open(path,'rt') as f:
        for line in f:
            if line.startswith('!series_matrix_table_begin'):break
            if line.startswith('!Sample_'):rows.append(next(csv.reader([line],delimiter='\t')))
    m={};seen={}
    for row in rows:
        name=row[0].removeprefix('!Sample_');seen[name]=seen.get(name,0)+1
        key=name if seen[name]==1 else name+'_'+str(seen[name]);m[key]=row[1:]
    return pd.DataFrame(m)
def unpack(m):
    for c in [c for c in m if c.startswith('characteristics_ch1')]:
        for i,v in m[c].items():
            if ': ' in v:
                k,val=v.split(': ',1);m.loc[i,'cov_'+k.lower().strip()]=val
    return m
for name in ['frozen_CRRM388.csv','frozen_functional_branches.csv','score_gene_sets.json']:
    shutil.copy2(R4/'parameters'/name,OUT/'parameters'/name)
inputs=[];clinical=[]
for platform in ['GPL6884','GPL10558']:
    p=ROOT/f'data/raw/GSE38900/GSE38900-{platform}_series_matrix.txt.gz'
    m=unpack(header(p));m['dataset']='GSE38900_'+platform;m['platform']=platform
    src=m.source_name_ch1.str.lower();m['site']=np.select([src.str.contains('dall'),src.str.contains('finn'),src.str.contains('columbus')],['Dallas','Finland','Columbus'],default='Unknown')
    m['condition']=np.select([src.str.contains('follow'),src.str.contains('healthy'),src.str.contains('influenza'),src.str.contains('hrv'),src.str.contains('rsv')],['Recovery','Healthy','Influenza','Rhinovirus','RSV'],default='Other')
    m['age_months']=pd.to_numeric(m.get('cov_age (months)',m.get('cov_age')),errors='coerce');m['sex']=m.cov_gender.str.lower().map({'female':'F','male':'M'})
    m['participant']=m.description.str.strip();m['primary']=m.site.eq('Dallas')&m.condition.isin(['Healthy','RSV','Influenza','Rhinovirus'])
    m['role']=np.where(m.primary,'Primary same-site multivirus comparison',np.where(m.condition.eq('Recovery'),'Excluded repeated follow-up','RSV-only context; no direct multivirus comparison'))
    clinical.append(m);inputs.append(dict(path=str(p),sha256=sha(p)))
p=OUT/'data/GSE42026_series_matrix.txt.gz';m=unpack(header(p));m['dataset']='GSE42026';m['site']='UK';m['platform']='GPL6947'
src=m.source_name_ch1.str.lower();text=m.title.str.lower()+' '+src+' '+m[[c for c in m if c.startswith('characteristics')]].astype(str).agg(' '.join,axis=1).str.lower()
m['condition']=np.select([text.str.contains('healthy|control'),text.str.contains('h1n1|influenza'),text.str.contains('rsv'),text.str.contains('bacter')],['Healthy','Influenza','RSV','Bacterial'],default='Other')
m['age_months']=np.nan;m['sex']='Unknown';m['participant']=m.geo_accession;m['primary']=m.condition.isin(['Healthy','RSV','Influenza']);m['role']=np.where(m.primary,'Independent external replication','Excluded non-target infection')
clinical.append(m);inputs.append(dict(path=str(p),sha256=sha(p)))
cm=pd.concat(clinical,ignore_index=True);cm.to_csv(OUT/'data/clinical_sample_manifest.csv',index=False)
assert cm.geo_accession.is_unique
assert cm[(cm.dataset=='GSE38900_GPL6884')&cm.primary].groupby('condition').size().to_dict()=={'Healthy':27,'Influenza':16,'RSV':91,'Rhinovirus':30}

# Exact source-native sample/GSM crosswalk, including COVID participants unresolved in the descriptor's RSV-focused registry.
geo=pd.read_csv(R4/'data/Infant_PBMC_public_GEO_covariates.csv');geo['participant']=geo.title.str.extract(r'(p(?:RSV|HC|CoV)\d+)')[0]
geo['age_months']=pd.to_numeric(geo['age (months):ch1'],errors='coerce');geo['sex']=geo['Sex:ch1'].str.lower().map({'female':'F','male':'M'})
def decode(h,name):
    a=h['obs'][name]
    if isinstance(a,h5py.Group):
        cats=a['categories'][:].astype(str);codes=a['codes'][:];out=np.full(len(codes),'',dtype=object);ok=codes>=0;out[ok]=cats[codes[ok]];return out
    return a[:].astype(str)
with h5py.File(RESOURCE/'DATA/expression/Infant_PBMC.h5ad') as h:
    cross=pd.DataFrame({c:decode(h,c) for c in ['sample_id','library_id','infection_status']}).drop_duplicates()
cross=cross.merge(geo[['geo_accession','participant','title','age_months','sex']],left_on='library_id',right_on='geo_accession',validate='one_to_one')
assert len(cross)==68 and cross.library_id.is_unique
known=cross.participant.notna();assert (cross.loc[known,'sample_id']==cross.loc[known,'participant']).all()
cross['analysis_group']=np.select([cross.infection_status.eq('Healthy'),cross.infection_status.isin(['Mild RSV','Moderate RSV','Severe RSV']),cross.infection_status.isin(['Moderate COVID-19','Severe COVID-19'])],['Healthy','RSV','COVID'],default='Excluded')
cross['primary']=cross.analysis_group.ne('Excluded');cross['exclusion_reason']=np.select([cross.infection_status.str.contains('Convalescent'),cross.infection_status.str.contains('Steroid')],['Convalescent/subacute COVID; timing differs','Steroid-treated COVID; treatment differs'],default='Other or no exclusion')
assert cross[cross.primary].groupby('analysis_group').size().to_dict()=={'COVID':17,'Healthy':17,'RSV':19}
cross.to_csv(OUT/'data/PBMC_sample_crosswalk.csv',index=False)
gp=pd.read_csv(R4/'data/profile_metadata.csv');gp=gp[gp.source_label.isin(['Infant_PBMC','CD1c_DC'])].copy()
libmap=cross.set_index('library_id');partmap=cross.set_index('participant',drop=False)
for i,r in gp[gp.source_label.eq('Infant_PBMC')].iterrows():
    tail=r.analysis_unit.split('::')[-1];q=libmap.loc[tail] if tail.startswith('GSM') else partmap.loc[tail]
    assert q.infection_status==r.infection_status
    gp.loc[i,'participant']=q.participant;gp.loc[i,'analysis_group']=q.analysis_group;gp.loc[i,'age_months']=q.age_months;gp.loc[i,'sex']=q.sex
    if q.primary:gp.loc[i,'analysis_unit']='Infant_PBMC::'+q.participant;gp.loc[i,'unit_resolved']=True
gp.loc[gp.source_label.eq('CD1c_DC'),'analysis_group']=gp.loc[gp.source_label.eq('CD1c_DC'),'infection_status']
gp.to_csv(OUT/'data/profile_metadata.csv',index=False)

specs=[]
def add(key,source,cell,case,control,paired=False,family='Direct virus total'):
    specs.append(dict(key=key,source=source,cell=cell,case=case,control=control,paired=paired,total_family=family))
for v in ['RSV','Influenza','Rhinovirus']:add('Dallas_'+v+'_Healthy','GSE38900_GPL6884','Whole blood',v,'Healthy',family='Clinical shared induction')
for v in ['Influenza','Rhinovirus']:add('Dallas_RSV_'+v,'GSE38900_GPL6884','Whole blood','RSV',v)
for cell in ['B_cell','T_cell','NK_cell','Monocyte']:
    add('PBMC_'+cell+'_RSV_COVID','Infant_PBMC',cell,'RSV','COVID')
    for v in ['RSV','COVID']:add('PBMC_'+cell+'_'+v+'_Healthy','Infant_PBMC',cell,v,'Healthy',family='PBMC shared induction')
add('DC_RSV_Flu','CD1c_DC','Dendritic_cell','RSV','Flu',True)
for v in ['RSV','Flu']:add('DC_'+v+'_Untreated','CD1c_DC','Dendritic_cell',v,'Untreated',True,family='DC shared induction')
add('DC_RSV_HK','CD1c_DC','Dendritic_cell','RSV','RSV_HK',True,family='Heat-killed control')
for a,b in [('RSV','Healthy'),('Influenza','Healthy'),('RSV','Influenza')]:add('UK_'+a+'_'+b,'GSE42026','Whole blood',a,b,family='Independent replication')
direct=[s['key'] for s in specs if s['total_family']=='Direct virus total'];assert len(direct)==7
plan=dict(status='DECLARED_BEFORE_NEW_CROSS_VIRUS_CRRM388_SCORING',declared_at=datetime.now(timezone.utc).isoformat(),question='Shared respiratory-virus response versus virus-associated differences in the frozen CRRM; no diagnostic model or gene-set rebuilding',retrospective=True,known_results='Results1–4 and historical70-gene analysis already known; no new388 multivirus scores/effects inspected before this declaration',core_members=388,up=280,down=108,branches=7,signatures='Unchanged frozen members/signs/overlapping annotations and Hallmark membership components',specifications=specs,
    primary_clinical='GSE38900 Dallas GPL6884 only:91 acuteRSV,16influenza,30rhinovirus,27healthy; Finland/Columbus not mixed into direct contrasts;21follow-ups excluded',
    independent_replication='GSE42026 UK whole blood, independent recruitment institution and years; use exact deposited group labels; no forced pooling with Dallas; GSE34205 not needed as independent validation because age/batch and potential study-network overlap require separate resolution',
    PBMC='Exact68-sample native GSM/sample identity reconciled to public GEO titles; acute untreated COVID=12moderate+5severe;19RSV;17healthy. Exclude10convalescent,4steroid-treated,1Other. Keep descriptor RNA scope; author-intent subset sensitivity explicitly retained.',
    DC='Primary canonical Dendritic_cell,3donors; within-donor and same4/8/18h subtraction, equal average times per donor; common complete3-time pairs primary. Donors are n, not9donor-time strata. Monocyte supplementary if eligible; candidate basophil not used.',
    scoring=dict(arrays='Deposited normalized expression; log2(max(x,0)+1)if99th percentile>100, median over unambiguous probes per gene; no outcome-based probe selection',single_cell='Reuse Result4 mean log1p(10000*native_count/total_host_counts) and exact profile membership; do not recompute cell annotation or alter raw resources',reference='One shared within-study gene mean/SD across primary groups, equal total weight per biological unit; canonical eligible profiles>=30cells. PBMC scaling now includes acute untreatedCOVID together with HC/RSV, a deliberate cross-virus reference distinct from Result4. No separately fitted virus scales.',eligibility='SD>1e-8;total>=70%388,>=50%each arm;subsets>=5and>=50%nominal. Missing features never zero. Same feature set/reference for all contrasts within a source.',contributions='Signed means;up/down member-count fractions add exactly to total. PositiveRTK/microtubule means movement in original down-defined direction, not pathway activation.'),
    inference=dict(primary='Welch participant contrasts with95%tCI; paired donor t test on time-balanced donor deltas with95%tCI; two-sided',multiplicity='Total direct-virus family BH7 across2Dallas+4PBMC+1DC contrasts. Branch differences separateBH49,IFN-list componentsBH14,up/downcontributionsBH14. Clinical shared-inductionBH3,PBMCsharedBH8,DCsharedBH2; independent replication totalBH3 and all7branchesx3contrasts BH21. Unavailable planned hypotheses retainedP=1 in family.',resampling='19999studentized participant-label permutations(seed20261113)for primary direct total; exact sign flips for3donor DC; small-n exact minimumP=.25 is reported alongside model inference',gene_level='All388 frozen hypotheses retained in BH776 family for2Dallasdirectvirus contrasts; unmeasuredP=1; no significance-in-one-group specificity rule. Lock entire discovery directional gene findings for UKRSVvsflu test before computing UKgene effects; validationBH388 retains full frozen family. PBMC/DC genes descriptive context only unless separate declared family.',sensitivities='Age/sexHC3 within clinical group-pairs where individual data exist;common-age support; signed transcriptome-rank score without gene variance scaling; common-measured core;cell thresholds15/50; author-intent PBMC subset; four-lineage composition standardization;leave-one-participant/donor-out and gene-deletion leverage; detection-supported/PIGR check. Complete families retained.'),
    composition='Observed versus fixed equal-healthy-participant four-lineage PBMC composition,restrictedB/T/NK/Monocyte;stratified participant bootstrap2000percentileCI,not cell bootstrap and not whole-blood mediation',interpretation='No equivalence from nonsignificance, no IFN-independent term, no direct degranulation activity claim, no cross-study virus effect ranking, no claiming universal virus or RSV specificity from limited comparators',input_metadata_summary=cm.groupby(['dataset','site','condition','primary']).size().rename('samples').reset_index().to_dict('records'),inputs=inputs,frozen_sha256=sha(OUT/'parameters/frozen_CRRM388.csv'))
dest=OUT/'parameters/analysis_plan_v1.json';assert not dest.exists(),'Do not replace a declared plan';dump(dest,plan)
dump(OUT/'parameters/declaration_receipt.json',dict(status='FROZEN',declared_at=plan['declared_at'],plan_sha256=sha(dest),new_cross_virus_scores_computed=False))
dump(OUT/'audit/sample_audit.json',dict(status='PASS',clinical_counts=plan['input_metadata_summary'],PBMC_counts=cross.groupby(['infection_status','analysis_group']).size().rename('participants').reset_index().to_dict('records'),PBMC_exact_native_to_GEO_links=68,primary_direct_hypotheses=7,no_source_raw_data_modified=True))
print(json.dumps(dict(clinical=plan['input_metadata_summary'],PBMC_primary=cross[cross.primary].groupby('analysis_group').size().to_dict(),status='FROZEN'),indent=2))
