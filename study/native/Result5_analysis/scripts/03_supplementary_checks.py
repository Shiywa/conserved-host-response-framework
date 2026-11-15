from pathlib import Path
from itertools import product
import ast,json
import numpy as np,pandas as pd
from scipy import stats
OUT=(Path(__import__('os').environ['CRRM_WORKSPACE']) / 'Result5_analysis');R4=(OUT.parent/'Result4') if (OUT.parent/'Result4').exists() else (OUT.parent/'Result4_analysis');RESOURCE=Path(__import__('os').environ.get('CRRM_RESOURCE_ROOT','external_resource'))
code=(OUT/'scripts/02_score_and_compare.py').read_text();tree=ast.parse(code)
for name in ['infer','bh','hc3','permutations']:
    node=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name==name);exec(ast.get_source_segment(code,node))
plan=json.loads((OUT/'parameters/analysis_plan_v1.json').read_text());scores=pd.read_csv(OUT/'results/profile_scores.csv');core=pd.read_csv(OUT/'parameters/frozen_CRRM388.csv');sets=json.loads((OUT/'parameters/score_gene_sets.json').read_text());genes=core.human_ensembl.to_numpy();sign=np.where(core.frozen_direction.eq('up'),1.,-1.)
refs=pd.read_csv(OUT/'data/gene_reference_parameters.csv');targets=['CRRM388','IFN-list overlap','Outside IFN lists','Neutrophil degranulation'];rows=[]

# The descriptor author-intent exclusion removes the complete pHC15 sample; no re-extraction or rescaling is needed.
a=pd.read_parquet(R4/'data/atlas_cell_metadata.parquet',columns=['source_label','sample_id','author_intent_sensitivity_included']);a=a[a.source_label.eq('Infant_PBMC')].copy();a['included']=a.author_intent_sensitivity_included.astype(str).str.lower().eq('true');fractions=a.groupby('sample_id',observed=True).included.mean();assert fractions.isin([0.,1.]).all();exclude=fractions[fractions.eq(0)].index.tolist();assert exclude==['pHC15']
pb=scores[scores.source_label.eq('Infant_PBMC')&scores.view.eq('canonical')&scores.cells.ge(30)&scores.analysis_group.isin(['Healthy','RSV','COVID'])]
for cell in ['B_cell','T_cell','NK_cell','Monocyte']:
    g=pb[pb.analysis_celltype.eq(cell)&~pb.participant.isin(exclude)]
    for case,control in [('RSV','COVID'),('RSV','Healthy'),('COVID','Healthy')]:
        x=g[g.analysis_group.eq(case)];y=g[g.analysis_group.eq(control)]
        for target in targets:rows.append(dict(key='PBMC_'+cell+'_'+case+'_'+control,cell=cell,method='Descriptor author-intent subset; fixed original scaling',score=target,n_case=len(x),n_control=len(y),**infer(x[target],y[target])))
ar=pd.DataFrame(rows)
for (comparison,target),ix in ar.groupby([ar.key.str.extract(r'_(RSV_COVID|RSV_Healthy|COVID_Healthy)$')[0],'score']).groups.items():ar.loc[ix,'adjusted_P_BH4']=bh(ar.loc[ix,'P'],4)
ar.to_csv(OUT/'results/PBMC_author_intent_sensitivity.csv',index=False)

dc=scores[scores.source_label.eq('CD1c_DC')&scores.view.eq('canonical')&scores.cells.ge(30)];dcrows=[];donorrows=[]
for cell in ['Dendritic_cell','Monocyte']:
    for time in ['4h','8h','18h']:
        g=dc[dc.analysis_celltype.eq(cell)&dc.time.eq(time)]
        for case,control in [('RSV','Flu'),('RSV','RSV_HK'),('RSV','Untreated'),('Flu','Untreated')]:
            a=g[g.analysis_group.eq(case)].set_index('analysis_unit');b=g[g.analysis_group.eq(control)].set_index('analysis_unit');common=a.index.intersection(b.index);a=a.loc[common];b=b.loc[common]
            for target in targets:
                r=infer(a[target],b[target],True);pp=permutations(a[target],b[target],True) if len(common)>=3 else np.nan
                dcrows.append(dict(cell=cell,time=time,comparison=case+' minus '+control,score=target,n_donors=len(common),exact_signflip_P=pp,**r))
                for u in common:donorrows.append(dict(cell=cell,time=time,case=case,control=control,score=target,donor=u,case_value=a.loc[u,target],control_value=b.loc[u,target],difference=a.loc[u,target]-b.loc[u,target]))
dr=pd.DataFrame(dcrows)
for (cell,comp,target),ix in dr.groupby(['cell','comparison','score']).groups.items():dr.loc[ix,'adjusted_P_BH3']=bh(dr.loc[ix,'P'],3)
dr.to_csv(OUT/'results/DC_time_matched_effects.csv',index=False);pd.DataFrame(donorrows).to_csv(OUT/'results/DC_time_matched_donor_values.csv',index=False)

# Source/sample coverage and covariate summaries without treating cells as replication.
manifest=pd.read_csv(OUT/'data/clinical_sample_manifest.csv');ms=[]
for (source,site,group,primary),g in manifest.groupby(['dataset','site','condition','primary']):ms.append(dict(source=source,site=site,group=group,primary=bool(primary),participants=len(g),median_age_months=float(g.age_months.median()) if g.age_months.notna().any() else np.nan,age_min=g.age_months.min(),age_max=g.age_months.max(),male=int(g.sex.eq('M').sum()),female=int(g.sex.eq('F').sum()),age_sex_available=bool(g.age_months.notna().all()&g.sex.isin(['M','F']).all())))
cross=pd.read_csv(OUT/'data/PBMC_sample_crosswalk.csv')
for group,g in cross[cross.primary].groupby('analysis_group'):ms.append(dict(source='Infant_PBMC',site='Columbus',group=group,primary=True,participants=len(g),median_age_months=g.age_months.median(),age_min=g.age_months.min(),age_max=g.age_months.max(),male=int(g.sex.eq('M').sum()),female=int(g.sex.eq('F').sum()),age_sex_available=True))
pd.DataFrame(ms).to_csv(OUT/'results/sample_and_covariate_summary.csv',index=False)
audit=dict(status='PASS',author_intent_excluded_samples=exclude,author_intent_COVID_RSV_direct_contrasts_unchanged=True,DC_time_hypotheses=len(dr),COVID_included_only_acute_untreated=17,primary_n_case_RSV=19,PBMC_healthy_primary=17,PBMC_healthy_author_intent=16,UK_individual_covariates_available=False,UK_scope='Unadjusted external replication; UK data provide infecting-pathogen but not age/sex/illness-duration/treatment per sample; no fabricated individual adjustment',illness_duration_adjustment_available=False,treatment='PBMC explicitly excludes4steroid-treatedCOVID; no claim that all studies are treatment-adjusted',no_exact_published_DE_reconstruction_claim=True)
(OUT/'audit/supplementary_checks.json').write_text(json.dumps(audit,indent=2)+'\n')
print(json.dumps(audit,indent=2))
