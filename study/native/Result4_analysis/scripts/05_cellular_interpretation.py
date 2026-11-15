from pathlib import Path
from datetime import datetime,timezone
import ast,json
import numpy as np,pandas as pd
from scipy import stats
from itertools import combinations
from math import comb

OUT=(Path(__import__('os').environ['CRRM_WORKSPACE']) / 'Result4_analysis');plan=json.loads((OUT/'parameters/analysis_plan_v1.json').read_text());scores=pd.read_csv(OUT/'results/profile_scores.csv');core=pd.read_csv(OUT/'parameters/frozen_CRRM388.csv');sets=json.loads((OUT/'parameters/score_gene_sets.json').read_text());branches=pd.read_csv(OUT/'parameters/frozen_functional_branches.csv').branch.drop_duplicates().tolist();permutation_cache={}
src=(OUT/'scripts/04_score_and_test_profiles.py').read_text();tree=ast.parse(src)
for name in ['unpaired_allocations','infer','contrast_rows','bh']:
    node=next(x for x in tree.body if isinstance(x,ast.FunctionDef) and x.name==name);exec(ast.get_source_segment(src,node))
def dump(p,x):p.write_text(json.dumps(x,ensure_ascii=False,indent=2,default=lambda x:x.item() if isinstance(x,np.generic) else str(x))+'\n')
def available_mean(x,axis=0):
    count=np.isfinite(x).sum(axis=axis);total=np.nansum(x,axis=axis);out=np.full_like(total,np.nan,dtype=float);np.divide(total,count,out=out,where=count>0);return out

sens=[]
for threshold in [15,50]:
    for spec in plan['primary_comparisons']:
        for target in ['CRRM388','Outside IFN lists','Degranulation without LCN2']+branches:
            row,_=contrast_rows(spec,target,min_cells=threshold);row['scope']='Minimum-cell sensitivity';sens.append(row)
s=pd.DataFrame(sens)
for (threshold,targetkind),ix in s.groupby(['minimum_cells',s.score.map(lambda x:'branches' if x in branches else x)]).groups.items():s.loc[ix,'correction_family']=f'Threshold{threshold} | '+str(targetkind);s.loc[ix,'adjusted_P']=bh(s.loc[ix,'P'])
s.to_csv(OUT/'results/minimum_cell_sensitivities.csv',index=False)

def ols_hc3(y,x,coefficient=1):
    if len(y)-x.shape[1]<5 or np.linalg.matrix_rank(x)<x.shape[1]:return dict(status='UNIDENTIFIABLE_MODEL',estimate=np.nan,lower=np.nan,upper=np.nan,P=np.nan)
    beta=np.linalg.lstsq(x,y,rcond=None)[0];xx=np.linalg.inv(x.T@x);res=y-x@beta;h=np.einsum('ij,jk,ik->i',x,xx,x);u=res/np.maximum(1-h,1e-8);cov=xx@((x.T*(u*u))@x)@xx;se=np.sqrt(cov[coefficient,coefficient]);df=len(y)-x.shape[1];estimate=beta[coefficient];ci=stats.t.ppf(.975,df)*se
    return dict(status='OK',estimate=float(estimate),lower=float(estimate-ci),upper=float(estimate+ci),P=float(2*stats.t.sf(abs(estimate/se),df)))
cov=pd.read_csv(OUT/'data/patient_covariates_internal.csv');adjusted=[]
for spec in plan['primary_comparisons']:
    if spec['source']!='Infant_PBMC':continue
    for endpoint,case,control in [('RSV_vs_Healthy',['Mild RSV','Moderate RSV','Severe RSV'],['Healthy']),('Severe_vs_Mild',['Severe RSV'],['Mild RSV'])]:
        g=scores[(scores.source_label=='Infant_PBMC')&(scores.view=='canonical')&(scores.analysis_celltype==spec['celltype'])&scores.unit_resolved&(scores.cells>=30)&scores.infection_status.isin(case+control)].merge(cov,on='analysis_unit',validate='one_to_one')
        for target in list(dict.fromkeys(['CRRM388','Outside IFN lists','Degranulation without LCN2']+branches)):
            gg=g.dropna(subset=[target,'age_months','sex']);x=np.column_stack([np.ones(len(gg)),gg.infection_status.isin(case).astype(int),gg.age_months.to_numpy(),(gg.sex=='M').astype(int)]);r=ols_hc3(gg[target].to_numpy(),x);adjusted.append(dict(endpoint=endpoint,cell_identity=spec['celltype'],score=target,case_units=int(gg.infection_status.isin(case).sum()),control_units=int(gg.infection_status.isin(control).sum()),covariates='public GEO age in months and sex',**r))
a=pd.DataFrame(adjusted)
for (endpoint,targetkind),ix in a.groupby(['endpoint',a.score.map(lambda x:'branches' if x in branches else x)]).groups.items():a.loc[ix,'correction_family']=endpoint+' | age/sex adjusted | '+str(targetkind);a.loc[ix,'adjusted_P']=bh(a.loc[ix,'P'])
a.to_csv(OUT/'results/PBMC_covariate_adjusted_effects.csv',index=False)

types=['B_cell','T_cell','NK_cell','Monocyte'];targets=list(dict.fromkeys(['CRRM388','IFN-list overlap','Outside IFN lists','Degranulation without LCN2']+branches))
pb=scores[(scores.source_label=='Infant_PBMC')&(scores.view=='canonical')&scores.unit_resolved&scores.infection_status.isin(['Healthy','Mild RSV','Moderate RSV','Severe RSV'])]
units=sorted(pb.analysis_unit.unique());gg=pb[pb.analysis_celltype.isin(types)].set_index(['analysis_unit','analysis_celltype']).reindex(pd.MultiIndex.from_product([units,types],names=['analysis_unit','analysis_celltype']));assert len(units)==36 and gg.cells.notna().all() and (gg.cells>=30).all()
counts=gg.cells.to_numpy().reshape(36,4);fractions=counts/counts.sum(axis=1,keepdims=True);patient=pd.DataFrame({'analysis_unit':units});patient=patient.merge(pb[['analysis_unit','infection_status']].drop_duplicates(),validate='one_to_one');conditions=patient.infection_status.to_numpy();healthy=np.flatnonzero(conditions=='Healthy');pref=fractions[healthy].mean(axis=0);cube=gg[targets].to_numpy().reshape(36,4,len(targets));observed=np.einsum('ik,ikt->it',fractions,cube);fixed=np.einsum('k,ikt->it',pref,cube)
allcounts=pb.groupby('analysis_unit').cells.sum().reindex(units).to_numpy();patient['four_lineage_cells']=counts.sum(axis=1);patient['all_observed_PBMC_cells']=allcounts;patient['four_lineage_fraction']=counts.sum(axis=1)/allcounts
for k,name in enumerate(types):patient[name+'_fraction']=fractions[:,k]
for j,target in enumerate(targets):patient[target+'_observed']=observed[:,j];patient[target+'_fixed_composition']=fixed[:,j]
patient.to_csv(OUT/'results/PBMC_patient_composition_and_scores.csv',index=False)
pd.DataFrame({'cell_identity':types,'reference_healthy_fraction':pref}).to_csv(OUT/'results/PBMC_fixed_reference_composition.csv',index=False)
composition=[];rng=np.random.default_rng(20261113)
for endpoint,case,control in [('RSV_vs_Healthy',['Mild RSV','Moderate RSV','Severe RSV'],['Healthy']),('Severe_vs_Mild',['Severe RSV'],['Mild RSV'])]:
    ia=np.flatnonzero(np.isin(conditions,case));ib=np.flatnonzero(np.isin(conditions,control));boot_observed=[];boot_fixed=[]
    for _ in range(2000):
        aa=rng.choice(ia,len(ia),replace=True);bb=rng.choice(ib,len(ib),replace=True);hh=bb if endpoint=='RSV_vs_Healthy' else rng.choice(healthy,len(healthy),replace=True);pr=fractions[hh].mean(axis=0);sf=np.einsum('k,ikt->it',pr,cube);boot_observed.append(observed[aa].mean(axis=0)-observed[bb].mean(axis=0));boot_fixed.append(sf[aa].mean(axis=0)-sf[bb].mean(axis=0))
    bo=np.array(boot_observed);bf=np.array(boot_fixed)
    for j,target in enumerate(targets):
        for label,matrix,boots in [('Observed composition',observed,bo),('Fixed healthy composition',fixed,bf)]:
            va=matrix[ia,j];vb=matrix[ib,j]
            if not np.isfinite(np.r_[va,vb]).all():continue
            r=infer(va,vb);low,high=np.quantile(boots[:,j],[.025,.975]);composition.append(dict(endpoint=endpoint,score=target,model=label,case_units=len(ia),control_units=len(ib),estimate=r['estimate'],lower=low,upper=high,interval='95% stratified-patient bootstrap;2000 draws',P=r['P'],welch_lower=r['lower'],welch_upper=r['upper'],four_lineage_population='B/T/NK/Monocyte',minimum_four_lineage_fraction=float(patient.four_lineage_fraction.min())))
        diff=(observed[ia,j]-fixed[ia,j]).mean()-(observed[ib,j]-fixed[ib,j]).mean();lo,hi=np.quantile((bo-bf)[:,j],[.025,.975]);composition.append(dict(endpoint=endpoint,score=target,model='Observed minus fixed effect; descriptive',case_units=len(ia),control_units=len(ib),estimate=diff,lower=lo,upper=hi,interval='95% joint stratified-patient bootstrap',P=np.nan,four_lineage_population='B/T/NK/Monocyte'))
pd.DataFrame(composition).to_csv(OUT/'results/PBMC_composition_standardization.csv',index=False)

# Fixed-gene contribution descriptions are explanatory follow-up, not a new gene-set selection.
references=pd.read_csv(OUT/'data/gene_reference_parameters.csv');contrib=[];sgn=np.where(core.frozen_direction=='up',1.,-1.);degran=set(sets['Neutrophil degranulation']);outside=set(sets['Outside IFN lists'])
for spec in plan['primary_comparisons']:
    source=spec['source']
    with np.load(OUT/'data'/('native_profiles_'+source+'.npz'),allow_pickle=False) as z:
        gp=scores.set_index('profile_id').loc[z['profile_ids']].reset_index();x=z['mean_log_expression'];detect=z['expressing_cells']/z['cells'][:,None];measured=z['common_measured'];detect[:,~measured]=np.nan
    selected=(gp.view==spec['view'])&(gp.analysis_celltype==spec['celltype'])&gp.unit_resolved&(gp.cells>=30)
    ref=references[references.source_label==source].set_index('human_ensembl').reindex(core.human_ensembl);assert ref.index.is_unique
    variable=ref.reference_variable.to_numpy();zz=np.full_like(x,np.nan);zz[:,variable]=(x[:,variable]-ref.reference_mean.to_numpy()[variable])/ref.reference_SD.to_numpy()[variable]*sgn[variable]
    ai=np.flatnonzero(selected&gp.infection_status.isin(spec['case']));bi=np.flatnonzero(selected&gp.infection_status.isin(spec['control']))
    if spec['design'] in ['within_donor_time_balanced','paired_donors']:
        fa=gp.loc[ai,['analysis_unit','time']].copy();fa['ix']=ai;fb=gp.loc[bi,['analysis_unit','time']].copy();fb['ix']=bi;paired=fa.merge(fb,on=['analysis_unit','time'],suffixes=('_case','_control'));rows=[]
        for unit,g in paired.groupby('analysis_unit'):
            rows.append((available_mean(zz[g.ix_case.to_numpy()]),available_mean(zz[g.ix_control.to_numpy()]),available_mean(detect[g.ix_case.to_numpy()]),available_mean(detect[g.ix_control.to_numpy()])))
    elif spec['design']=='paired_donors':
        aa=gp.loc[ai,['analysis_unit']].copy();aa['ix']=ai;bb=gp.loc[bi,['analysis_unit']].copy();bb['ix']=bi;rows=[]
        for unit in sorted(set(aa.analysis_unit)&set(bb.analysis_unit)):
            ia=aa.loc[aa.analysis_unit==unit,'ix'].to_numpy();ib=bb.loc[bb.analysis_unit==unit,'ix'].to_numpy();rows.append((available_mean(zz[ia]),available_mean(zz[ib]),available_mean(detect[ia]),available_mean(detect[ib])))
    else:
        case_ix=[g.index.to_numpy() for _,g in gp.loc[ai].groupby('analysis_unit')];control_ix=[g.index.to_numpy() for _,g in gp.loc[bi].groupby('analysis_unit')]
        rows=[(available_mean(np.array([available_mean(zz[v]) for v in case_ix])),available_mean(np.array([available_mean(zz[v]) for v in control_ix])),available_mean(np.array([available_mean(detect[v]) for v in case_ix])),available_mean(np.array([available_mean(detect[v]) for v in control_ix])))]
    if len(rows)<(3 if spec['design']!='unpaired_participants' else 1):continue
    casez=available_mean(np.array([r[0] for r in rows]));ctrlz=available_mean(np.array([r[1] for r in rows]));case_detect=available_mean(np.array([r[2] for r in rows]));ctrl_detect=available_mean(np.array([r[3] for r in rows]))
    totaln=int(variable.sum());dn=int(sum(g in degran and v for g,v in zip(core.human_ensembl,variable)))
    for j,r in enumerate(core.itertuples()):contrib.append(dict(comparison_id=spec['id'],source_label=source,cell_identity=spec['celltype'],human_ensembl=r.human_ensembl,display_symbol=r.display_symbol,frozen_direction=r.frozen_direction,feature_measured=bool(measured[j]),scored_variable=bool(variable[j]),signed_gene_effect=casez[j]-ctrlz[j],total_score_contribution=(casez[j]-ctrlz[j])/totaln if variable[j] else np.nan,in_degranulation=r.human_ensembl in degran,degranulation_score_contribution=(casez[j]-ctrlz[j])/dn if r.human_ensembl in degran and variable[j] else np.nan,outside_IFN_lists=r.human_ensembl in outside,case_detection_fraction=case_detect[j],control_detection_fraction=ctrl_detect[j]))
genes=pd.DataFrame(contrib);genes.to_csv(OUT/'results/fixed_gene_contributions.csv',index=False)
primary=pd.read_csv(OUT/'results/primary_cellular_effects.csv');checks=[]
for comparison,g in genes.groupby('comparison_id'):
    for target,col in [('CRRM388','total_score_contribution'),('Neutrophil degranulation','degranulation_score_contribution')]:
        expected=primary[(primary.comparison_id==comparison)&(primary.score==target)].estimate.iloc[0];actual=g[col].sum()
        if np.isfinite(expected):assert abs(expected-actual)<1e-10,(comparison,target,expected,actual);checks.append((comparison,target))
dump(OUT/'audit/cellular_interpretation_checks.json',dict(status='PASS',completed_at=datetime.now(timezone.utc).isoformat(),clinical_participants=36,fixed_composition_lineages=types,reference_fraction_sum=float(pref.sum()),minimum_observed_cell_coverage=float(patient.four_lineage_fraction.min()),patient_bootstraps=2000,gene_contribution_additivity_checks=len(checks),gene_members_changed=False))
print(pd.DataFrame(composition).query('endpoint == "RSV_vs_Healthy" and score in ["CRRM388","Outside IFN lists","Neutrophil degranulation"]')[['score','model','estimate','lower','upper']].to_string(index=False),flush=True)
