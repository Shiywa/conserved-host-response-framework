from pathlib import Path
from datetime import datetime,timezone
from itertools import combinations
from math import comb
import hashlib,json
import numpy as np,pandas as pd
from scipy import stats

OUT=(Path(__import__('os').environ['CRRM_WORKSPACE']) / 'Result4_analysis')
plan=json.loads((OUT/'parameters/analysis_plan_v1.json').read_text());receipt=json.loads((OUT/'parameters/declaration_receipt.json').read_text())
assert receipt['analysis_plan_sha256']==hashlib.sha256((OUT/'parameters/analysis_plan_v1.json').read_bytes()).hexdigest()
assert json.loads((OUT/'audit/expression_extraction_receipt.json').read_text())['status']=='PASS'
meta=pd.read_csv(OUT/'data/profile_metadata.csv');core=pd.read_csv(OUT/'parameters/frozen_CRRM388.csv');sets=json.loads((OUT/'parameters/score_gene_sets.json').read_text());branches=pd.read_csv(OUT/'parameters/frozen_functional_branches.csv').branch.drop_duplicates().tolist()
sign=np.where(core.frozen_direction=='up',1.,-1.);geneid=core.human_ensembl.to_numpy();score_frames=[];measurement=[];references=[]
for label in sorted(meta.source_label.unique()):
    with np.load(OUT/'data'/('native_profiles_'+label+'.npz'),allow_pickle=False) as z:
        gp=meta.set_index('profile_id').loc[z['profile_ids']].reset_index();x=z['mean_log_expression'];measured=z['common_measured'];depth=z['mean_log_host_depth']
        for study,ix in gp.groupby('study_id',sort=False).groups.items():
            ix=np.array(list(ix));gg=gp.loc[ix].copy();xx=x[ix];eligible=(gg.view=='canonical')&(gg.cells>=30)
            if label in ['Berg_ALI','CD1c_DC','Infant_PBMC']:eligible&=gg.unit_resolved
            elif (eligible&gg.unit_resolved).any():eligible&=gg.unit_resolved
            if label=='Infant_PBMC':eligible&=gg.infection_status.isin(['Healthy','Mild RSV','Moderate RSV','Severe RSV'])
            rr=gg.loc[eligible];weights=1./rr.groupby('analysis_unit').analysis_unit.transform('size').to_numpy() if len(rr) else np.array([])
            if len(rr)>=2:
                weights/=weights.sum();ref=xx[eligible.to_numpy()];ref0=np.nan_to_num(ref,nan=0.);mu=np.sum(weights[:,None]*ref0,axis=0);sd=np.sqrt(np.sum(weights[:,None]*(ref0-mu)**2,axis=0));variable=measured&np.isfinite(sd)&(sd>1e-8)
            else:mu=np.zeros(388);sd=np.zeros(388);variable=np.zeros(388,dtype=bool)
            zz=np.full_like(xx,np.nan);zz[:,variable]=(xx[:,variable]-mu[variable])/sd[variable];signed=zz*sign
            for j,g in enumerate(geneid):references.append(dict(source_label=label,study_id=study,human_ensembl=g,display_symbol=core.iloc[j].display_symbol,measured=bool(measured[j]),reference_variable=bool(variable[j]),reference_mean=mu[j],reference_SD=sd[j],reference_profiles=len(rr),reference_units=rr.analysis_unit.nunique()))
            for target,genes in sets.items():
                mask=np.isin(geneid,genes)&variable;n=int(mask.sum());nominal=len(genes)
                if target=='CRRM388':eligible_target=n>=np.ceil(.7*388) and int((mask&(sign>0)).sum())>=140 and int((mask&(sign<0)).sum())>=54
                else:eligible_target=n>=max(5,np.ceil(.5*nominal))
                gg[target]=np.mean(signed[:,mask],axis=1) if eligible_target else np.nan
                measurement.append(dict(source_label=label,study_id=study,score=target,nominal_members=nominal,measured_members=int((np.isin(geneid,genes)&measured).sum()),scored_variable_members=n,eligible=bool(eligible_target),reference_profiles=len(rr),reference_units=rr.analysis_unit.nunique()))
            nu=int((variable&(sign>0)).sum());nd=int((variable&(sign<0)).sum());total=nu+nd
            gg['Up contribution']=gg.Up*nu/total if total else np.nan;gg['Down contribution']=gg.Down*nd/total if total else np.nan
            ok=gg.CRRM388.notna();assert np.allclose((gg.loc[ok,'Up contribution']+gg.loc[ok,'Down contribution']),gg.loc[ok,'CRRM388'],atol=1e-12),label
            gg['mean_log_host_depth']=depth[ix];score_frames.append(gg)
scores=pd.concat(score_frames,ignore_index=True).sort_values('profile_id');scores.to_csv(OUT/'results/profile_scores.csv',index=False);pd.DataFrame(measurement).to_csv(OUT/'results/score_measurement_coverage.csv',index=False);pd.DataFrame(references).to_csv(OUT/'data/gene_reference_parameters.csv',index=False)

permutation_cache={}
def unpaired_allocations(n,na):
    key=(n,na)
    if key not in permutation_cache:
        total=comb(n,na);exact=total<=100000
        indices=np.array(list(combinations(range(n),na))) if exact else np.argsort(np.random.default_rng(20261113).random((19999,n)),axis=1)[:,:na]
        w=np.zeros((len(indices),n));w[np.arange(len(indices))[:,None],indices]=1.;permutation_cache[key]=(w,exact)
    return permutation_cache[key]
def infer(a,b=None,paired=False):
    a=np.asarray(a,dtype=float);b=np.asarray(b,dtype=float) if b is not None else None
    if paired:
        delta=a-b if b is not None else a;n=len(delta)
        if n<3:return dict(status='INSUFFICIENT_UNITS',estimate=np.nan,lower=np.nan,upper=np.nan,P=np.nan,resampling_P=np.nan)
        estimate=delta.mean();se=stats.sem(delta);df=n-1;p=stats.ttest_1samp(delta,0).pvalue
        signs=((np.arange(2**n)[:,None]>>np.arange(n))&1)*2-1 if n<=15 else np.random.default_rng(20261113).choice([-1,1],size=(19999,n))
        perm=signs@delta/n;resample=float(np.mean(np.abs(perm)>=abs(estimate)-1e-12)) if n<=15 else float((1+np.sum(np.abs(perm)>=abs(estimate)-1e-12))/(len(perm)+1))
    else:
        n,m=len(a),len(b)
        if min(n,m)<3:return dict(status='INSUFFICIENT_UNITS',estimate=np.nan,lower=np.nan,upper=np.nan,P=np.nan,resampling_P=np.nan)
        va=a.var(ddof=1)/n;vb=b.var(ddof=1)/m;se=np.sqrt(va+vb);df=(va+vb)**2/(va**2/(n-1)+vb**2/(m-1));estimate=a.mean()-b.mean();p=stats.ttest_ind(a,b,equal_var=False).pvalue
        y=np.r_[a,b];w,exact=unpaired_allocations(n+m,n);ma=w@y/n;mb=(y.sum()-w@y)/m;sa=(w@(y*y)-n*ma*ma)/(n-1);sb=((y*y).sum()-w@(y*y)-m*mb*mb)/(m-1);t=(ma-mb)/np.sqrt(np.maximum(sa/n+sb/m,1e-20));obs=estimate/max(se,1e-12);count=np.sum(np.abs(t)>=abs(obs)-1e-10);resample=float(count/len(t)) if exact else float((count+1)/(len(t)+1))
    ci=stats.t.ppf(.975,df)*se
    return dict(status='OK',estimate=float(estimate),lower=float(estimate-ci),upper=float(estimate+ci),P=float(p),resampling_P=resample)

def contrast_rows(spec,target,min_cells=30):
    g=scores[(scores.source_label==spec['source'])&(scores.view==spec['view'])&(scores.analysis_celltype==spec['celltype'])&scores.unit_resolved&(scores.cells>=min_cells)].dropna(subset=[target])
    a=g[g.infection_status.isin(spec['case'])];b=g[g.infection_status.isin(spec['control'])]
    if spec['design'] in ['within_donor_time_balanced','paired_donors']:
        aa=a.groupby(['analysis_unit','time'])[target].mean();bb=b.groupby(['analysis_unit','time'])[target].mean();pair=pd.concat([aa.rename('case'),bb.rename('control')],axis=1).dropna();units=pair.groupby(level=0)[['case','control']].mean();timecounts=pair.groupby(level=0).size()
    else:
        aa=a.groupby('analysis_unit')[target].mean();bb=b.groupby('analysis_unit')[target].mean();units=pd.concat([aa.rename('case'),bb.rename('control')],axis=1);timecounts=pd.Series(dtype=int)
    paired=spec['design']!='unpaired_participants'
    if paired:units=units.dropna();result=infer(units['case'],units['control'],paired=True);ncase=ncontrol=len(units)
    else:result=infer(units['case'].dropna(),units['control'].dropna());ncase=int(units['case'].notna().sum());ncontrol=int(units['control'].notna().sum())
    row=dict(comparison_id=spec['id'],source_label=spec['source'],view=spec['view'],cell_identity=spec['celltype'],score=target,comparison=' + '.join(spec['case'])+' minus '+' + '.join(spec['control']),design=spec['design'],minimum_cells=min_cells,case_units=ncase,control_units=ncontrol,**result)
    unitrows=[dict(comparison_id=spec['id'],score=target,analysis_unit=u,case_score=r['case'],control_score=r['control'],paired_difference=r['case']-r['control'] if paired else np.nan,matched_time_points=int(timecounts.get(u,0))) for u,r in units.iterrows()]
    return row,unitrows

def bh(values):
    p=np.array(values,dtype=float);valid=np.isfinite(p);p0=np.where(valid,p,1.);order=np.argsort(p0);q=np.minimum.accumulate((p0[order]*len(p0)/np.arange(1,len(p0)+1))[::-1])[::-1];out=np.empty(len(p0));out[order]=np.minimum(q,1);out[~valid]=np.nan;return out

out=[];unitdata=[]
targets=list(sets)+['Up contribution','Down contribution']
for spec in plan['primary_comparisons']:
    for target in targets:
        row,units=contrast_rows(spec,target);row['scope']='Primary contrast';out.append(row);unitdata+=units
    print('Primary contrast completed',spec['id'],flush=True)
results=pd.DataFrame(out)
families={'Primary total (BH8)':['CRRM388'],'Fixed branches (BH56)':branches,'IFN-list components (BH16)':['IFN-list overlap','Outside IFN lists'],'LCN2 omission (BH8)':['Degranulation without LCN2'],'Up/down contributions (BH16)':['Up contribution','Down contribution'],'Direction means (BH16)':['Up','Down']}
for name,targets_ in families.items():
    mask=results.score.isin(targets_);results.loc[mask,'correction_family']=name;results.loc[mask,'adjusted_P']=bh(results.loc[mask,'P']);results.loc[mask,'resampling_adjusted_P']=bh(results.loc[mask,'resampling_P'])
results.to_csv(OUT/'results/primary_cellular_effects.csv',index=False);pd.DataFrame(unitdata).to_csv(OUT/'results/primary_unit_effects.csv',index=False)

secondary=[];secondary_units=[]
for spec in plan['primary_comparisons']:
    if spec['source']=='Berg_ALI':
        changes=[('bystander_vs_control',['bystdr'],['ctrl']),('infected_vs_bystander',['inf'],['bystdr'])]
    elif spec['source']=='Infant_PBMC':changes=[('Severe_vs_Mild',['Severe RSV'],['Mild RSV'])]
    else:continue
    for suffix,case,control in changes:
        ss=dict(spec,id=spec['id']+'_'+suffix,case=case,control=control)
        for target in ['CRRM388','Outside IFN lists','Degranulation without LCN2']+branches:
            row,units=contrast_rows(ss,target);row['scope']=suffix;secondary.append(row);secondary_units+=units
sec=pd.DataFrame(secondary)
for (scope,targetkind),ix in sec.groupby(['scope',sec.score.map(lambda x:'branches' if x in branches else x)]).groups.items():sec.loc[ix,'correction_family']=str(scope)+' | '+str(targetkind);sec.loc[ix,'adjusted_P']=bh(sec.loc[ix,'P']);sec.loc[ix,'resampling_adjusted_P']=bh(sec.loc[ix,'resampling_P'])
sec.to_csv(OUT/'results/secondary_cellular_effects.csv',index=False);pd.DataFrame(secondary_units).to_csv(OUT/'results/secondary_unit_effects.csv',index=False)

qc={'status':'PASS','completed_at':datetime.now(timezone.utc).isoformat(),'physical_human_cells':2119831,'source_labels':13,'profile_rows':len(scores),'primary_total_tests':8,'primary_branch_tests':56,'frozen_gene_members':388,'new_score_gene_selection':False,'reference_gene_SD_guard':1e-8,'primary_total_results':results[results.score=='CRRM388'][['comparison_id','case_units','control_units','estimate','P','adjusted_P','resampling_P','status']].to_dict('records')}
(OUT/'audit/scoring_and_primary_tests.json').write_text(json.dumps(qc,indent=2)+'\n');print(json.dumps(qc,indent=2),flush=True)
