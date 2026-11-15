"""Frozen CRRM388 scoring, whole-unit inference, and prespecified clinical bridge."""
from pathlib import Path
from datetime import datetime, timezone
from itertools import combinations
import hashlib, json, math
import numpy as np
import pandas as pd
from scipy import stats

OUT=(Path(__import__('os').environ['CRRM_WORKSPACE']) / 'Result6_analysis')
REV=Path(str(Path(__import__('os').environ['CRRM_WORKSPACE']) / ''))
R3=REV/'phase20_Result3_clinical_analysis_20261112'
def dump(path,value): path.write_text(json.dumps(value,indent=2,ensure_ascii=False,default=lambda a:a.item() if isinstance(a,np.generic) else str(a))+'\n')
plan=json.loads((OUT/'parameters/analysis_plan_v1.json').read_text())
assert hashlib.sha256((OUT/'parameters/analysis_plan_v1.json').read_bytes()).hexdigest()==json.loads((OUT/'parameters/declaration_receipt.json').read_text())['plan_sha256']
core=pd.read_csv(OUT/'parameters/frozen_CRRM388.csv');sets=json.loads((OUT/'parameters/score_gene_sets.json').read_text())
genes=core.human_ensembl.to_numpy();sign=np.where(core.frozen_direction.eq('up'),1.,-1.)
branches=pd.read_csv(OUT/'parameters/frozen_functional_branches.csv').branch.drop_duplicates().tolist()
for target in ['TNF / NF-kB','Neutrophil degranulation']:
    sets[target+' outside IFN lists']=sorted(set(sets[target])-set(sets['IFN-list overlap']))
for target in ['TNF / NF-kB','Neutrophil degranulation']:
    sets[target+' minus IFN-list overlap']=[]  # Computed as within-profile score differences below.
contrast_targets=[t+' minus IFN-list overlap' for t in ['TNF / NF-kB','Neutrophil degranulation']]
disjoint_targets=[t+' outside IFN lists minus IFN-list overlap' for t in ['TNF / NF-kB','Neutrophil degranulation']]

def bh(p,n=None):
    p=np.array(p,float);valid=np.isfinite(p);n=len(p) if n is None else n
    q=np.where(valid,p,1.);ix=np.argsort(q);v=np.minimum.accumulate((q[ix]*n/np.arange(1,len(q)+1))[::-1])[::-1]
    out=np.empty(len(q));out[ix]=np.minimum(v,1);out[~valid]=np.nan;return out

def infer(a,b,exact=True):
    a=np.asarray(a,float);b=np.asarray(b,float);a=a[np.isfinite(a)];b=b[np.isfinite(b)]
    n,m=len(a),len(b)
    if min(n,m)<3:return dict(estimate=np.nan,lower=np.nan,upper=np.nan,SE=np.nan,df=np.nan,P=np.nan,exact_P=np.nan,n_case=n,n_control=m,status='INSUFFICIENT_UNITS')
    d=float(a.mean()-b.mean())
    if np.ptp(a)<=1e-12 and np.ptp(b)<=1e-12:
        return dict(estimate=d,lower=np.nan,upper=np.nan,SE=0.,df=np.nan,P=1.,exact_P=1.,n_case=n,n_control=m,status='BOTH_GROUPS_CONSTANT')
    va=a.var(ddof=1)/n;vb=b.var(ddof=1)/m;se=np.sqrt(va+vb);df=(va+vb)**2/(va**2/(n-1)+vb**2/(m-1));q=stats.t.ppf(.975,df)
    p=float(2*stats.t.sf(abs(d/se),df));pe=np.nan
    if exact:
        y=np.r_[a,b];assign=np.array(list(combinations(range(n+m),n)))
        a0=y[assign];sa=a0.sum(1);ss=(a0*a0).sum(1);ma=sa/n;mb=(y.sum()-sa)/m
        v1=np.maximum((ss-n*ma*ma)/(n-1),0);v2=np.maximum(((y*y).sum()-ss-m*mb*mb)/(m-1),0)
        null=(ma-mb)/np.sqrt(np.maximum(v1/n+v2/m,1e-30))
        pe=float(np.mean(abs(null)>=abs(d/se)-1e-10))
    return dict(estimate=d,lower=float(d-q*se),upper=float(d+q*se),SE=float(se),df=float(df),P=p,exact_P=pe,n_case=n,n_control=m,status='OK')

coverage=[];references=[];data={}
def score(m,values,measured,source,method):
    z=m.copy()
    for target, members in sets.items():
        if not members:continue
        mask=measured&np.isin(genes,members);n=int(mask.sum());nominal=len(members)
        eligible=n>=max(5,math.ceil(.5*nominal))
        if target=='CRRM388':eligible=n>=math.ceil(.7*388) and (mask&(sign>0)).sum()>=140 and (mask&(sign<0)).sum()>=54
        z[target]=values[:,mask].mean(1) if eligible else np.nan
        coverage.append(dict(source=source,method=method,score=target,nominal_genes=nominal,scored_genes=n,eligible=bool(eligible)))
    up=measured&(sign>0);down=measured&(sign<0);n=up.sum()+down.sum()
    z['Up contribution']=values[:,up].mean(1)*up.sum()/n
    z['Down contribution']=values[:,down].mean(1)*down.sum()/n
    assert np.allclose(z.CRRM388,z['Up contribution']+z['Down contribution'],equal_nan=True)
    for target in ['TNF / NF-kB','Neutrophil degranulation']:
        z[target+' minus IFN-list overlap']=z[target]-z['IFN-list overlap']
        z[target+' outside IFN lists minus IFN-list overlap']=z[target+' outside IFN lists']-z['IFN-list overlap']
    z['method']=method
    return z

def standardize(x,measured,weights=None):
    if weights is None:weights=np.full(len(x),1/len(x));ddof=True
    else:ddof=False
    x0=np.nan_to_num(x,nan=0.)
    mu=(x0*weights[:,None]).sum(0)
    sd=x.std(0,ddof=1) if ddof else np.sqrt((((x0-mu)**2)*weights[:,None]).sum(0))
    variable=measured&np.isfinite(sd)&(sd>1e-8)
    z=np.full_like(x,np.nan);z[:,variable]=(x[:,variable]-mu[variable])/sd[variable]*sign[variable]
    return z,variable,mu,sd

hm=pd.read_csv(OUT/'data/human_sample_manifest.csv');hm['cells']=np.nan
h=pd.read_csv(OUT/'data/human_author_TPM.csv.gz').set_index('gene_id')
full=h[hm.column].to_numpy().T;native=h.index.to_numpy();lookup={g:i for i,g in enumerate(native)}
measured=np.isin(genes,native);xx=np.full((len(hm),388),np.nan);xx[:,measured]=np.log2(1+full[:,[lookup[g] for g in genes[measured]]])
signed,variable,mu,sd=standardize(xx,measured)
hs=score(hm,signed,variable,'BRD4_hSAEC','Primary')
rr=(stats.rankdata(full,axis=1)-1)/(full.shape[1]-1)-.5
rank=np.full_like(xx,np.nan);rank[:,measured]=rr[:,[lookup[g] for g in genes[measured]]]*sign[measured]
hr=score(hm,rank,measured,'BRD4_hSAEC','Rank')
scaled=full/full.sum(1)[:,None]*1e6;scaledx=np.full_like(xx,np.nan);scaledx[:,measured]=np.log2(1+scaled[:,[lookup[g] for g in genes[measured]]])
zs,vs,_,_=standardize(scaledx,measured)
ha=score(hm,zs,vs,'BRD4_hSAEC','Abundance rescaling')
data['BRD4_hSAEC']=dict(meta=hm,expression=xx,measured=measured,signed=signed,variable=variable,rankgene=rank,scores={'Primary':hs,'Rank':hr,'Abundance rescaling':ha})
for j,g in enumerate(genes):references.append(dict(source='BRD4_hSAEC',human_ensembl=g,reference_mean=mu[j],reference_SD=sd[j],measured=bool(measured[j]),variable=bool(variable[j])))
np.savez_compressed(OUT/'data/human_score_arrays.npz',human_ensembl=genes.astype(str),unit=hm.unit.to_numpy(dtype=str),expression=xx,signed=signed,variable=variable,rank=rank,measured=measured)

mm=pd.read_csv(OUT/'data/mouse_profile_metadata.csv');projection=pd.read_csv(OUT/'parameters/frozen_mouse_member_projection.csv')
with np.load(OUT/'data/mouse_expression_profiles.npz',allow_pickle=False) as z:
    assert np.array_equal(mm.profile_id.to_numpy(),z['profile_ids'])
    mx=z['mean_log_expression'];full=z['transcriptome_mean_log_expression'];measured=z['common_measured'];native=z['native_gene_ids'];detection=z['expressing_cells']/z['cells'][:,None]
eligible=mm.celltype.isin(plan['primary_lineages'])&mm.cells.ge(30)
ref=mm[eligible];weights=1/ref.groupby('unit').unit.transform('size').to_numpy();weights/=weights.sum()
zs,variable,mu,sd=standardize(mx[eligible],measured,weights)
signed=np.full_like(mx,np.nan);signed[:,variable]=(mx[:,variable]-mu[variable])/sd[variable]*sign[variable]
ms=score(mm,signed,variable,'RelA_mouse','Primary')
nr=(stats.rankdata(full,axis=1)-1)/(full.shape[1]-1)-.5;index={g:i for i,g in enumerate(native)}
rank=np.full_like(mx,np.nan);rank[:,measured]=nr[:,[index[g] for g in projection.loc[measured,'ensembl']]]*sign[measured]
mr=score(mm,rank,measured,'RelA_mouse','Rank')
supported=np.zeros(388,bool)
for ct in plan['primary_lineages']:
    use=eligible&mm.celltype.eq(ct)
    supported|=(detection[use]>=.1).sum(0)>=3
supported&=variable
md=score(mm,signed,supported,'RelA_mouse','Detection-supported')
data['RelA_mouse']=dict(meta=mm,expression=mx,measured=measured,signed=signed,variable=variable,rankgene=rank,detection=detection,scores={'Primary':ms,'Rank':mr,'Detection-supported':md})
for j,g in enumerate(genes):references.append(dict(source='RelA_mouse',human_ensembl=g,reference_mean=mu[j],reference_SD=sd[j],measured=bool(measured[j]),variable=bool(variable[j]),detection_supported=bool(supported[j])))
np.savez_compressed(OUT/'data/mouse_score_arrays.npz',human_ensembl=genes.astype(str),profile_id=mm.profile_id.to_numpy(),signed=signed,variable=variable,rank=rank,measured=measured)
pd.DataFrame(coverage).to_csv(OUT/'results/measurement_coverage.csv',index=False)
pd.DataFrame(references).to_csv(OUT/'data/gene_reference_parameters.csv',index=False)
pd.concat([s for d in data.values() for s in d['scores'].values()],ignore_index=True).to_csv(OUT/'results/profile_scores.csv',index=False)

targets=['CRRM388']+branches+['IFN-list overlap','Outside IFN lists','Up contribution','Down contribution','Degranulation without LCN2']+contrast_targets+disjoint_targets
def target_family(t):
    if t=='CRRM388':return 'Total',12
    if t in branches:return 'Fixed branches',84
    if t in ['IFN-list overlap','Outside IFN lists']:return 'IFN components',24
    if t in ['Up contribution','Down contribution']:return 'Direction contributions',24
    if t in contrast_targets:return 'Direct component differences',24
    if t in disjoint_targets:return 'Disjoint component differences',24
    return 'LCN2 omission',12

effects=[];unitrows=[];gene_effects=[];gene_units=[];omissions=[]
for spec in plan['comparisons']:
    d=data[spec['source']];m=d['meta']
    for method,gp in d['scores'].items():
        for scope,case,control in [('Intervention',spec['case'],spec['control']),('Infection context',spec['context_case'],spec['context_control'])]:
            use=gp.celltype.eq(spec['celltype'])
            if spec['source']=='RelA_mouse':use&=gp.cells.ge(30)
            a=gp[use&gp.group.eq(case)];b=gp[use&gp.group.eq(control)]
            for target in targets:
                family,nfam=target_family(target)
                key='|'.join([scope,method,spec['key'],target])
                row=dict(test_id=key,key=spec['key'],source=spec['source'],celltype=spec['celltype'],method=method,scope=scope,case=case,control=control,score=target,family='|'.join([scope,method,family]),family_size=nfam,**infer(a[target],b[target]))
                effects.append(row)
                for role,g in [('case',a),('control',b)]:
                    for u,v in zip(g.unit,g[target]):unitrows.append(dict(test_id=key,role=role,unit=u,value=v))
    gp=d['scores']['Primary'];a=gp.celltype.eq(spec['celltype'])&gp.group.eq(spec['case']);b=gp.celltype.eq(spec['celltype'])&gp.group.eq(spec['control'])
    if spec['source']=='RelA_mouse':a&=gp.cells.ge(30);b&=gp.cells.ge(30)
    for j,g in enumerate(genes):
        values=d['signed'][:,j];gx=d['expression'][:,j]
        r=infer(values[a],values[b],exact=False)
        key=spec['key']+'|'+g
        gene_effects.append(dict(test_id=key,key=spec['key'],source=spec['source'],celltype=spec['celltype'],human_ensembl=g,display_symbol=core.iloc[j].display_symbol,frozen_direction=core.iloc[j].frozen_direction,log_expression_effect=float(np.nanmean(gx[a])-np.nanmean(gx[b])) if d['measured'][j] else np.nan,**r))
        if d['variable'][j]:
            for role,mask in [('case',a),('control',b)]:
                for u,v in zip(m.loc[mask,'unit'],values[mask]):gene_units.append(dict(test_id=key,role=role,unit=u,value=v))
    for target in ['CRRM388','Neutrophil degranulation']:
        mask=d['variable']&np.isin(genes,sets[target]);ids=np.flatnonzero(mask);delta=np.nanmean(d['signed'][a],0)-np.nanmean(d['signed'][b],0)
        total=float(delta[mask].sum())
        for j in ids:omissions.append(dict(key=spec['key'],score=target,omitted_gene=core.iloc[j].display_symbol,human_ensembl=genes[j],estimate_without_gene=float((total-delta[j])/(len(ids)-1)),full_estimate=total/len(ids),members_after_omission=len(ids)-1))
    print('Completed',spec['key'],flush=True)

for spec in plan['comparisons']:
    for method in ['Primary','Rank','Abundance rescaling','Detection-supported']:
        if method in data[spec['source']]['scores']:continue
        for scope,case,control in [('Intervention',spec['case'],spec['control']),('Infection context',spec['context_case'],spec['context_control'])]:
            for target in targets:
                family,nfam=target_family(target)
                effects.append(dict(test_id='|'.join([scope,method,spec['key'],target]),key=spec['key'],source=spec['source'],celltype=spec['celltype'],method=method,scope=scope,case=case,control=control,score=target,family='|'.join([scope,method,family]),family_size=nfam,estimate=np.nan,lower=np.nan,upper=np.nan,SE=np.nan,df=np.nan,P=np.nan,exact_P=np.nan,n_case=0,n_control=0,status='METHOD_INAPPLICABLE'))
ef=pd.DataFrame(effects)
for family,ix in ef.groupby('family').groups.items():
    n=int(ef.loc[ix,'family_size'].iloc[0]);assert len(ix)==n
    ef.loc[ix,'adjusted_P']=bh(ef.loc[ix,'P'],n);ef.loc[ix,'exact_adjusted_P']=bh(ef.loc[ix,'exact_P'],n)
ef.to_csv(OUT/'results/score_effects.csv',index=False)
pd.DataFrame(unitrows).to_csv(OUT/'data/inference_units.csv',index=False)
ge=pd.DataFrame(gene_effects);assert len(ge)==388*12
ge['adjusted_P_BH4656']=bh(ge.P,388*12);ge.to_csv(OUT/'results/intervention_gene_effects.csv',index=False)
pd.DataFrame(gene_units).to_csv(OUT/'data/gene_inference_units.csv.gz',index=False,compression='gzip')
pd.DataFrame(omissions).to_csv(OUT/'results/every_gene_omission.csv',index=False)

# Thresholds change eligibility only; retain primary reference and target families.
thresholds=[]
for threshold in [15,50]:
    for spec in plan['comparisons']:
        gp=data[spec['source']]['scores']['Primary'];use=gp.celltype.eq(spec['celltype'])
        if spec['source']=='RelA_mouse':use&=gp.cells.ge(threshold)
        a=gp[use&gp.group.eq(spec['case'])];b=gp[use&gp.group.eq(spec['control'])]
        for target in ['CRRM388']+branches+['IFN-list overlap','Outside IFN lists']:
            family,n=target_family(target)
            thresholds.append(dict(key=spec['key'],source=spec['source'],celltype=spec['celltype'],method=f'Minimum{threshold}cells',minimum_cells=threshold,score=target,family=family,family_size=n,**infer(a[target],b[target])))
th=pd.DataFrame(thresholds)
for (threshold,family),ix in th.groupby(['minimum_cells','family']).groups.items():th.loc[ix,'adjusted_P']=bh(th.loc[ix,'P'],int(th.loc[ix,'family_size'].iloc[0]))
th.to_csv(OUT/'results/threshold_sensitivities.csv',index=False)

# Human intervention-unit permutations preserve each complete correlated gene vector.
cg=pd.read_csv(R3/'tables/12_clinical_gene_effects.csv')
degran=sets['Neutrophil degranulation'];assert len(degran)==25
pairs=pd.read_csv(R3/'tables/03_verified_recovery_pairs.csv');assert len(pairs)==15
cx=pd.read_csv(R3/'data/processed/GSE69606_log_expression.csv.gz',index_col=0)
recovery=(cx[pairs.recovery_sample].to_numpy()-cx[pairs.acute_sample].to_numpy()).mean(1)
rec=pd.DataFrame(dict(human_ensembl=cx.index,log_expression_effect=recovery)).set_index('human_ensembl')
human=data['BRD4_hSAEC'];hm=human['meta'];intervention=hm.group.isin(['RSV','RSV + BRD4i']);x=human['expression'][intervention]
group=hm.loc[intervention,'group'].to_numpy();case=group=='RSV + BRD4i';aids=np.array(list(combinations(range(8),4)))
null=np.array([x[ii].mean(0)-x[np.setdiff1d(np.arange(8),ii)].mean(0) for ii in aids]);effect=x[case].mean(0)-x[~case].mean(0)
clinical_vectors={key:cg[cg.contrast_key.eq(key)].set_index('human_ensembl') for key in ['GSE246622_Severe_Mild','GSE69606_Severe_Mild']}
clinical_vectors['GSE69606_Recovery_Acute']=rec
alignment=[];points=[]
for key,cv in clinical_vectors.items():
    common=[g for g in degran if g in cv.index and human['variable'][np.where(genes==g)[0][0]] and np.isfinite(cv.loc[g,'log_expression_effect'])]
    ix=np.array([np.where(genes==g)[0][0] for g in common]);yy=cv.loc[common,'log_expression_effect'].to_numpy();xx=effect[ix]
    r=float(stats.spearmanr(xx,yy).statistic)
    rv=np.array([stats.spearmanr(n[ix],yy).statistic for n in null]);pe=float(np.mean(np.abs(rv)>=abs(r)-1e-12))
    direction='same as recovery' if 'Recovery' in key else 'opposite to severity'
    desired=(xx*yy>0) if 'Recovery' in key else (xx*yy<0)
    alignment.append(dict(clinical_contrast=key,members=len(ix),rho=r,permutation_P=pe,allocations=70,expected_alignment=direction,directionally_aligned=int(desired.sum()),nonzero_pairs=int(((xx!=0)&(yy!=0)).sum())))
    for j,g,xc,yc,yes in zip(ix,common,xx,yy,desired):points.append(dict(clinical_contrast=key,human_ensembl=g,display_symbol=core.iloc[j].display_symbol,intervention_log_effect=xc,clinical_log_effect=yc,aligned=bool(yes)))
al=pd.DataFrame(alignment);al['adjusted_P_BH3']=bh(al.permutation_P,3)
al.to_csv(OUT/'results/clinical_gene_alignment.csv',index=False);pd.DataFrame(points).to_csv(OUT/'results/clinical_gene_alignment_points.csv',index=False)
pd.DataFrame(dict(human_ensembl=cx.index,recovery_minus_acute_log_effect=recovery,n_pairs=15)).to_csv(OUT/'results/frozen_recovery_gene_vector.csv',index=False)

primary=ef[ef.method.eq('Primary')&ef.scope.eq('Intervention')&ef.score.eq('CRRM388')]
receipt=dict(status='PASS',completed_at=datetime.now(timezone.utc).isoformat(),frozen_members=388,primary_contrasts=12,score_tests=len(ef),gene_tests=len(ge),clinical_alignment_tests=3,mouse_cells=84967,mouse_units=11,primary_total=primary[['key','estimate','lower','upper','P','adjusted_P','exact_P','n_case','n_control']].to_dict('records'),clinical_alignment=al.to_dict('records'))
dump(OUT/'audit/scoring_receipt.json',receipt)
print('PRIMARY TOTAL\n'+primary[['key','estimate','lower','upper','adjusted_P','exact_P']].to_string(index=False))
print('HUMAN COMPONENTS\n'+ef[ef.key.eq('BRD4_Airway')&ef.method.eq('Primary')&ef.scope.eq('Intervention')][['score','estimate','adjusted_P','exact_P']].to_string(index=False))
print('CLINICAL ALIGNMENT\n'+al.to_string(index=False))
