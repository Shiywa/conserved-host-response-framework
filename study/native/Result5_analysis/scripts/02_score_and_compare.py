from pathlib import Path
from datetime import datetime,timezone
from itertools import product
import json,hashlib,math
import numpy as np,pandas as pd
from scipy import stats

OUT=(Path(__import__('os').environ['CRRM_WORKSPACE']) / 'Result5_analysis');R4=(OUT.parent/'Result4') if (OUT.parent/'Result4').exists() else (OUT.parent/'Result4_analysis');RESOURCE=Path(__import__('os').environ.get('CRRM_RESOURCE_ROOT','external_resource'))
plan=json.loads((OUT/'parameters/analysis_plan_v1.json').read_text());receipt=json.loads((OUT/'parameters/declaration_receipt.json').read_text())
assert hashlib.sha256((OUT/'parameters/analysis_plan_v1.json').read_bytes()).hexdigest()==receipt['plan_sha256']
core=pd.read_csv(OUT/'parameters/frozen_CRRM388.csv');sets=json.loads((OUT/'parameters/score_gene_sets.json').read_text());branches=pd.read_csv(OUT/'parameters/frozen_functional_branches.csv').branch.drop_duplicates().tolist()
sign=np.where(core.frozen_direction.eq('up'),1.,-1.);genes=core.human_ensembl.to_numpy();rng=np.random.default_rng(20261113)
targets=['CRRM388','Up contribution','Down contribution','IFN-list overlap','Outside IFN lists']+branches
def dump(p,x):p.write_text(json.dumps(x,indent=2,ensure_ascii=False,default=lambda x:x.item() if isinstance(x,np.generic) else str(x))+'\n')
def bh(p,n=None):
    p=np.array(p,float);valid=np.isfinite(p);n=n or len(p);pp=np.where(valid,p,1.);order=np.argsort(pp);v=np.minimum.accumulate((pp[order]*n/np.arange(1,len(p)+1))[::-1])[::-1];out=np.empty(len(p));out[order]=np.minimum(v,1);out[~valid]=np.nan;return out
def infer(a,b=None,paired=False):
    a=np.asarray(a,float);b=np.asarray(b,float) if b is not None else None
    if paired:
        d=a-b if b is not None else a;d=d[np.isfinite(d)];n=len(d)
        if n<3:return dict(estimate=np.nan,lower=np.nan,upper=np.nan,P=np.nan,df=n-1)
        eff=d.mean();se=stats.sem(d);df=n-1
    else:
        a=a[np.isfinite(a)];b=b[np.isfinite(b)];n,m=len(a),len(b)
        if min(n,m)<3:return dict(estimate=np.nan,lower=np.nan,upper=np.nan,P=np.nan,df=np.nan)
        va=a.var(ddof=1)/n;vb=b.var(ddof=1)/m;se=np.sqrt(va+vb);df=(va+vb)**2/(va*va/(n-1)+vb*vb/(m-1)) if va+vb>1e-30 else n+m-2;eff=a.mean()-b.mean()
    p=2*stats.t.sf(abs(eff)/se,df) if se>1e-15 else float(abs(eff)<1e-12);q=stats.t.ppf(.975,df)
    return dict(estimate=float(eff),lower=float(eff-q*se),upper=float(eff+q*se),P=float(p),df=float(df))
def hc3(y,group,age,sex):
    X=np.column_stack([np.ones(len(y)),group,np.log1p(age),np.asarray(sex)=='M']).astype(float);y=np.asarray(y,float)
    if not np.isfinite(X).all() or np.linalg.matrix_rank(X)<4:return dict(estimate=np.nan,lower=np.nan,upper=np.nan,P=np.nan,df=np.nan)
    inv=np.linalg.inv(X.T@X);coef=inv@X.T@y;res=y-X@coef;h=np.sum((X@inv)*X,axis=1);w=res**2/(1-h)**2;cov=inv@(X.T@(w[:,None]*X))@inv;se=np.sqrt(cov[1,1]);df=len(y)-4;q=stats.t.ppf(.975,df)
    return dict(estimate=float(coef[1]),lower=float(coef[1]-q*se),upper=float(coef[1]+q*se),P=float(2*stats.t.sf(abs(coef[1]/se),df)),df=df)
def permutations(a,b=None,paired=False):
    if paired:
        d=np.asarray(a)-np.asarray(b);n=len(d);s=np.array(list(product([-1,1],repeat=n)));null=s@d/n;return float(np.mean(np.abs(null)>=abs(d.mean())-1e-12))
    a=np.asarray(a);b=np.asarray(b);y=np.r_[a,b];n=len(a);N=len(y);idx=np.argsort(np.random.default_rng(20261113).random((19999,N)),axis=1)[:,:n];sa=y[idx].sum(1);sq=(y[idx]**2).sum(1);ma=sa/n;mb=(y.sum()-sa)/(N-n);va=(sq-n*ma**2)/(n-1);vb=((y*y).sum()-sq-(N-n)*mb**2)/(N-n-1);null=(ma-mb)/np.sqrt(np.maximum(va/n+vb/(N-n),1e-20));obs=(a.mean()-b.mean())/np.sqrt(a.var(ddof=1)/n+b.var(ddof=1)/len(b));return float((1+np.sum(abs(null)>=abs(obs)-1e-10))/20000)

data={};coverage=[];references=[]
def scores_from_signed(m,signed,variable,method='Primary signed gene z'):
    gp=m.copy();up=variable&(sign>0);down=variable&(sign<0);total=variable.sum()
    for target,gs in sets.items():
        mask=variable&np.isin(genes,gs);eligible=int(mask.sum())>=max(5,math.ceil(.5*len(gs)))
        if target=='CRRM388':
            eligible=total>=math.ceil(.7*388) and up.sum()>=140 and down.sum()>=54
            if method=='Common detection-supported genes':eligible=up.sum()>=5 and down.sum()>=5
        gp[target]=np.mean(signed[:,mask],axis=1) if eligible else np.nan
        coverage.append(dict(source=m.source_label.iloc[0],method=method,score=target,nominal=len(gs),scored=int(mask.sum()),eligible=eligible))
    gp['Up contribution']=np.mean(signed[:,up],axis=1)*up.sum()/total;gp['Down contribution']=np.mean(signed[:,down],axis=1)*down.sum()/total
    ok=gp.CRRM388.notna();assert np.allclose(gp.loc[ok,'CRRM388'],gp.loc[ok,'Up contribution']+gp.loc[ok,'Down contribution'],atol=1e-12)
    return gp

for source in ['GSE38900_GPL6884','GSE42026']:
    x=pd.read_csv(OUT/'data'/f'{source}_log_expression.csv.gz',index_col=0);m=pd.read_csv(OUT/'data'/f'{source}_sample_metadata.csv');assert x.columns.tolist()==m.geo_accession.tolist()
    full=x.to_numpy().T;ids=x.index.to_numpy();lookup={g:i for i,g in enumerate(ids)};measured=np.isin(genes,ids);xx=np.full((len(m),388),np.nan);xx[:,measured]=full[:,[lookup[g] for g in genes[measured]]]
    mu=np.nanmean(xx,axis=0);sd=np.nanstd(xx,axis=0,ddof=1);variable=measured&np.isfinite(sd)&(sd>1e-8);signed=np.full_like(xx,np.nan);signed[:,variable]=(xx[:,variable]-mu[variable])/sd[variable]*sign[variable]
    m['source_label']=source;m['analysis_group']=m.condition;m['analysis_unit']=source+'::'+m.participant;m['analysis_celltype']='Whole blood';m['time']='acute';m['cells']=np.nan;m['view']='canonical';m['unit_resolved']=True;m['profile_id']=np.arange(len(m))
    gp=scores_from_signed(m,signed,variable);rank=(stats.rankdata(full,axis=1)-1)/(full.shape[1]-1)-.5;rr=np.full_like(xx,np.nan);rr[:,measured]=rank[:,[lookup[g] for g in genes[measured]]]*sign[measured];rankgp=scores_from_signed(m,rr,measured,'Signed transcriptome rank')
    data[source]=dict(meta=m,scores=gp,signed=signed,variable=variable,rank=rankgp,measured=measured,expression=xx,rankgene=rr)
    for j,g in enumerate(genes):references.append(dict(source=source,human_ensembl=g,mean=mu[j],SD=sd[j],measured=measured[j],variable=variable[j]))

profilemeta=pd.read_csv(OUT/'data/profile_metadata.csv')
for source in ['Infant_PBMC','CD1c_DC']:
    with np.load(R4/'data'/f'native_profiles_{source}.npz',allow_pickle=False) as z:
        m=profilemeta.set_index('profile_id').loc[z['profile_ids']].reset_index();xx=z['mean_log_expression'];measured=z['common_measured'];detection=z['expressing_cells']/z['cells'][:,None]
    eligible=m.view.eq('canonical')&m.unit_resolved&m.cells.ge(30)
    if source=='Infant_PBMC':eligible&=m.analysis_group.isin(['Healthy','RSV','COVID'])
    else:eligible&=m.analysis_group.isin(['Untreated','RSV','Flu'])&m.time.isin(['4h','8h','18h'])
    ref=m.loc[eligible];w=1/ref.groupby('analysis_unit').analysis_unit.transform('size').to_numpy();w/=w.sum();x0=np.nan_to_num(xx[eligible],nan=0.);mu=np.sum(x0*w[:,None],axis=0);sd=np.sqrt(np.sum((x0-mu)**2*w[:,None],axis=0));variable=measured&np.isfinite(sd)&(sd>1e-8)
    signed=np.full_like(xx,np.nan);signed[:,variable]=(xx[:,variable]-mu[variable])/sd[variable]*sign[variable];gp=scores_from_signed(m,signed,variable)
    oldrank=pd.read_csv(R4/'results/rank_score_profiles.csv');oldrank=oldrank[oldrank.source_label.eq(source)].set_index('profile_id').loc[m.profile_id];rankgp=m.copy()
    for target in sets:rankgp[target]=oldrank[target].to_numpy()
    nu=int((measured&(sign>0)).sum());nd=int((measured&(sign<0)).sum());rankgp['Up contribution']=rankgp.Up*nu/(nu+nd);rankgp['Down contribution']=rankgp.Down*nd/(nu+nd)
    supported=np.zeros(388,bool)
    ctypes=['B_cell','T_cell','NK_cell','Monocyte'] if source=='Infant_PBMC' else ['Dendritic_cell']
    for cell in ctypes:
        use=eligible&m.analysis_celltype.eq(cell);unitdet=[detection[np.array(list(ix))].mean(0) for _,ix in m.loc[use].groupby('analysis_unit').groups.items()]
        if len(unitdet):supported|=(np.array(unitdet)>=.1).sum(0)>=3
    supported&=variable;supportgp=scores_from_signed(m,signed,supported,'Common detection-supported genes')
    data[source]=dict(meta=m,scores=gp,signed=signed,variable=variable,rank=rankgp,measured=measured,expression=xx,supported=supportgp,detection=detection)
    for j,g in enumerate(genes):references.append(dict(source=source,human_ensembl=g,mean=mu[j],SD=sd[j],measured=measured[j],variable=variable[j],detection_supported=supported[j]))

pd.DataFrame(references).to_csv(OUT/'data/gene_reference_parameters.csv',index=False);pd.DataFrame(coverage).to_csv(OUT/'results/measurement_coverage.csv',index=False)
public_profiles=pd.concat([d['scores'] for d in data.values()],ignore_index=True)
public_profiles.drop(columns=[c for c in public_profiles if c.startswith('contact_')]).to_csv(OUT/'results/profile_scores.csv',index=False)

def units(spec,gp=None,min_cells=30,omit=None):
    d=data[spec['source']];g=(d['scores'] if gp is None else gp).copy();keep=g.analysis_celltype.eq(spec['cell'])&g.view.eq('canonical')&g.unit_resolved
    if spec['source'] in ['Infant_PBMC','CD1c_DC']:keep&=g.cells.ge(min_cells)
    if spec['source']=='CD1c_DC':keep&=g.time.isin(['4h','8h','18h'])
    if omit is not None:keep&=g.analysis_unit.ne(omit)
    g=g[keep];numeric=targets if all(t in g for t in targets) else [t for t in targets if t in g]
    a=g[g.analysis_group.eq(spec['case'])];b=g[g.analysis_group.eq(spec['control'])]
    if spec['paired']:
        aa=a.groupby(['analysis_unit','time'])[numeric].mean();bb=b.groupby(['analysis_unit','time'])[numeric].mean();common=aa.index.intersection(bb.index);aa=aa.loc[common];bb=bb.loc[common]
        # Restrict primary DC to donors retaining all three matched times; thresholds remain separate sensitivities.
        good=aa.groupby(level=0).size();good=set(good[good==3].index);aa=aa[aa.index.get_level_values(0).isin(good)];bb=bb[bb.index.get_level_values(0).isin(good)]
        return aa.groupby(level=0).mean(),bb.groupby(level=0).mean(),len(common)
    return a.set_index('analysis_unit')[numeric],b.set_index('analysis_unit')[numeric],0

rows=[];ur=[]
for s in plan['specifications']:
    a,b,nt=units(s)
    for t in targets:
        r=dict(key=s['key'],source=s['source'],cell=s['cell'],comparison=s['case']+' minus '+s['control'],score=t,family=s['total_family'],paired=s['paired'],n_case=len(a),n_control=len(b),matched_strata=nt,**infer(a[t],b[t],s['paired']))
        if t=='CRRM388' and s['total_family']=='Direct virus total':r['resampling_P']=permutations(a[t],b[t],s['paired']) if min(len(a),len(b))>=3 else np.nan
        rows.append(r)
        for group,z in [('case',a),('control',b)]:
            for unit,v in z[t].items():ur.append(dict(key=s['key'],score=t,role=group,analysis_unit=unit,value=v))
effects=pd.DataFrame(rows)
for family in effects.family.unique():
    direct=family=='Direct virus total';ix=effects.family.eq(family)&effects.score.eq('CRRM388');effects.loc[ix,'adjusted_P']=bh(effects.loc[ix,'P']);effects.loc[ix,'correction_family']=family+' total BH'+str(ix.sum())
    if direct:effects.loc[ix,'resampling_adjusted_P']=bh(effects.loc[ix,'resampling_P'])
    for label,tt in [('branches',branches),('IFN components',['IFN-list overlap','Outside IFN lists']),('direction contributions',['Up contribution','Down contribution'])]:
        ix=effects.family.eq(family)&effects.score.isin(tt);effects.loc[ix,'adjusted_P']=bh(effects.loc[ix,'P']);effects.loc[ix,'correction_family']=family+' '+label+' BH'+str(ix.sum())
effects.to_csv(OUT/'results/virus_comparison_effects.csv',index=False);pd.DataFrame(ur).to_csv(OUT/'results/participant_donor_values.csv',index=False)
print(effects[effects.score.eq('CRRM388')][['key','n_case','n_control','estimate','lower','upper','adjusted_P','resampling_P']].to_string(index=False),flush=True)

# All frozen gene hypotheses: direct tests, never significant-in-only-one-group classification.
gene_rows=[]
def gene_welch(source,case,control,key):
    d=data[source];m=d['meta'];x=d['signed'];a=x[m.analysis_group.eq(case)];b=x[m.analysis_group.eq(control)];p=np.full(388,np.nan);est=np.full(388,np.nan);lo=est.copy();hi=est.copy()
    base=d['variable'].copy();est[base]=a[:,base].mean(0)-b[:,base].mean(0);mask=base.copy();mask&=(np.ptp(a,axis=0)>1e-12)|(np.ptp(b,axis=0)>1e-12);t=stats.ttest_ind(a[:,mask],b[:,mask],axis=0,equal_var=False);p[mask]=t.pvalue;est[mask]=a[:,mask].mean(0)-b[:,mask].mean(0)
    va=a[:,mask].var(0,ddof=1)/len(a);vb=b[:,mask].var(0,ddof=1)/len(b);se=np.sqrt(va+vb);df=(va+vb)**2/(va*va/(len(a)-1)+vb*vb/(len(b)-1));q=stats.t.ppf(.975,df);lo[mask]=est[mask]-q*se;hi[mask]=est[mask]+q*se
    return [dict(key=key,source=source,case=case,control=control,human_ensembl=g,display_symbol=core.iloc[j].display_symbol,direction=core.iloc[j].frozen_direction,measured=bool(base[j]),contrast_estimable=bool(mask[j]),aligned_gene_effect=est[j],lower=lo[j],upper=hi[j],P=p[j]) for j,g in enumerate(genes)]
for v in ['Influenza','Rhinovirus']:gene_rows+=gene_welch('GSE38900_GPL6884','RSV',v,'Dallas_RSV_'+v)
gd=pd.DataFrame(gene_rows);gd['adjusted_P']=bh(gd.P,776);gd.to_csv(OUT/'results/discovery_gene_differences.csv',index=False)
locked=gd[gd.key.eq('Dallas_RSV_Influenza')&gd.adjusted_P.lt(.05)][['human_ensembl','display_symbol','aligned_gene_effect','adjusted_P']]
locked.to_csv(OUT/'parameters/locked_RSV_influenza_gene_findings.csv',index=False)
if not (OUT/'parameters/gene_replication_lock.json').exists():dump(OUT/'parameters/gene_replication_lock.json',dict(locked_at=datetime.now(timezone.utc).isoformat(),discovery_family=776,replication_family=388,discovery_findings=len(locked),sha256=hashlib.sha256((OUT/'parameters/locked_RSV_influenza_gene_findings.csv').read_bytes()).hexdigest(),UK_gene_contrasts_computed=False,scope='Freeze full discovery findings, not selected examples; no new CRRM members'))
if (OUT/'parameters/constant_gene_numerical_correction.json').exists():
    dump(OUT/'parameters/gene_replication_lock.json',dict(status='DETERMINISTIC_NUMERICAL_REISSUE_AFTER_UK_RESULTS_KNOWN',locked_at=datetime.now(timezone.utc).isoformat(),original_pre_UK_lock='audit/original_pre_UK_gene_lock.json',original_discovery_count=118,discovery_family=776,replication_family=388,discovery_findings=len(locked),sha256=hashlib.sha256((OUT/'parameters/locked_RSV_influenza_gene_findings.csv').read_bytes()).hexdigest(),UK_gene_results_already_known=True,outcome_based_selection=False,correction='Uniform constant-group guard; recompute same BHfamilies; frozen388unchanged'))
gu=pd.DataFrame(gene_welch('GSE42026','RSV','Influenza','UK_RSV_Influenza'));gu['adjusted_P']=bh(gu.P,388);rep=locked.merge(gu,on=['human_ensembl','display_symbol'],suffixes=('_discovery','_replication'),validate='one_to_one');rep['same_direction']=np.sign(rep.aligned_gene_effect_discovery)==np.sign(rep.aligned_gene_effect_replication);rep['replicated']=rep.same_direction&rep.adjusted_P_replication.lt(.05);rep.to_csv(OUT/'results/locked_gene_replication.csv',index=False);gu.to_csv(OUT/'results/replication_all_gene_differences.csv',index=False)
for v in ['RSV','Influenza','Rhinovirus']:gene_rows+=gene_welch('GSE38900_GPL6884',v,'Healthy','Dallas_'+v+'_Healthy')
allg=pd.DataFrame(gene_rows);allg.to_csv(OUT/'results/gene_induction_and_difference.csv',index=False)

# Uniform sensitivity families: all7 direct-virus contrasts, retain failed hypotheses.
sens=[];direct=[s for s in plan['specifications'] if s['total_family']=='Direct virus total'];common=np.logical_and.reduce([d['variable'] for d in data.values()]);np.savetxt(OUT/'parameters/common_measured_ensembl.txt',genes[common],fmt='%s')
for s in direct:
    d=data[s['source']];m=d['meta'];methods=[('Signed transcriptome rank',d['rank'],30)]
    common_gp=scores_from_signed(m,d['signed'],common,'Common measured members');methods.append(('Common measured members',common_gp,30))
    if s['source'] in ['Infant_PBMC','CD1c_DC']:methods += [('Minimum15cells',d['scores'],15),('Minimum50cells',d['scores'],50),('Detection-supported members',d['supported'],30)]
    for method,gp,threshold in methods:
        a,b,nt=units(s,gp,threshold)
        for t in ['CRRM388','IFN-list overlap','Outside IFN lists','Neutrophil degranulation']:
            sens.append(dict(key=s['key'],method=method,score=t,n_case=len(a),n_control=len(b),**infer(a[t],b[t],s['paired'])))
    if not s['paired']:
        mm=m[m.analysis_celltype.eq(s['cell'])&m.view.eq('canonical')&m.analysis_group.isin([s['case'],s['control']])].copy()
        if s['source']=='Infant_PBMC':mm=mm[mm.cells.ge(30)]
        for t in ['CRRM388','IFN-list overlap','Outside IFN lists','Neutrophil degranulation']:
            y=d['scores'].loc[mm.index,t];ok=mm.age_months.notna()&mm.sex.isin(['M','F'])
            r=hc3(y[ok],mm.loc[ok,'analysis_group'].eq(s['case']).to_numpy(),mm.loc[ok,'age_months'].to_numpy(),mm.loc[ok,'sex'].to_numpy())
            sens.append(dict(key=s['key'],method='Age/sex HC3',score=t,n_case=int((mm.analysis_group.eq(s['case'])&ok).sum()),n_control=int((mm.analysis_group.eq(s['control'])&ok).sum()),**r))
        aa=mm[mm.analysis_group.eq(s['case'])].age_months;bb=mm[mm.analysis_group.eq(s['control'])].age_months
        if aa.notna().all() and bb.notna().all():
            low=max(aa.min(),bb.min());high=min(aa.max(),bb.max());mask=mm.age_months.between(low,high);gp=d['scores'].copy();gp=gp[~gp.analysis_celltype.eq(s['cell'])|gp.index.isin(mm.index[mask])]
            a,b,_=units(s,gp)
            sens.append(dict(key=s['key'],method='Common age support',score='CRRM388',n_case=len(a),n_control=len(b),age_min=low,age_max=high,**infer(a.CRRM388,b.CRRM388)))
    # Unit influence: remove exactly one participant/donor; recompute same estimand/reference without refitting scaling.
    a,b,_=units(s)
    for unit in sorted(set(a.index)|set(b.index)):
        aa,bb,_=units(s,omit=unit);r=infer(aa.CRRM388,bb.CRRM388,s['paired']);sens.append(dict(key=s['key'],method='Leave one biological unit out',score='CRRM388',removed_unit=unit,n_case=len(aa),n_control=len(bb),**r))
sf=pd.DataFrame(sens)
for (method,target),ix in sf[~sf.method.eq('Leave one biological unit out')].groupby(['method','score']).groups.items():sf.loc[ix,'adjusted_P_BH7']=bh(sf.loc[ix,'P'],7)
sf.to_csv(OUT/'results/sensitivity_effects.csv',index=False)

# Fixed four-lineage PBMC composition, all53 eligible participants.
d=data['Infant_PBMC'];g=d['scores'];four=g[g.view.eq('canonical')&g.analysis_celltype.isin(['B_cell','T_cell','NK_cell','Monocyte'])&g.analysis_group.isin(['Healthy','RSV','COVID'])&g.cells.ge(30)].copy()
complete=four.groupby('analysis_unit').analysis_celltype.nunique();complete=set(complete[complete.eq(4)].index);four=four[four.analysis_unit.isin(complete)].copy();assert four.groupby('analysis_unit').analysis_celltype.nunique().eq(4).all();four['fraction']=four.cells/four.groupby('analysis_unit').cells.transform('sum');ref=four[four.analysis_group.eq('Healthy')].groupby('analysis_celltype').fraction.mean();assert abs(ref.sum()-1)<1e-12
composition=[];comp_effect=[]
for target in ['CRRM388','IFN-list overlap','Outside IFN lists']:
    means=four.pivot(index='analysis_unit',columns='analysis_celltype',values=target);fractions=four.pivot(index='analysis_unit',columns='analysis_celltype',values='fraction').reindex(columns=means.columns);groups=four.drop_duplicates('analysis_unit').set_index('analysis_unit').loc[means.index,'analysis_group'];w=ref.reindex(means.columns).to_numpy()
    for method,value in [('Observed composition',(means*fractions).sum(axis=1)),('Fixed healthy composition',means@w)]:
        for unit in means.index:composition.append(dict(analysis_unit=unit,group=groups.loc[unit],score=target,method=method,value=value.loc[unit]))
        a=np.flatnonzero(groups.eq('RSV'));b=np.flatnonzero(groups.eq('COVID'));h=np.flatnonzero(groups.eq('Healthy'));boot=[];random=np.random.default_rng(20261113)
        for _ in range(2000):
            ia=random.choice(a,len(a),replace=True);ib=random.choice(b,len(b),replace=True);ih=random.choice(h,len(h),replace=True);ww=fractions.iloc[ih].mean().to_numpy()
            aa=(means.iloc[ia].to_numpy()*(fractions.iloc[ia].to_numpy() if method=='Observed composition' else ww)).sum(1).mean();bb=(means.iloc[ib].to_numpy()*(fractions.iloc[ib].to_numpy() if method=='Observed composition' else ww)).sum(1).mean();boot.append(aa-bb)
        comp_effect.append(dict(score=target,method=method,estimate=float(value[groups.eq('RSV')].mean()-value[groups.eq('COVID')].mean()),lower=float(np.quantile(boot,.025)),upper=float(np.quantile(boot,.975)),n_RSV=len(a),n_COVID=len(b),n_healthy_reference=len(h),CI='Stratified participant bootstrap2000; healthy reference re-estimated in every fixed-composition bootstrap'))
pd.DataFrame(composition).to_csv(OUT/'results/PBMC_composition_values.csv',index=False);pd.DataFrame(comp_effect).to_csv(OUT/'results/PBMC_composition_effects.csv',index=False);ref.to_csv(OUT/'results/PBMC_fixed_healthy_reference.csv')

# Detection and per-gene signed contribution to direct effects; no frozen member deletions in primary score.
contrib=[]
for s in direct:
    d=data[s['source']];m=d['meta'];gg=pd.concat([m.copy(),pd.DataFrame(d['signed'],columns=genes,index=m.index)],axis=1)
    temp=d['scores'];d['scores']=gg
    # Same unit extraction for genes via the original target keys temporarily replaced by gene columns.
    oldtargets=targets;targets=list(genes);a,b,_=units(s);targets=oldtargets;d['scores']=temp
    eff=a.to_numpy().mean(0)-b.to_numpy().mean(0);n=d['variable'].sum();dv=eff/n
    aa=m.analysis_celltype.eq(s['cell'])&m.analysis_group.eq(s['case'])&m.view.eq('canonical');bb=m.analysis_celltype.eq(s['cell'])&m.analysis_group.eq(s['control'])&m.view.eq('canonical')
    if s['source'] in ['Infant_PBMC','CD1c_DC']:
        aa&=m.cells.ge(30);bb&=m.cells.ge(30)
    if s['source']=='CD1c_DC':
        aa&=m.time.isin(['4h','8h','18h']);bb&=m.time.isin(['4h','8h','18h'])
    for j,gene in enumerate(genes):
        row=dict(key=s['key'],human_ensembl=gene,display_symbol=core.iloc[j].display_symbol,gene_effect=eff[j],contribution=dv[j],remaining_total_after_gene_omission=(np.nansum(eff)-eff[j])/(n-1) if d['variable'][j] else np.nan)
        if 'detection' in d:row.update(case_detection_fraction=float(d['detection'][aa,j].mean()),control_detection_fraction=float(d['detection'][bb,j].mean()))
        contrib.append(row)
pd.DataFrame(contrib).to_csv(OUT/'results/gene_contributions_and_detection.csv',index=False)
pd.DataFrame(coverage).to_csv(OUT/'results/measurement_coverage.csv',index=False)
dump(OUT/'audit/analysis_receipt.json',dict(status='PASS',completed_at=datetime.now(timezone.utc).isoformat(),direct_total_hypotheses=7,common_members=int(common.sum()),discovery_differential_genes=len(locked),replicated_genes=int(rep.replicated.sum()),profile_rows=sum(len(d['scores']) for d in data.values()),source_scaling_separate=True,unchanged_frozen_members=388,raw_resources_modified=False))
print('Gene discovery/replication:',len(locked),int(rep.replicated.sum()),'common_members',int(common.sum()),flush=True)
