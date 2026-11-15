from pathlib import Path
import hashlib,itertools,json,math,os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import numpy as np,pandas as pd
from scipy import stats

OUT=(Path(__import__('os').environ['CRRM_WORKSPACE']) / 'phase18_Result2_external_response_20261112')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,x):p.write_text(json.dumps(x,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
def write(df,name):df.to_csv(OUT/'tables'/name,index=False)
cfg=json.loads((OUT/'parameters/analysis_plan_v1.json').read_text());receipt=json.loads((OUT/'parameters/declaration_receipt.json').read_text())
assert sha(OUT/'parameters/analysis_plan_v1.json')==receipt['plan_sha256']
core=pd.read_csv(OUT/'parameters/frozen_CRRM388.csv').set_index('human_ensembl');branches=pd.read_csv(OUT/'parameters/frozen_functional_branches.csv');branch_names=branches.branch.drop_duplicates().tolist()
rng=np.random.default_rng(20261112)
def truth(x):return pd.Series(x).astype(str).str.lower().eq('true').to_numpy()
def padjust(p,method,n=None):
    p=np.asarray(p,dtype=float);valid=np.flatnonzero(np.isfinite(p));out=np.full(len(p),np.nan)
    n=max(n or len(valid),len(valid))
    if not len(valid):return out
    order=valid[np.argsort(p[valid],kind='stable')];rank=np.arange(1,len(order)+1)
    values=np.maximum.accumulate((n-rank+1)*p[order]) if method=='holm' else np.minimum.accumulate((p[order]*n/rank)[::-1])[::-1]
    out[order]=np.minimum(values,1);return out
def one_ols(y,X,c):
    y=np.asarray(y,dtype=float);X=np.asarray(X,dtype=float);c=np.asarray(c,dtype=float)
    if len(y)<=X.shape[1] or np.linalg.matrix_rank(X)<X.shape[1] or not np.isfinite(y).all():return dict(effect=np.nan,se=np.nan,df=np.nan,CI95_low=np.nan,CI95_high=np.nan,p_value=np.nan)
    pinv=np.linalg.pinv(X);coef=pinv@y;res=y-X@coef;df=len(y)-X.shape[1];var=float((res@res)/df);se=float(np.sqrt(max(0,var*(c@np.linalg.inv(X.T@X)@c))));effect=float(c@coef)
    if se<=1e-15:p=float(effect==0);lo=hi=effect
    else:p=float(2*stats.t.sf(abs(effect/se),df));q=stats.t.ppf(.975,df);lo=effect-q*se;hi=effect+q*se
    return dict(effect=effect,se=se,df=df,CI95_low=float(lo),CI95_high=float(hi),p_value=p)
def one_welch(y,g):
    y=np.asarray(y,dtype=float);g=np.asarray(g,dtype=bool);a=y[g];b=y[~g]
    if len(a)<2 or len(b)<2 or not np.isfinite(y).all():return dict(effect=np.nan,se=np.nan,df=np.nan,CI95_low=np.nan,CI95_high=np.nan,p_value=np.nan)
    va=np.var(a,ddof=1)/len(a);vb=np.var(b,ddof=1)/len(b);se=np.sqrt(va+vb);df=(va+vb)**2/(va**2/(len(a)-1)+vb**2/(len(b)-1));eff=float(a.mean()-b.mean());q=stats.t.ppf(.975,df);p=2*stats.t.sf(abs(eff/se),df)
    return dict(effect=eff,se=float(se),df=float(df),CI95_low=float(eff-q*se),CI95_high=float(eff+q*se),p_value=float(p))
def fit(spec,metric='Total',keep=None):
    y=spec['Y'][metric].to_numpy();keep=np.ones(len(y),bool) if keep is None else np.asarray(keep)
    return one_welch(y[keep],spec['group'][keep]) if spec['kind']=='welch' else one_ols(y[keep],spec['X'][keep],spec['c'])
def signflip(values):
    values=np.asarray(values,dtype=float);n=len(values);obs=abs(values.mean())
    if n<=15:
        signs=np.asarray(list(itertools.product([-1,1],repeat=n)));null=signs@values/n;return float(np.mean(np.abs(null)>=obs-1e-12)),'Exact paired sign flips',len(null)
    signs=rng.choice([-1,1],size=(19999,n));null=signs@values/n;return float((1+np.sum(np.abs(null)>=obs-1e-12))/(len(null)+1)),'Monte Carlo paired sign flips',len(null)
def permutation(spec):
    y=spec['Y']['Total'].to_numpy();obs=abs(fit(spec)['effect'])
    if spec['kind']=='paired':return signflip(y)
    if spec['kind']=='welch':
        n1=int(spec['group'].sum());N=len(y);comb=math.comb(N,n1)
        if comb<=100000:
            indices=np.asarray(list(itertools.combinations(range(N),n1)));sums=y[indices].sum(axis=1);null=sums/n1-(y.sum()-sums)/(N-n1)
            return float(np.mean(np.abs(null)>=obs-1e-12)),'Exact group-label permutation',len(null)
        indices=np.argsort(rng.random((19999,N)),axis=1)[:,:n1];sums=y[indices].sum(axis=1);null=sums/n1-(y.sum()-sums)/(N-n1)
        return float((1+np.sum(np.abs(null)>=obs-1e-12))/(len(null)+1)),'Monte Carlo group-label permutation',len(null)
    if spec['kind']=='strain':
        null=np.zeros(19999)
        for strain in sorted(set(spec['strata'])):
            mask=spec['strata']==strain;v=y[mask];nt=int(spec['group'][mask].sum());idx=np.argsort(rng.random((19999,len(v))),axis=1)[:,:nt];sums=v[idx].sum(axis=1);null+=(sums/nt-(v.sum()-sums)/(len(v)-nt))/5
        return float((1+np.sum(np.abs(null)>=obs-1e-12))/(len(null)+1)),'Monte Carlo strain-blocked permutation',len(null)
    # Donor-level Rademacher residual flips under the constrained primary null.
    X=spec['X'];c=spec['c'];_,_,Vt=np.linalg.svd(c[None,:],full_matrices=True);N=Vt[1:].T;X0=X@N;y0=X0@np.linalg.pinv(X0)@y;res=y-y0;signs=rng.choice([-1,1],size=(19999,len(y)));w=c@np.linalg.pinv(X);null=(signs*res)@w+float(y0@w)
    return float((1+np.sum(np.abs(null)>=obs-1e-12))/(len(null)+1)),'Wild donor residual sign flips (sensitivity)',len(null)

datasets={};score_rows=[];functional_cov=[]
for dataset in json.loads((OUT/'data/processed/datasets.json').read_text()):
    x=pd.read_csv(OUT/'data/processed'/f'{dataset}_log_expression.csv.gz',index_col=0,float_precision='round_trip');m=pd.read_csv(OUT/'data/processed'/f'{dataset}_metadata.csv');assert list(x.columns)==m.sample_id.tolist()
    # Numerical amendment: do not amplify machine-precision variation in otherwise constant genes.
    ref=truth(m.reference_for_scaling);sd=x.loc[:,ref].std(axis=1,ddof=1);usable=np.isfinite(sd)&sd.gt(1e-8);x=x.loc[usable];mu=x.loc[:,ref].mean(axis=1);sd=sd.loc[usable];z=x.sub(mu,axis=0).div(sd,axis=0)
    signs=core.loc[x.index,'frozen_direction'].map({'up':1,'down':-1});directed=z.mul(signs,axis=0);up=signs.eq(1);down=~up
    scores=pd.DataFrame({'Total':directed.mean(axis=0),'Up contribution':directed.loc[up].sum(axis=0)/len(x),'Down contribution':directed.loc[down].sum(axis=0)/len(x)})
    assert np.max(np.abs(scores.Total-scores['Up contribution']-scores['Down contribution']))<1e-12
    for branch in branch_names:
        frozen=branches[branches.branch==branch].human_ensembl.tolist();measured=[g for g in frozen if g in directed.index];eligible=len(measured)>=5 and len(measured)/len(frozen)>=.5
        scores[branch]=directed.loc[measured].mean(axis=0) if eligible else np.nan
        functional_cov.append(dict(dataset=dataset,branch=branch,n_frozen=len(frozen),n_measured=len(measured),fraction_measured=len(measured)/len(frozen),eligible=eligible))
    assert len(scores)==len(m);scores=scores.reset_index(drop=True);merged=pd.concat([m,scores.add_prefix('score_')],axis=1);score_rows.append(merged)
    datasets[dataset]=dict(x=x,m=m,z=z,scores=scores,signs=signs)
    z.to_csv(OUT/'data/processed'/f'{dataset}_gene_z.csv.gz',index_label='human_ensembl',compression='gzip')
    scores.to_csv(OUT/'data/processed'/f'{dataset}_scores.csv',index=False)
write(pd.concat(score_rows,ignore_index=True),'05_sample_scores.csv');write(pd.DataFrame(functional_cov),'06_functional_branch_coverage.csv')

def paired_spec(dataset,ids0,ids1,key,label,family,**extra):
    o=datasets[dataset];assert len(ids0)==len(ids1);Y=o['scores'].iloc[ids1].reset_index(drop=True)-o['scores'].iloc[ids0].reset_index(drop=True);gene=o['x'].iloc[:,ids1].to_numpy().T-o['x'].iloc[:,ids0].to_numpy().T
    return dict(key=key,dataset=dataset,label=label,family=family,kind='paired',Y=Y,X=np.ones((len(Y),1)),c=np.ones(1),gene=gene,genes=o['x'].index.to_numpy(),units=o['m'].iloc[ids1].unit_id.astype(str).to_numpy(),n_RSV=len(ids1),n_control=len(ids0),n_units=len(ids1),**extra)
def simple_spec(dataset,mask,key,label,family):
    o=datasets[dataset];ids=np.flatnonzero(mask);g=o['m'].iloc[ids].condition.eq('RSV').to_numpy();Y=o['scores'].iloc[ids].reset_index(drop=True)
    return dict(key=key,dataset=dataset,label=label,family=family,kind='welch',Y=Y,group=g,gene=o['x'].iloc[:,ids].to_numpy().T,genes=o['x'].index.to_numpy(),units=o['m'].iloc[ids].unit_id.astype(str).to_numpy(),n_RSV=int(g.sum()),n_control=int((~g).sum()),n_units=len(ids))
def airway_pairs(hour,background='Non-CF'):
    m=datasets['GSE247298']['m'];donors=sorted(m.loc[m.background==background,'donor_id'].unique());i0=[];i1=[]
    for d in donors:i0.append(int(np.flatnonzero((m.donor_id==d)&(m.time_h==0))[0]));i1.append(int(np.flatnonzero((m.donor_id==d)&(m.time_h==hour))[0]))
    return i0,i1
def challenge_spec(day,cname,key,label):
    o=datasets['Human_challenge'];m=o['m'];primary_donors=set(m.loc[truth(m.primary),'donor_id']);donors=sorted(d for d,v in m.groupby('donor_id') if {0,day}.issubset(set(v.time_d)) and (day==7 or d in primary_donors))
    assert len(donors)==(38 if day==3 else 13)
    i0=[int(np.flatnonzero((m.donor_id==d)&(m.time_d==0))[0]) for d in donors];i1=[int(np.flatnonzero((m.donor_id==d)&(m.time_d==day))[0]) for d in donors]
    s=paired_spec('Human_challenge',i0,i1,key,label,'GSE155237_GSE166161');inf=m.iloc[i1].condition.eq('Infected').astype(int).to_numpy();batch=m.batch.eq('new').astype(int).to_numpy();bd=batch[i1]-batch[i0]
    X=np.column_stack([np.ones(len(donors)),inf,bd]) if day==3 else np.column_stack([np.ones(len(donors)),inf]);c={'infected':np.array([1.,1.,0.]) if day==3 else np.array([1.,1.]),'uninfected':np.array([1.,0.,0.]) if day==3 else np.array([1.,0.]),'interaction':np.array([0.,1.,0.]) if day==3 else np.array([0.,1.])}[cname]
    s.update(kind='challenge',X=X,c=c,infected=inf,batch_transition=bd,n_RSV=int(inf.sum()),n_control=int((1-inf).sum()),i0=i0,i1=i1,day=day,target=cname);assert np.linalg.matrix_rank(X)==X.shape[1]
    return s

specs=[];i0,i1=airway_pairs(24);specs.append(paired_spec('GSE247298',i0,i1,'Airway24h','Human airway, 24 h','GSE247298'))
o=datasets['GSE112956'];mask=o['m'].time_d.eq(1).to_numpy();ids=np.flatnonzero(mask);m=o['m'].iloc[ids];strains=sorted(m.strain.unique());cell=m.strain+'_'+m.condition;levels=[s+'_'+c for s in strains for c in ['Control','RSV']];X=np.column_stack([cell.eq(v).astype(float) for v in levels]);c=np.tile([-.2,.2],5)
specs.append(dict(key='Mouse_multistrain_D1',dataset='GSE112956',label='Mouse strains, day 1',family='GSE112956',kind='strain',Y=o['scores'].iloc[ids].reset_index(drop=True),group=m.condition.eq('RSV').to_numpy(),strata=m.strain.to_numpy(),X=X,c=c,gene=o['x'].iloc[:,ids].to_numpy().T,genes=o['x'].index.to_numpy(),units=m.unit_id.astype(str).to_numpy(),n_RSV=15,n_control=14,n_units=29))
for dataset,key,label in [('GSE119849','Mouse119849','Mouse lung, GSE119849'),('GSE140226','Mouse140226_D5','Mouse lung, day 5'),('GSE18170','Mouse18170','Mouse lung, GSE18170'),('PRJNA1271255_Lung','Hamster_lung','Hamster lung, day 5')]:
    mask=datasets[dataset]['m'].time_d.eq(5) if dataset=='GSE140226' else np.ones(len(datasets[dataset]['m']),bool);specs.append(simple_spec(dataset,mask,key,label,dataset.split('_Lung')[0]))
specs.append(challenge_spec(3,'infected','Human_D3','Human challenge, day 3'))
assert [s['key'] for s in specs]==[s['key'] for s in cfg['primary_contrasts']]
eligibility=pd.read_csv(OUT/'tables/04_coverage_summary.csv').set_index('dataset');primary=[];functions=[];components=[];genes=[];influence=[];perm_rows=[]
for s in specs:
    common=dict(contrast_key=s['key'],dataset=s['dataset'],label=s['label'],study_family=s['family'],design=s['kind'],n_RSV=s['n_RSV'],n_control=s['n_control'],n_units=s['n_units'],n_measured=len(datasets[s['dataset']]['x']))
    eligible=bool(eligibility.loc[s['dataset'],'eligible_primary']);r=fit(s);primary.append(dict(**common,**r,eligible_primary=eligible))
    pp,method,nperm=permutation(s);perm_rows.append(dict(contrast_key=s['key'],permutation_p=pp,permutation_method=method,n_permutations=nperm))
    for b in branch_names:functions.append(dict(**common,branch=b,**fit(s,b)))
    for metric in ['Up contribution','Down contribution']:components.append(dict(**common,component=metric,**fit(s,metric)))
    if s['kind']=='welch':effect=s['gene'][s['group']].mean(axis=0)-s['gene'][~s['group']].mean(axis=0)
    else:effect=s['c']@np.linalg.pinv(s['X'])@s['gene']
    for gene,eff in zip(s['genes'],effect):
        sign=1 if core.loc[gene,'frozen_direction']=='up' else -1;genes.append(dict(contrast_key=s['key'],dataset=s['dataset'],human_ensembl=gene,display_symbol=core.loc[gene,'display_symbol'],frozen_direction=core.loc[gene,'frozen_direction'],log_expression_effect=float(eff),aligned_effect=float(eff*sign),direction_concordant=bool(eff*sign>0)))
    for i,unit in enumerate(s['units']):
        keep=np.ones(len(s['Y']),bool);keep[i]=False;rloo=fit(s,keep=keep);influence.append(dict(contrast_key=s['key'],removed_unit=unit,full_effect=r['effect'],**rloo))
primary=pd.DataFrame(primary);primary['primary_Holm_P']=padjust(primary.p_value.where(primary.eligible_primary),'holm',7);primary=primary.merge(pd.DataFrame(perm_rows),on='contrast_key',validate='one_to_one');primary['permutation_Holm_P']=padjust(primary.permutation_p.where(primary.eligible_primary),'holm',7)
functions=pd.DataFrame(functions);functions['FDR_BH49']=padjust(functions.p_value,'BH',49);components=pd.DataFrame(components);components['FDR_BH14']=padjust(components.p_value,'BH',14)
write(primary,'07_primary_response_effects.csv');write(functions,'08_primary_functional_effects.csv');write(components,'09_directional_component_effects.csv');write(pd.DataFrame(genes),'10_external_gene_effects.csv');write(pd.DataFrame(influence),'11_unit_deletion_influence.csv')
gene_df=pd.DataFrame(genes);concordance=gene_df.groupby(['contrast_key','frozen_direction']).agg(n_measured=('direction_concordant','size'),n_concordant=('direction_concordant','sum')).reset_index();concordance['fraction_concordant']=concordance.n_concordant/concordance.n_measured;write(concordance,'12_direction_concordance_summary.csv')

# Common measured-member sensitivity retains the original z scores and frozen signs.
common_genes=set.intersection(*(set(datasets[s['dataset']]['x'].index) for s in specs if bool(eligibility.loc[s['dataset'],'eligible_primary'])))
common_scores={}
for name,o in datasets.items():
    keep=[g for g in o['z'].index if g in common_genes];sg=core.loc[keep,'frozen_direction'].map({'up':1,'down':-1});common_scores[name]=o['z'].loc[keep].mul(sg,axis=0).mean(axis=0).to_numpy()
common_rows=[]
for s in specs:
    o=datasets[s['dataset']];oldY=s['Y'];cs=common_scores[s['dataset']]
    if s['kind'] in ['paired','challenge']:
        if 'i0' in s:i0=s['i0'];i1=s['i1']
        else:i0,i1=airway_pairs(24)
        y=cs[i1]-cs[i0]
    elif s['key']=='Mouse_multistrain_D1':y=cs[o['m'].time_d.eq(1)]
    elif s['key']=='Mouse140226_D5':y=cs[o['m'].time_d.eq(5)]
    else:y=cs
    ss={**s,'Y':pd.DataFrame({'Total':y})};r=fit(ss);common_rows.append(dict(contrast_key=s['key'],dataset=s['dataset'],n_common=len(common_genes),n_common_up=sum(core.loc[g,'frozen_direction']=='up' for g in common_genes),n_common_down=sum(core.loc[g,'frozen_direction']=='down' for g in common_genes),full_score_correlation=float(np.corrcoef(o['scores'].Total,cs)[0,1]),**r))
write(pd.DataFrame(common_rows),'13_common_member_sensitivity.csv');write(core.loc[sorted(common_genes)].reset_index()[['human_ensembl','display_symbol','frozen_direction']],'14_common_measured_members.csv')

time_rows=[];time_pairs=[];cf_rows=[]
for background in ['Non-CF','CF']:
    for hour in [12,24,48,72,120]:
        i0,i1=airway_pairs(hour,background);s=paired_spec('GSE247298',i0,i1,f'{background}_{hour}',f'{background}, {hour} h','GSE247298');r=fit(s);p,method,np0=signflip(s['Y'].Total)
        time_rows.append(dict(background=background,time_h=hour,n_donors=4,**r,exact_signflip_p=p))
        for donor,v in zip(s['units'],s['Y'].Total):time_pairs.append(dict(background=background,time_h=hour,donor_id=donor,delta=float(v)))
for hour in [12,24,48,72,120]:
    z=pd.DataFrame(time_pairs);z=z[z.time_h==hour];r=one_welch(z.delta,z.background.eq('CF'));cf_rows.append(dict(time_h=hour,contrast='CF change minus non-CF change',**r))
time_df=pd.DataFrame(time_rows);time_df['FDR_BH_within_background']=np.nan
for bg in ['Non-CF','CF']:
    ix=time_df.background.eq(bg);time_df.loc[ix,'FDR_BH_within_background']=padjust(time_df.loc[ix,'p_value'],'BH',5)
cf_df=pd.DataFrame(cf_rows);cf_df['FDR_BH5']=padjust(cf_df.p_value,'BH',5)
write(time_df,'15_airway_time_effects.csv');write(pd.DataFrame(time_pairs),'16_airway_donor_changes.csv');write(cf_df,'17_CF_background_effects.csv')

strain_rows=[]
o=datasets['GSE112956']
for strain in sorted(o['m'].strain.unique()):
    for day in [1,5]:
        mask=o['m'].strain.eq(strain)&(o['m'].condition.eq('Control')|o['m'].time_d.eq(day));s=simple_spec('GSE112956',mask,f'{strain}_{day}',f'{strain}, day {day}','GSE112956');strain_rows.append(dict(strain=strain,time_d=day,n_control=s['n_control'],n_RSV=s['n_RSV'],**fit(s)))
strain_df=pd.DataFrame(strain_rows);strain_df['FDR_BH10']=padjust(strain_df.p_value,'BH',10);write(strain_df,'18_mouse_strain_time_effects.csv')

human_rows=[];human_pairs=[]
for day in [3,7]:
    for target in ['infected','uninfected','interaction']:
        s=challenge_spec(day,target,f'Human_D{day}_{target}',target);human_rows.append(dict(time_d=day,target=target,n_infected=s['n_RSV'],n_uninfected=s['n_control'],**fit(s)))
        if target=='infected':
            for donor,delta,inf,bd in zip(s['units'],s['Y'].Total,s['infected'],s['batch_transition']):human_pairs.append(dict(time_d=day,donor_id=donor,eventual_outcome='Infected' if inf else 'Uninfected',delta=float(delta),batch_transition=int(bd)))
human_df=pd.DataFrame(human_rows);human_df['secondary_BH']=np.nan
for day in [3,7]:
    ix=human_df.time_d.eq(day);human_df.loc[ix,'secondary_BH']=padjust(human_df.loc[ix,'p_value'],'BH',3)
write(human_df,'19_human_challenge_effects.csv');write(pd.DataFrame(human_pairs),'20_human_challenge_donor_changes.csv')

secondary=[]
for dataset,day,key,label in [('GSE140226',7,'Mouse140226_D7','Mouse lung, day 7'),('PRJNA1271255_Blood',5,'Hamster_blood','Hamster blood, day 5')]:
    mask=datasets[dataset]['m'].time_d.eq(day);s=simple_spec(dataset,mask,key,label,dataset);secondary.append(dict(contrast_key=key,dataset=dataset,**fit(s)))
write(pd.DataFrame(secondary),'21_additional_context_effects.csv')

summary=dict(status='ANALYZED',frozen_members=388,primary_tests=7,primary_positive=int((primary.effect>0).sum()),primary_Holm_supported=int(((primary.effect>0)&(primary.primary_Holm_P<.05)).sum()),functional_tests=49,functional_BH_supported=int((functions.FDR_BH49<.05).sum()),n_common_genes=len(common_genes),n_datasets=len(datasets),no_CRRM_reselection=True,language='American English',primary_results=primary.replace({np.nan:None}).to_dict('records'))
save(OUT/'analysis_summary.json',summary)
print(primary[['contrast_key','n_measured','effect','CI95_low','CI95_high','p_value','primary_Holm_P','permutation_p']].to_string(index=False));print('Common measured genes:',len(common_genes));print(human_df.to_string(index=False))
