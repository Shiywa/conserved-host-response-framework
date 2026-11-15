#!/usr/bin/env python3
import itertools,hashlib,json,shutil,builtins
from pathlib import Path
from datetime import datetime
import numpy as np,pandas as pd
from scipy import stats
R=Path(str(Path(__import__('os').environ['CRRM_WORKSPACE']) / 'phase06_external_validation_figure1'));F=R.parent/'phase05_mapping_corrected_reconstruction'
checks=[]
def check(name,ok,detail=''):
 checks.append(dict(check=name,passed=bool(ok),detail=str(detail)))
 if not ok:raise AssertionError((name,detail))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
cfg=json.loads((R/'parameters/analysis_plan_v1.json').read_text());receipt=json.loads((R/'parameters/declaration_receipt.json').read_text());impl=json.loads((R/'parameters/count_implementation_receipt.json').read_text());proc=json.loads((R/'logs/external_scoring_process.json').read_text())
check('Analysis configuration unchanged since declaration',sha(R/'parameters/analysis_plan_v1.json')==receipt['config_sha256'])
check('Protocol precedes actual scoring',datetime.fromisoformat(cfg['declared_at'])<datetime.fromisoformat(proc['started_at']) and proc['returncode']==0)
check('Count implementation fixed before joint scoring',datetime.fromisoformat(impl['declared_at'])<datetime.fromisoformat(proc['started_at']))
for key,value in cfg['inputs'].items():check('Frozen input unchanged '+key,sha(Path(value['path']))==value['sha256'])
manifest=pd.read_csv(R/'tables/01_external_sample_manifest.csv');selected=manifest[manifest.included_primary]
check('Selected6 lung/12 nasal profiles',dict(selected.dataset.value_counts())=={'GSE309353':12,'GSE270463':6})
check('Distinct selected GSMs and BioSamples',selected.gsm.nunique()==18 and selected.biosample.nunique()==18)
scores=pd.read_csv(R/'tables/04_external_sample_scores.csv',dtype={'donor_id':str});pairs=pd.read_csv(R/'tables/05_external_donor_pairs.csv',dtype={'donor_id':str});cov=pd.read_csv(R/'tables/06_external_member_coverage.csv');norm=pd.read_csv(R/'tables/07_external_normalization.csv',dtype={'donor_id':str});summary=pd.read_csv(R/'tables/09_external_paired_statistics.csv')
for cohort in ['GSE270463','GSE309353']:
 rows=selected[selected.dataset.eq(cohort)].copy();rows.donor_id=rows.donor_id.astype(str)
 native=pd.read_csv(R/'data/processed'/f'{cohort}_native_counts.csv.gz').set_index('feature_id')
 mm=pd.read_csv(R/'data/processed'/f'{cohort}_selected_manifest.csv',dtype={'donor_id':str})
 check(cohort+' unique complete count matrix',not native.index.duplicated().any() and list(native.columns)==list(mm.gsm))
 for row in rows.to_dict('records'):
  p=Path(row['local_file']);rec=json.loads(p.with_name(p.name+'.receipt.json').read_text());check('Downloaded source hash '+row['gsm'],sha(p)==rec['sha256'])
  if cohort=='GSE270463':
   d=pd.read_csv(p,sep='\t');dd=d.drop_duplicates();check(row['gsm']+' exactly200 copied rows removed',len(d)-len(dd)==200)
   c=dd.set_index('Feature').Count
  else:
   d=pd.read_csv(p,sep='\t',header=None,names=['id','symbol','count']);c=d.set_index('id')['count'];c=c[c.index.to_series().str.fullmatch(r'ENSG\d+').to_numpy()]
  check(row['gsm']+' independent raw count reconstruction',c.index.equals(native.index) and np.array_equal(c.to_numpy(),native[row['gsm']].to_numpy()))
 pooled={}
 for (donor,condition),g in mm.groupby(['donor_id','condition']):pooled[donor+'__'+condition]=native[g.gsm].sum(axis=1)
 pc=pd.DataFrame(pooled)
 nn=norm[norm.cohort.eq(cohort)&norm.method.eq('primary')].set_index('profile_id');pc=pc[nn.index]
 check(cohort+' donor-condition aggregation retains six units',pc.shape[1]==6 and np.allclose(pc.sum(),nn['lib.size'],rtol=1e-14))
 check(cohort+' TMM geometric mean one',abs(np.log(nn['norm.factors']).sum())<1e-12)
 eff=(nn['lib.size']*nn['norm.factors']).to_numpy();prior=eff/eff.mean()
 for set,nmember,key in [('priority388',388,'priority388'),('broad1980',1980,'broad1980')]:
  members=pd.read_csv(cfg['inputs'][key]['path']).set_index('human_ensembl')
  cc=cov[cov.cohort.eq(cohort)&cov.method.eq('primary')&cov['set'].eq(set)].set_index('human_ensembl')
  check(cohort+set+' frozen member coverage',len(cc)==nmember and cc.index.equals(members.index))
  check(cohort+set+' measurable mask',np.array_equal(cc.scored,cc.present&np.isfinite(cc.sd)&cc.sd.gt(0)&cc.n_profiles_nonzero.gt(0)))
  x=pd.read_csv(R/'tables'/f'gene_logCPM_{cohort}_{set}.csv').set_index('human_ensembl');z=pd.read_csv(R/'tables'/f'gene_zscores_{cohort}_{set}.csv').set_index('human_ensembl')
  check(cohort+set+' measurable IDs reproduced',list(x.index)==list(cc[cc.scored].index) and list(z.index)==list(x.index))
  true_x=np.log2((pc.loc[x.index].to_numpy()+prior)/(eff+2*prior)*1e6)
  check(cohort+set+' independent edgeR prior-scaled logCPM',np.allclose(true_x,x[nn.index].to_numpy(),atol=3e-12,rtol=1e-12))
  xx=x[nn.index].to_numpy();true_z=(xx-xx.mean(axis=1)[:,None])/xx.std(axis=1,ddof=1)[:,None]
  check(cohort+set+' independent z normalization',np.allclose(true_z,z[nn.index],atol=3e-12,rtol=1e-12))
  sign=members.loc[x.index].frozen_direction.map({'up':1,'down':-1}).to_numpy();true_score=(true_z*sign[:,None]).mean(axis=0)
  ss=scores[scores.cohort.eq(cohort)&scores.method.eq('primary')&scores['set'].eq(set)].set_index('profile_id').loc[nn.index]
  check(cohort+set+' independent signed score',np.allclose(true_score,ss.score,atol=3e-13))
  for row in pairs[pairs.cohort.eq(cohort)&pairs['set'].eq(set)].to_dict('records'):
   yy=scores[scores.cohort.eq(cohort)&scores.method.eq(row['method'])&scores['set'].eq(set)&scores.donor_id.eq(row['donor_id'])]
   delta=yy[yy.condition.eq('RSV')].score.mean()-yy[yy.condition.eq('mock')].score.mean()
   check(cohort+set+row['method']+' donor delta '+row['donor_id'],abs(delta-row['delta'])<2e-13)
signs=np.array(list(itertools.product([-1.,1.],repeat=3)))
for row in summary.to_dict('records'):
 d=pairs[pairs.cohort.eq(row['cohort'])&pairs.method.eq(row['method'])&pairs['set'].eq(row['set'])].delta.to_numpy();p=np.mean(np.abs(signs@d/3)>=abs(d.mean())-1e-12)
 check('Exact8-assignment inference '+row['cohort']+row['method']+row['set'],len(d)==3 and abs(p-row['exact_signflip_p'])<1e-14 and row['min_attainable_exact_two_sided_P']==.25)
 lo,hi=stats.t.interval(.95,2,loc=d.mean(),scale=stats.sem(d))
 check('Independent paired CI '+row['cohort']+row['method']+row['set'],np.allclose([lo,hi],[row['CI95_low'],row['CI95_high']],atol=1e-12))
for set,method in [('priority388','holm'),('broad1980','BH')]:
 x=summary[summary.method.eq('primary')&summary['set'].eq(set)];ps=x.exact_signflip_p.to_numpy();o=np.argsort(ps)
 q=np.maximum.accumulate((2-np.arange(2))*ps[o]).clip(0,1) if method=='holm' else np.minimum.accumulate((ps[o]*2/np.arange(1,3))[::-1])[::-1].clip(0,1)
 expected=np.empty_like(q);expected[o]=q;check('Independent '+method+' multiplicity '+set,np.allclose(expected,x.adjusted_p,atol=1e-14))
first=pd.read_csv(R/'evidence/GSE270463_initial_completed/tables/04_external_sample_scores.csv');current=scores[scores.cohort.eq('GSE270463')]
check('Lung scores unchanged after nasal data arrived',np.allclose(first.score,current.score,atol=1e-14))
go=pd.read_csv(R/'tables/16_GO_BP_ORA_full.csv');gmeta=json.loads((R/'GO_summary.json').read_text())
check('GO declared joint-tested background',gmeta['n_joint_tested_genes']==8607 and gmeta['n_GO_background']==8079 and len(go)==2*gmeta['n_eligible_terms'])
gp=stats.hypergeom.sf(go.k_overlap.to_numpy()-1,go.N_background.to_numpy(),go.K_term.to_numpy(),go.n_query.to_numpy())
check('Independent GO hypergeometric tail',np.allclose(gp,go.p_value,atol=1e-14,rtol=1e-10))
p=gp;o=np.argsort(p);q=np.empty_like(p);q[o]=np.minimum.accumulate((p[o]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1].clip(0,1)
check('Independent joint GO BH including zero hits',np.allclose(q,go.FDR_joint_BH,rtol=1e-10,atol=1e-14))
check('Zero-overlap GO tests retained',go.k_overlap.eq(0).any() and go.loc[go.k_overlap.eq(0),'p_value'].eq(1).all())
check('GO term size gate and query sizes',go.K_term.between(10,500).all() and dict(go.groupby('direction').n_query.first())=={'down':102,'up':271})
display=pd.read_csv(R/'tables/17_GO_display_terms.csv');check('Displayed GO terms from full significant table',builtins.set(display.GOID)<=builtins.set(go.loc[go.FDR_joint_BH.lt(.05),'GOID']) and len(display)==10)
gene=pd.read_csv(R/'tables/08_external_gene_direction.csv');check('Descriptive gene concordance matches frozen signs',np.array_equal(gene.mean_paired_logCPM_delta*gene.direction_multiplier>0,gene.direction_concordant))
directions=gene.groupby(['cohort','set','direction_multiplier']).agg(n_measurable=('direction_concordant','size'),n_concordant=('direction_concordant','sum')).reset_index();directions['direction']=directions.direction_multiplier.map({1:'up',-1:'down'});directions['fraction_concordant']=directions.n_concordant/directions.n_measurable;directions.to_csv(R/'tables/14_external_direction_summary.csv',index=False)
result=dict(status='PASS',checks=len(checks),passed=sum(c['passed'] for c in checks),details=checks)
(R/'logs/independent_analysis_checks.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
shutil.copy2(Path(__file__),R/'scripts/05_independent_analysis_checks.py')
print(json.dumps({k:v for k,v in result.items() if k!='details'},indent=2));print(directions.to_string(index=False))
