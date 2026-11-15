"""Update only the unadjusted whole-blood clinical vector and original BH3 family."""
from pathlib import Path
from itertools import combinations
import json
import numpy as np,pandas as pd
from scipy import stats
BASE=(Path(__import__('os').environ['CRRM_WORKSPACE']) / 'clinical_otherwisehealthy_262_r1')
REV=Path(str(Path(__import__('os').environ['CRRM_WORKSPACE']) / ''))
OLD=REV/'Result6_analysis';OUT=BASE/'host'
core=pd.read_csv(OLD/'parameters/frozen_CRRM388.csv');genes=core.human_ensembl.to_numpy()
sets=json.loads((OLD/'parameters/score_gene_sets.json').read_text())
arr=np.load(OLD/'data/human_score_arrays.npz',allow_pickle=False)
assert np.array_equal(arr['human_ensembl'],genes)
meta=pd.read_csv(OLD/'data/human_sample_manifest.csv')
assert np.array_equal(arr['unit'],meta.unit.to_numpy())
keep=meta.group.isin(['RSV','RSV + BRD4i']).to_numpy()
xx=arr['expression'][keep];case=meta.loc[keep,'group'].eq('RSV + BRD4i').to_numpy()
assign=np.array(list(combinations(range(8),4)));assert len(assign)==70
null=np.array([xx[ii].mean(0)-xx[np.setdiff1d(np.arange(8),ii)].mean(0) for ii in assign])
effect=xx[case].mean(0)-xx[~case].mean(0)
cv=pd.read_csv(BASE/'main/tables/12_clinical_gene_effects.csv').query("contrast_key=='GSE246622_Severe_Mild'").set_index('human_ensembl')
common=[g for g in sets['Neutrophil degranulation'] if g in cv.index and arr['variable'][np.where(genes==g)[0][0]] and np.isfinite(cv.loc[g,'log_expression_effect'])]
ix=np.array([np.where(genes==g)[0][0] for g in common]);y=cv.loc[common,'log_expression_effect'].to_numpy();x=effect[ix]
rho=stats.spearmanr(x,y).statistic
rv=np.array([stats.spearmanr(n[ix],y).statistic for n in null]);p=float(np.mean(np.abs(rv)>=abs(rho)-1e-12))
key='GSE246622_Severe_Mild';aligned=x*y<0
row=dict(clinical_contrast=key,members=len(common),rho=rho,permutation_P=p,allocations=70,expected_alignment='opposite to severity',directionally_aligned=int(aligned.sum()),nonzero_pairs=int(((x!=0)&(y!=0)).sum()))
old=pd.read_csv(OLD/'results/clinical_gene_alignment.csv',float_precision='round_trip')
al=pd.concat([pd.DataFrame([row]),old[old.clinical_contrast.ne(key)]],ignore_index=True)
pv=al.permutation_P.to_numpy();order=np.argsort(pv);q=np.empty(3);q[order]=np.minimum(1,np.minimum.accumulate((pv[order]*3/np.arange(1,4))[::-1])[::-1]);al['adjusted_P_BH3']=q
al.to_csv(OUT/'results/clinical_gene_alignment.csv',index=False)
pts=pd.DataFrame(dict(clinical_contrast=key,human_ensembl=common,display_symbol=core.set_index('human_ensembl').loc[common,'display_symbol'].to_numpy(),intervention_log_effect=x,clinical_log_effect=y,aligned=aligned))
op=pd.read_csv(OLD/'results/clinical_gene_alignment_points.csv',float_precision='round_trip')
pd.concat([pts,op[op.clinical_contrast.ne(key)]],ignore_index=True).to_csv(OUT/'results/clinical_gene_alignment_points.csv',index=False)
pd.DataFrame({'allocation':np.arange(1,71),'rho':rv}).to_csv(BASE/'audit/S7E_whole_blood_permutation_null.tsv',sep='\t',index=False)
ex=pd.DataFrame(xx[:,ix],columns=common);ex.insert(0,'BRD4i',case.astype(int));ex.insert(0,'culture_unit',meta.loc[keep,'unit'].to_numpy());ex.to_csv(BASE/'audit/S7E_eight_culture_expression.tsv',sep='\t',index=False)
pts.to_csv(BASE/'audit/S7E_whole_blood_vector.tsv',sep='\t',index=False)
(BASE/'audit/S7E_scope.json').write_text(json.dumps(dict(clinical_vector='Unadjusted native log-expression Severe–Mild,not the Fig3E adjusted vector',whole_blood_Severe=36,whole_blood_Mild=79,annotation_members=25,measured_intersection=len(common),independent_culture_units=8,label_allocations=70,adjustment_family='BH3 across all three original clinical correspondence tests',other_clinical_raw_vectors_and_P_unchanged=True,host_perturbation_effects_unchanged=True),ensure_ascii=False,indent=2)+'\n')
print(al.to_string(index=False))
