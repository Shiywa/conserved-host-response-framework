from pathlib import Path
import json
import numpy as np
import pandas as pd
OUT=(Path(__import__('os').environ['CRRM_WORKSPACE']) / 'phase18_Result2_external_response_20261112')
read=lambda f:pd.read_csv(OUT/'tables'/f)
spec=[]
def add(sheet,title,d,notes):
 file=sheet+'.csv';d.to_csv(OUT/'tables'/file,index=False);spec.append(dict(sheet=sheet,title=title,file=file,notes=notes))
p=read('07_primary_response_effects.csv');keys=p.contrast_key.tolist();context=dict(zip(keys,p.label));ds=dict(zip(keys,p.dataset))
units=['4 donor pairs','15 RSV / 14 vehicle; equal 5-strain contrast','10 RSV / 10 control','3 RSV / 3 mock','3 RSV / 3 mock','4 live-RSV (MOCK-immunized) / 4 NC','38 donor pairs: 22 infected / 16 uninfected']
d=p[['dataset','label','n_measured','effect','CI95_low','CI95_high','p_value','primary_Holm_P','permutation_p','permutation_Holm_P','permutation_method']].copy();d.insert(2,'Analysis units',units)
d=d.rename(columns={'dataset':'Source','label':'Primary context','n_measured':'Scored CRRM members','effect':'Mean score change','CI95_low':'95% CI lower','CI95_high':'95% CI upper','p_value':'Model P','primary_Holm_P':'Primary Holm P (7 tests)','permutation_p':'Resampling sensitivity P','permutation_Holm_P':'Resampling Holm P (7 tests)','permutation_method':'Resampling method'})
add('S2A_Primary_effects','Table S2A. Primary external responses of the frozen 388-gene CRRM',d,[
 'Seven planned primary contrasts; two-sided model-based tests and 95% t confidence intervals. Scores use each dataset\'s fixed scaling reference and should not be ranked across platforms.',
 'Human challenge is one study family spanning GSE155237 and GSE166161. Thirteen repeated baseline records were collapsed; the primary analysis uses 38 unique donor pairs and adjusts for deposition batch.',
 'Human resampling P uses added wild donor-residual sign flips for the adjusted infected-donor change; it is a sensitivity analysis, not an exact permutation test. Other comparisons use paired sign flips, group labels, or strain-blocked labels.',
 'GSE247298 uses each donor\'s preinfection culture baseline; GSE112956 has 14 deposited vehicle profiles. The hamster MOCK arm received mock immunization followed by live RSV. Hamster mapping is a qualified author-symbol projection.'
])
core=pd.read_csv(OUT/'parameters/frozen_CRRM388.csv');g=read('10_external_gene_effects.csv');wide=g.pivot(index='human_ensembl',columns='contrast_key',values='log_expression_effect').reindex(index=core.human_ensembl,columns=keys)
d=core[['display_symbol','human_ensembl','frozen_direction']].copy().set_index('human_ensembl');d=d.join(wide);d=d.reset_index()[['display_symbol','human_ensembl','frozen_direction']+keys];d.columns=['Gene symbol','Human Ensembl ID','Frozen direction']+[context[k]+' log-expression effect' for k in keys]
add('S2B_Gene_directions','Table S2B. External effects for all frozen CRRM members',d,[
 'All 388 frozen members are retained. NA denotes an unavailable or unscored member, not zero expression. Effects are RSV-minus-control, paired post-minus-baseline, or the corresponding adjusted contrast on each source\'s native normalized log-expression scale.',
 'A positive external effect agrees with an up member; a negative external effect agrees with a down member. Directional concordance includes every measured member without an external differential-expression filter. Genes are correlated and were not treated as independent replicate units.',
 'GSE112956 is the equal-strain day 1 contrast, GSE140226 is day 5, and the human challenge effect is the batch-adjusted infected-donor day 3 change.'
])
b=read('08_primary_functional_effects.csv');cov=read('06_functional_branch_coverage.csv');b=b.drop(columns='n_measured').merge(cov[['dataset','branch','n_frozen','n_measured']],on=['dataset','branch'],validate='many_to_one')
d=b[['label','branch','n_frozen','n_measured','effect','CI95_low','CI95_high','p_value','FDR_BH49']].rename(columns={'label':'Primary context','branch':'Fixed functional branch','n_frozen':'Frozen branch members','n_measured':'Scored branch members','effect':'Mean branch score change','CI95_low':'95% CI lower','CI95_high':'95% CI upper','p_value':'Model P','FDR_BH49':'BH-adjusted P (49 tests)'})
add('S2C_Functional_responses','Table S2C. Responses of seven fixed CRRM functional branches',d,[
 'Branches are fixed, potentially overlapping intersections between CRRM and annotations established in Figure 1; they are not newly inferred disjoint modules.',
 'IFN-alpha and IFN-gamma, TNF/NF-kB, TLR signaling, and neutrophil degranulation use up members; RTK signaling and microtubule use sign-reversed down members. Eligibility requires at least five scored genes and at least 50% coverage.',
 'BH adjustment includes all 49 planned branch-by-primary-context tests. Functional estimates and P values use the same biological units and contrasts as Table S2A.'
])
rows=[]
def row(source,analysis,time,group,units,r,adj=np.nan,family='',note=''):
 rows.append({'Source':source,'Analysis':analysis,'Time':time,'Subgroup or contrast':group,'Analysis units':units,'Mean score change':r['effect'],'95% CI lower':r.get('CI95_low',np.nan),'95% CI upper':r.get('CI95_high',np.nan),'Model P':r.get('p_value',np.nan),'Adjusted P':adj,'Adjustment family':family,'Interpretation note':note})
for _,r in read('15_airway_time_effects.csv').iterrows():row('GSE247298','Airway time course',str(int(r.time_h))+' h',r.background,'4 donor pairs',r,r.FDR_BH_within_background,'BH: 5 time points per background','Baseline comparison; no parallel mock time course')
for _,r in read('17_CF_background_effects.csv').iterrows():row('GSE247298','Background difference',str(int(r.time_h))+' h','CF change minus non-CF change','4 CF and 4 non-CF donors',r,r.FDR_BH5,'BH: 5 background comparisons','Observational background comparison; common non-CF scaling reference')
for _,r in read('18_mouse_strain_time_effects.csv').iterrows():row('GSE112956','Strain-specific response','Day '+str(int(r.time_d)),r.strain,f'{int(r.n_RSV)} RSV / {int(r.n_control)} vehicle',r,r.FDR_BH10,'BH: 10 strain-by-time comparisons','Vehicle sampled on day 1; day 5 is not matched to a day 5 vehicle group')
for _,r in read('19_human_challenge_effects.csv').iterrows():row('GSE155237 / GSE166161','Human challenge','Day '+str(int(r.time_d)),{'infected':'Within infected-outcome donors','uninfected':'Within uninfected-outcome donors','interaction':'Infected change minus uninfected change'}[r.target],f'{int(r.n_infected)} infected / {int(r.n_uninfected)} uninfected donor pairs in model',r,r.secondary_BH,'BH: 3 contrasts per follow-up day','Day 3 within-infected primary Holm P is in S2A; day 7 is the same study family')
for _,r in read('21_additional_context_effects.csv').iterrows():row(r.dataset,'Additional context','Day 7' if r.contrast_key=='Mouse140226_D7' else 'Day 5','Untreated RSV minus mock' if r.contrast_key=='Mouse140226_D7' else 'Blood: live RSV minus NC','3 / 3 animals' if r.contrast_key=='Mouse140226_D7' else '4 / 4 animals',r,note='Secondary, unadjusted context; not an additional independent primary study')
for _,r in read('09_directional_component_effects.csv').iterrows():row(r.dataset,'Additive score contribution',context[r.contrast_key],r.component,units[keys.index(r.contrast_key)],r,r.FDR_BH14,'BH: 14 contribution comparisons','Weighted by the measured direction-specific member count; up plus down contributions equal the total')
add('S2D_Context_and_components','Table S2D. Key contextual responses and additive score contributions',pd.DataFrame(rows),[
 'Secondary analyses retain their declared comparison families; they do not replace an unsupported primary contrast. All confidence intervals are model-based, pointwise 95% intervals.',
 'Human baseline is the prechallenge nasal-curettage specimen (7–14 days before inoculation), even where deposited metadata labels it day 0. Eventual infection outcome is a donor attribute, not a baseline infection status.',
 'Hamster blood retains 350 numerically variable members after excluding 20 machine-precision-variance genes (SD <= 1e-8); this numerical correction affects no primary contrast or frozen member definition.'
])
(OUT/'parameters/workbook_spec.json').write_text(json.dumps(spec,ensure_ascii=False,indent=2)+'\n');print([(s['sheet'],len(pd.read_csv(OUT/'tables'/s['file']))) for s in spec])
