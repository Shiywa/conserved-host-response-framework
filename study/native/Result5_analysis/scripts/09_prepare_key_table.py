from pathlib import Path
import pandas as pd,numpy as np,json
OUT=(Path(__import__('os').environ['CRRM_WORKSPACE']) / 'Result5_analysis');e=pd.read_csv(OUT/'results/virus_comparison_effects.csv');rows=[]
cols=['Evidence block','Study','Comparison','Cell context / gene','Score / endpoint','n case','n comparator','Estimate','95% CI lower','95% CI upper','P','Adjusted P','Correction family / interval','Measurement / note']
def add(r,block,method='Primary signed gene-z',q='adjusted_P'):
    rows.append(dict(zip(cols,[block,r.get('source',r.get('key','')),r.get('comparison',r.get('key','')),r.get('cell',''),r.get('score',''),r.get('n_case',np.nan),r.get('n_control',np.nan),r.get('estimate',np.nan),r.get('lower',np.nan),r.get('upper',np.nan),r.get('P',np.nan),r.get(q,np.nan),r.get('correction_family','BH7 sensitivity; unavailable hypotheses retained'),method])))
for _,r in e[e.score.eq('CRRM388')].iterrows():add(r,'Total response')
for _,r in e[e.key.str.startswith('Dallas_RSV_')&e.score.isin(['Up contribution','Down contribution'])].iterrows():add(r,'Signed contributions')
branches=['IFN-alpha','IFN-gamma','TNF / NF-kB','TLR signaling','Neutrophil degranulation','RTK signaling','Microtubule']
for _,r in e[e.key.isin(['Dallas_RSV_Influenza','Dallas_RSV_Rhinovirus','UK_RSV_Influenza'])&e.score.isin(branches)].iterrows():add(r,'Functional subsets')
for _,r in e[e.family.eq('Direct virus total')&e.key.str.startswith('PBMC_')&e.score.isin(['IFN-list overlap','Outside IFN lists'])].iterrows():add(r,'PBMC IFN-list components')
for _,r in e[e.family.eq('Direct virus total')&e.key.str.startswith('PBMC_')&e.score.eq('Neutrophil degranulation')].iterrows():add(r,'PBMC annotation sensitivity')
g=pd.read_csv(OUT/'results/replicated_genes_age_sex_sensitivity.csv')
for _,r in g.iterrows():
    s=pd.Series(dict(source='GSE42026',comparison='RSV minus Influenza',cell=r.display_symbol,score='Frozen-member signed expression',n_case=22,n_control=19,estimate=r.aligned_gene_effect_replication,lower=r.lower_x if 'lower_x' in r else r.lower,upper=r.upper_x if 'upper_x' in r else r.upper,P=r.P_x if 'P_x' in r else r.P,adjusted_P=r.adjusted_P_replication,correction_family='UK BH388; discovery BH776'))
    # Explicit original UK CI/P rather than the HC3 sensitivity columns after the merge.
    original=pd.read_csv(OUT/'results/locked_gene_replication.csv');z=original[original.human_ensembl.eq(r.human_ensembl)].iloc[0];s['lower']=z.lower;s['upper']=z.upper;s['P']=z.P
    add(s,'Replicated frozen members','Same direction in Dallas/UK; Dallas age/sex HC3 BH388='+format(r.adjusted_P_BH388,'.6g'))
# Retain the initially discussed boundary member explicitly; it is NOT counted as formal replication.
gd=pd.read_csv(OUT/'results/discovery_gene_differences.csv');boundary=gd[gd.key.eq('Dallas_RSV_Influenza')&gd.display_symbol.eq('TLR2')].iloc[0];uk=pd.read_csv(OUT/'results/replication_all_gene_differences.csv');uk=uk[uk.display_symbol.eq('TLR2')].iloc[0]
z=pd.Series(dict(source='GSE38900_GPL6884',comparison='RSV minus Influenza',cell='TLR2',score='Frozen-member signed expression',n_case=91,n_control=16,estimate=boundary.aligned_gene_effect,lower=boundary.lower,upper=boundary.upper,P=boundary.P,adjusted_P=boundary.adjusted_P,correction_family='Discovery BH776'))
add(z,'Boundary member; not formal replication','Fails discovery FDR 0.05; UK same direction BH388='+format(uk.adjusted_P,'.6g')+'; not counted among 4 replicated genes')
sens=pd.read_csv(OUT/'results/sensitivity_effects.csv');sens=sens[sens.score.eq('CRRM388')&sens.method.isin(['Signed transcriptome rank','Age/sex HC3'])]
for _,r in sens.iterrows():
    spec=next(s for s in json.loads((OUT/'parameters/analysis_plan_v1.json').read_text())['specifications'] if s['key']==r.key);r['source']=spec['source'];r['cell']=spec['cell'];r['comparison']=spec['case']+' minus '+spec['control'];add(r,'Total-score robustness',r.method,'adjusted_P_BH7')
rank=pd.read_csv(OUT/'results/rank_induction_and_replication_effects.csv');r=rank[rank.key.eq('UK_RSV_Influenza')&rank.score.eq('CRRM388')].iloc[0];r['source']='GSE42026';r['cell']='Whole blood';r['comparison']='RSV minus Influenza';r['correction_family']='Independent replication total BH3';add(r,'Total-score robustness','Signed transcriptome rank')
nk=rank[rank.key.eq('PBMC_NK_cell_RSV_Healthy')&rank.score.eq('CRRM388')].iloc[0];nk['source']='Infant_PBMC';nk['cell']='NK_cell';nk['comparison']='RSV minus Healthy';nk['correction_family']='PBMC shared induction BH8';add(nk,'Method-dependent induction','Signed transcriptome rank; corrected threshold not met')
comp=pd.read_csv(OUT/'results/PBMC_composition_effects.csv')
for _,r in comp[comp.score.isin(['CRRM388','Outside IFN lists'])].iterrows():
    s=pd.Series(dict(source='Infant_PBMC',comparison='RSV minus COVID',cell='Fixed four-lineage compartment',score=r.score,n_case=r.n_RSV,n_control=r.n_COVID,estimate=r.estimate,lower=r.lower,upper=r.upper,correction_family='2000 participant-bootstrap percentile CI; no P/FDR'))
    add(s,'Composition standardization',r.method+';17 healthy reference participants')
df=pd.DataFrame(rows,columns=cols);df.to_csv(OUT/'results/TableS6_key_results.csv',index=False)
(OUT/'audit/TableS6_scope.json').write_text(json.dumps(dict(rows=len(df),sheets=1,title_row=1,header_row=2,first_data_row=3,all_testing_values_reused=True,complete_tests='results/virus_comparison_effects.csv;discovery/replication gene tables;sensitivity tables',missing_or_inapplicable_P='Blank, never zero',CI='Pointwise t intervals except explicitly labelled bootstrap rows',blocks=df['Evidence block'].value_counts().to_dict()),indent=2)+'\n');print('Table S6:',len(df),'key rows;1sheet')
