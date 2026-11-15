from pathlib import Path
import csv, hashlib, json, re, shutil

OUT = (Path(__import__('os').environ['CRRM_WORKSPACE']) / 'phase24_streamlined_tables_and_citations_20261112')
REV = Path(str(Path(__import__('os').environ['CRRM_WORKSPACE']) / ''))
OLD = REV / 'phase23_Result3_integrated_and_table_catalog_20261112'
PHASES = {'S2':'phase18_Result2_external_response_20261112', 'S3':'phase20_Result3_clinical_analysis_20261112', 'S4':'phase22_Result3_incremental_information_20261112'}
SHA = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
def dump(p, value): p.write_text(json.dumps(value, ensure_ascii=False, indent=2)+'\n')
def write_csv(p, rows, headers=None):
    with p.open('w', newline='') as f:
        w=csv.DictWriter(f, fieldnames=headers or list(rows[0])); w.writeheader(); w.writerows(rows)
def read_table(table, sheet):
    p=REV/PHASES[table]/'tables'/(sheet+'.csv')
    rows=list(csv.DictReader(p.open()))
    shutil.copy2(p, OUT/'evidence/technical_details'/p.name)
    inputs.append({'path':str(p),'sha256':SHA(p)})
    return rows

HEADERS=['Analysis','Study','Target or contrast','Model or metric','Analysis units','Measured members','Estimate','Lower bound','Upper bound','Interval meaning','Model P','Adjusted P','Correction family','Sensitivity P','Sensitivity adjusted P','Notes']
NUMERIC=['Measured members','Estimate','Lower bound','Upper bound','Model P','Adjusted P','Sensitivity P','Sensitivity adjusted P']
inputs=[]; provenance=[]; specs=[]
def add(rows, table, old_sheet, source_index, values, original):
    row={k:'' for k in HEADERS}; row.update(values)
    assert set(row)==set(HEADERS)
    rows.append(row)
    provenance.append({'table':table,'new_sheet':table+'A_Summary' if table!='S4' else 'S4_Key_results','new_data_row':len(rows)+2,'old_sheet':old_sheet,'old_data_row':source_index+3,'original':original,'output':row})

primary=read_table('S2','S2A_Primary_effects'); contexts={r['Primary context']:r for r in primary}
s2=[]
for i,r in enumerate(primary):
    add(s2,'S2','S2A_Primary_effects',i,{'Analysis':'Primary total response','Study':r['Source'],'Target or contrast':r['Primary context'],'Model or metric':'Mean signed CRRM score change','Analysis units':r['Analysis units'],'Measured members':r['Scored CRRM members'],'Estimate':r['Mean score change'],'Lower bound':r['95% CI lower'],'Upper bound':r['95% CI upper'],'Interval meaning':'Pointwise model-based 95% t CI','Model P':r['Model P'],'Adjusted P':r['Primary Holm P (7 tests)'],'Correction family':'Holm: 7 primary total-score tests','Sensitivity P':r['Resampling sensitivity P'],'Sensitivity adjusted P':r['Resampling Holm P (7 tests)'],'Notes':'Resampling: '+r['Resampling method']+'; sensitivity Holm family: 7 tests'},r)
for i,r in enumerate(read_table('S2','S2C_Functional_responses')):
    context=contexts[r['Primary context']]
    add(s2,'S2','S2C_Functional_responses',i,{'Analysis':'Fixed functional branch','Study':context['Source'],'Target or contrast':r['Primary context']+' | '+r['Fixed functional branch'],'Model or metric':'Mean signed branch score change','Analysis units':context['Analysis units'],'Measured members':r['Scored branch members'],'Estimate':r['Mean branch score change'],'Lower bound':r['95% CI lower'],'Upper bound':r['95% CI upper'],'Interval meaning':'Pointwise model-based 95% t CI','Model P':r['Model P'],'Adjusted P':r['BH-adjusted P (49 tests)'],'Correction family':'BH: 49 primary branch-by-context tests','Notes':'Frozen branch members: '+r['Frozen branch members']},r)
for i,r in enumerate(read_table('S2','S2D_Context_and_components')):
    add(s2,'S2','S2D_Context_and_components',i,{'Analysis':r['Analysis'],'Study':r['Source'],'Target or contrast':r['Time']+' | '+r['Subgroup or contrast'],'Model or metric':'Weighted score contribution' if 'contribution' in r['Analysis'] else 'Mean signed CRRM score change','Analysis units':r['Analysis units'],'Estimate':r['Mean score change'],'Lower bound':r['95% CI lower'],'Upper bound':r['95% CI upper'],'Interval meaning':'Pointwise model-based 95% t CI','Model P':r['Model P'],'Adjusted P':r['Adjusted P'],'Correction family':r['Adjustment family'],'Notes':r['Interpretation note']},r)
assert len(s2)==103

s3=[]
clinical=read_table('S3','S3A_Clinical_effects')
units_by_study={r['Source']:r['Analysis units'] for r in clinical if r['Analysis'] in ['Primary strict severity','Hospitalization sensitivity']}
recovery_units=next(r['Analysis units'] for r in clinical if r['Analysis']=='Paired recovery')
for i,r in enumerate(clinical):
    resample_note=('Resampling: separate Holm2 family across the two primary strict-severity tests' if r['Analysis']=='Primary strict severity' else ('Separate single paired sign-flip recovery test' if r['Analysis']=='Paired recovery' else 'Separate descriptive hospitalization resampling sensitivity')) if r['Resampling P'] else ''
    add(s3,'S3','S3A_Clinical_effects',i,{'Analysis':r['Analysis'],'Study':r['Source'],'Target or contrast':r['Comparison'],'Model or metric':'Mean signed CRRM score change','Analysis units':r['Analysis units'],'Measured members':r['Scored members'],'Estimate':r['Mean score change'],'Lower bound':r['95% CI lower'],'Upper bound':r['95% CI upper'],'Interval meaning':'Pointwise HC3 95% t CI' if 'HC3' in r['Analysis'] else ('Paired t 95% CI' if 'recovery' in r['Analysis'].lower() else 'Pointwise Welch 95% CI'),'Model P':r['Model P'],'Adjusted P':r['Adjusted P'],'Correction family':r['Adjustment family'],'Sensitivity P':r['Resampling P'],'Sensitivity adjusted P':r['Resampling adjusted P'],'Notes':resample_note},r)
for i,r in enumerate(read_table('S3','S3B_Functional_responses')):
    add(s3,'S3','S3B_Functional_responses',i,{'Analysis':'Functional branch | '+r['Analysis'],'Study':r['Source'],'Target or contrast':r['Comparison']+' | '+r['Fixed branch'],'Model or metric':'Covariate-adjusted branch score difference' if 'HC3' in r['Analysis'] else 'Mean signed branch score change','Analysis units':recovery_units if 'recovery' in r['Comparison'].lower() else units_by_study[r['Source']],'Measured members':r['Scored branch members'],'Estimate':r['Mean branch score change'],'Lower bound':r['95% CI lower'],'Upper bound':r['95% CI upper'],'Interval meaning':'Paired t 95% CI' if 'recovery' in r['Comparison'].lower() else ('Pointwise HC3 95% t CI' if 'HC3' in r['Analysis'] else 'Pointwise Welch 95% CI'),'Model P':r['Model P'],'Adjusted P':r['BH-adjusted P'],'Correction family':r['Adjustment family'],'Sensitivity P':r['Exact sign-flip P'],'Sensitivity adjusted P':r['Exact sign-flip BH P'],'Notes':'Frozen branch members: '+r['Frozen branch members']+('; exact paired sign-flip sensitivity: separate BH7 family' if r['Exact sign-flip P'] else '')},r)
for i,r in enumerate(read_table('S3','S3D_Context_sensitivity')):
    add(s3,'S3','S3D_Context_sensitivity',i,{'Analysis':r['Analysis'],'Study':r['Source'],'Target or contrast':r['Subgroup or contrast'],'Model or metric':r['Metric'],'Analysis units':r['Analysis units'],'Estimate':r['Estimate'],'Lower bound':r['Lower bound'],'Upper bound':r['Upper bound'],'Interval meaning':r['Interval meaning'],'Model P':r['Model P'],'Adjusted P':r['Adjusted P'],'Correction family':'See Notes for the original test family' if r['Adjusted P'] else 'Descriptive sensitivity; no adjusted P reported','Notes':r['Interpretation note']},r)
assert len(s3)==77

s4=[]
for i,r in enumerate(read_table('S4','S4A_Models')):
    groups='Severe/Mild' if 'Hospital' not in r['Role'] else 'Inpatient/Outpatient'
    add(s4,'S4','S4A_Models',i,{'Analysis':'Conditional association | '+r['Role'],'Study':r['Study'],'Target or contrast':r['Target'],'Model or metric':r['Model'],'Analysis units':r['Cases']+' cases; '+groups+' = '+r['Severe or inpatient']+'/'+r['Mild or outpatient'],'Estimate':r['Estimate'],'Lower bound':r['95% CI lower'],'Upper bound':r['95% CI upper'],'Interval meaning':'Pointwise HC3 95% t CI','Model P':r['Model P'],'Adjusted P':r['Adjusted P'],'Correction family':r['Correction family'],'Notes':'Exploratory; published-panel scores were constructed here from the published gene list'},r)
for i,r in enumerate(read_table('S4','S4B_Stability')):
    exact=re.search(r'exact P = ([0-9.]+)',r['Role'])
    add(s4,'S4','S4B_Stability',i,{'Analysis':'Member omission and recovery','Study':r['Study'],'Target or contrast':r['Analysis'],'Model or metric':'Clinical-model coefficient' if 'clinical model' in r['Analysis'] else 'Mean subset score change','Analysis units':r['Cases or pairs']+(' verified donor pairs' if 'recovery' in r['Analysis'] else ' clinical cases'),'Measured members':r['Measured members'],'Estimate':r['Estimate'],'Lower bound':r['Lower bound'],'Upper bound':r['Upper bound'],'Interval meaning':r['Interval type'],'Model P':r['Model P'],'Adjusted P':r['Adjusted P'],'Correction family':'Holm2: LCN2 ablation' if 'Holm2' in r['Role'] else ('Separate secondary paired test' if exact else 'Descriptive context or influence summary'),'Sensitivity P':exact[1] if exact else '', 'Notes':r['Role']},r)
coverage=read_table('S4','S4C_Measurement')
for i,r in enumerate(read_table('S4','S4D_Heldout')):
    for label,prefix,base,aug in [('Log loss','Log loss','Baseline log loss','Augmented log loss'),('Brier loss','Brier','Baseline Brier','Augmented Brier')]:
        add(s4,'S4','S4D_Heldout',i,{'Analysis':'Held-out comparison | '+r['Role'],'Study':r['Study'],'Target or contrast':'Reference model versus model adding the LCN2-omitted subset','Model or metric':label+' reduction (baseline minus augmented)','Analysis units':r['Cases']+' cases','Estimate':r[prefix+' reduction'],'Interval meaning':'Descriptive loss difference; no interval','Correction family':'No performance P value or multiplicity test','Notes':'Baseline loss = '+r[base]+'; augmented loss = '+r[aug]+'; positive reduction favors augmentation'},r)
assert len(s4)==41

summary_specs=[('S2','S2A_Summary',s2,'Table S2A. External CRRM responses: primary effects, functional branches, and context'),('S3','S3A_Summary',s3,'Table S3A. Clinical CRRM results: score and branch associations, paired recovery, and sensitivities'),('S4','S4_Key_results',s4,'Table S4. Key exploratory robustness and incremental-information results')]
for table,sheet,rows,title in summary_specs:
    file=sheet+'.csv';write_csv(OUT/'tables'/file,rows,HEADERS)
    specs.append(dict(table=table,sheet=sheet,file=file,title=title,numeric_columns=NUMERIC,n_rows=len(rows)))
for table,old_sheet,new_sheet,title in [('S2','S2B_Gene_directions','S2B_Gene_directions','Table S2B. External gene effects for all 388 frozen CRRM members'),('S3','S3C_Gene_effects','S3B_Gene_effects','Table S3B. Clinical gene effects for all 388 frozen CRRM members')]:
    rows=read_table(table,old_sheet);file=new_sheet+'.csv';write_csv(OUT/'tables'/file,rows)
    nums=[k for k in rows[0] if all(not x[k] or re.fullmatch(r'-?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?',x[k]) for x in rows)]
    specs.append(dict(table=table,sheet=new_sheet,file=file,title=title,numeric_columns=nums,n_rows=len(rows),unchanged_gene_effects=True))
filenames={'S1':'SupplementaryTableS1_CRRM_key_results.xlsx','S2':'SupplementaryTableS2_CRRM_key_results.xlsx','S3':'SupplementaryTableS3_CRRM_clinical_key_results.xlsx','S4':'SupplementaryTableS4_CRRM_incremental_key_results.xlsx'}
shutil.copy2(OLD/'submission_tables'/filenames['S1'],OUT/'submission_tables'/filenames['S1'])
assert SHA(OUT/'submission_tables'/filenames['S1'])==SHA(OLD/'submission_tables'/filenames['S1'])
for spec in specs:spec['workbook']=filenames[spec['table']]
order={'S2':['S2A_Summary','S2B_Gene_directions'],'S3':['S3A_Summary','S3B_Gene_effects'],'S4':['S4_Key_results']}
specs.sort(key=lambda s:(s['table'],order[s['table']].index(s['sheet'])))
dump(OUT/'parameters/workbook_spec.json',specs);dump(OUT/'evidence/row_migration_map.json',provenance)
dump(OUT/'parameters/editorial_scope.json',dict(scope='Keep S1 unchanged; streamline S2/S3 to two sheets and S4 to one; update bilingual Result2/Result3 citations',new_scientific_analysis=False,source_csv_inputs=inputs,table_S1_sha256=SHA(OUT/'submission_tables'/filenames['S1']),scientific_values_recalculated=False,measurement_coverage_rows_archived=len(coverage),summary_rows={'S2':103,'S3':77,'S4':41}))
mapping=[('S1','A–D','A–D','Unchanged'),('S2','A','A','Primary effects'),('S2','B','B','Gene effects, unchanged'),('S2','C','A','Functional branches'),('S2','D','A','Context and contributions'),('S3','A','A','Clinical effects'),('S3','B','A','Functional responses'),('S3','C','B','Gene effects, unchanged'),('S3','D','A','Context and sensitivities'),('S4','A','Table S4','Conditional models'),('S4','B','Table S4','Stability and recovery'),('S4','C','Internal analysis archive','Measurement coverage; retained in full'),('S4','D','Table S4','Held-out losses')]
write_csv(OUT/'tables/sheet_migration_catalog.csv',[dict(Table=t,Previous_section=a,Current_section=b,Content=c) for t,a,b,c in mapping])
print(json.dumps({'S1_sheets':4,'S2_sheets':2,'S3_sheets':2,'S4_sheets':1,'new_summary_rows':[len(s2),len(s3),len(s4)]}))
