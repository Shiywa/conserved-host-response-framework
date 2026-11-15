"""Compare regenerated outputs with independently saved current scientific results."""
from __future__ import annotations
from pathlib import Path
import gzip,hashlib,io,json,platform,sys
import numpy as np,pandas as pd,scipy
from .crrm388 import reconstruct,verify_definition
from .reference import FrozenReference,score
from .statistics import paired,exact_signflip,adjust,welch,hc3
from .study import CLINICAL,input_study_directory

TOLERANCE=1e-10

def _numeric_difference(a,b,columns=None):
    if len(a)!=len(b):raise AssertionError(f"Different result row counts:{len(a)}vs{len(b)}")
    columns=columns or [c for c in a if c in b and pd.api.types.is_numeric_dtype(a[c]) and not pd.api.types.is_bool_dtype(a[c])]
    maximum=0.
    for c in columns:
        av=pd.to_numeric(a[c],errors='coerce').to_numpy(float);bv=pd.to_numeric(b[c],errors='coerce').to_numpy(float)
        if not np.array_equal(np.isfinite(av),np.isfinite(bv)):raise AssertionError('Missingness differs:'+c)
        finite=np.isfinite(av)
        if finite.any():maximum=max(maximum,float(np.max(np.abs(av[finite]-bv[finite]))))
    if maximum>TOLERANCE:raise AssertionError(f"Saved-result regression differs by{maximum}")
    return maximum

def verify(bundle,workspace=None):
    bundle=Path(bundle);genes=pd.read_csv(bundle/'frozen/crrm388_genes.tsv',sep='\t');subsets=pd.read_csv(bundle/'frozen/crrm388_subsets.tsv',sep='\t')
    definition=verify_definition(genes,subsets)
    inputs=input_study_directory()
    ev=pd.read_csv(inputs/'evidence/common_tested_discovery.tsv.gz',sep='\t');rebuilt=reconstruct(ev)
    expected=genes.set_index('human_ensembl').frozen_direction.sort_index();actual=rebuilt.set_index('human_ensembl').frozen_direction.sort_index()
    if len(ev)!=8607 or not expected.equals(actual):raise AssertionError('Frozen discovery replay mismatch')
    x=pd.read_csv(inputs/'examples/GSE231666_logCPM.csv.gz',index_col=0)
    reference=FrozenReference.fit(x,x.columns.tolist(),sd_threshold=0.)
    scores,coverage,z=score(x,reference,genes.set_index('human_ensembl'))
    saved=pd.read_csv(inputs/'examples/GSE231666_metadata_expected.tsv',sep='\t')
    delta=np.max(np.abs(scores.loc[saved.expression_column,'CRRM388'].to_numpy()-saved.score.to_numpy()))
    if delta>TOLERANCE:raise AssertionError('Primary-human saved-score mismatch')
    primary_expected=pd.read_csv(inputs/'reference_results/primary_human_statistics.tsv',sep='\t')
    rows=[]
    for age in ['Y0','Y1']:
        m=saved[saved.chronological_age.eq(age)].copy();m['score_new']=scores.loc[m.expression_column,'CRRM388'].to_numpy()
        pivot=m.pivot(index='pair_id',columns='condition',values='score_new');changes=(pivot.RSV-pivot.mock).to_numpy()
        r=paired(changes);r.update(age=age,exact_signflip_p=exact_signflip(changes));rows.append(r)
    new=pd.DataFrame(rows);new['adjusted_p']=adjust(new.exact_signflip_p,'holm',family_size=2)
    old=primary_expected.set_index('age').loc[new.age]
    human_stat_error=max(float(np.max(abs(new.effect.to_numpy()-old.mean_delta.to_numpy()))),float(np.max(abs(new.adjusted_p.to_numpy()-old.adjusted_p.to_numpy()))))
    if human_stat_error>TOLERANCE:raise AssertionError('Primary-human exact inference mismatch')
    report=dict(status='PASS',scope='Discovery evidence and processed-expression replay;not a fresh raw-count or raw-single-cell rerun',tolerance=TOLERANCE,
      definition=definition,discovery_background=8607,discovery_selected=len(rebuilt),primary_human_profiles=44,
      primary_human_score_max_error=float(delta),primary_human_exact_inference_max_error=human_stat_error,
      software=dict(python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__,scipy=scipy.__version__),clinical_replayed=False)
    if workspace is None:return report
    root=Path(workspace);manifest=pd.read_csv(inputs/'asset_manifest.tsv',sep='\t').set_index('logical_path')
    def saved_frame(logical):
        r=manifest.loc[logical];content=(inputs.parent/r.repository_path).read_bytes()
        content=gzip.decompress(content) if r.compression=='gzip' else content
        return pd.read_csv(io.BytesIO(content),float_precision='round_trip')
    table_errors={}
    for part in ['main','incremental']:
        path=root/CLINICAL/part/'tables'
        for p in sorted(path.glob('*.csv')):
            logical=f'{CLINICAL}/{part}/tables/{p.name}'
            if logical not in manifest.index or p.name.startswith('S'):continue
            generated=pd.read_csv(p,float_precision='round_trip');expected=saved_frame(logical)
            if 'dataset' in generated and 'dataset' in expected:
                # Regenerated native tables retain original ordered comparisons.
                common=[k for k in ['dataset','branch','component','target','model','sample_id','group','visit','variable','level','human_ensembl','draw','split','omitted_gene','omitted_sample'] if k in generated and k in expected]
                if common:
                    generated=generated.sort_values(common,kind='stable').reset_index(drop=True);expected=expected.sort_values(common,kind='stable').reset_index(drop=True)
            table_errors[part+'/'+p.name]=_numeric_difference(generated,expected)
    values=pd.read_csv(root/CLINICAL/'main/tables/07_sample_scores.csv');wb=values[values.dataset.eq('GSE246622')]
    if len(wb)!=262 or wb.group.value_counts().to_dict()!={'Moderate':91,'Mild':79,'Healthy':56,'Severe':36}:raise AssertionError('Wrong whole-blood population')
    if max(abs(wb.score_Total-wb['score_Up contribution']-wb['score_Down contribution']))>TOLERANCE:raise AssertionError('Weighted contributions fail')
    expr=pd.read_csv(root/CLINICAL/'main/data/processed/GSE246622_log_expression.csv.gz',index_col=0)
    ref=FrozenReference.fit(expr,expr.columns.tolist());computed,cov,z=score(expr,ref,genes.set_index('human_ensembl'))
    api_error=float(max(abs(computed.loc[wb.sample_id,'CRRM388'].to_numpy()-wb.score_Total.to_numpy())))
    if api_error>TOLERANCE or cov.scored.iloc[0]!=369 or cov.up_scored.iloc[0]!=267 or cov.down_scored.iloc[0]!=102:raise AssertionError('Whole-blood signed API/reference mismatch')
    evaluation=pd.read_csv(root/CLINICAL/'incremental/tables/19_leave_one_unit_out_predictions.csv')
    for ds in ['GSE246622','GSE246622_RMA']:
        v=evaluation[evaluation.dataset.eq(ds)]
        if len(v)!=230 or not v.training_n.eq(114).all() or not v.reference_n.eq(261).all():raise AssertionError('Wrong leave-one-profile reference')
    align=pd.read_csv(root/CLINICAL/'host/results/clinical_gene_alignment.csv');old=pd.read_csv(inputs/'reference_results/clinical_gene_alignment.csv')
    alignment_error=_numeric_difference(align,old)
    if not align.allocations.eq(70).all():raise AssertionError('Wrong clinical-bridge permutation unit')
    total=pd.read_csv(root/CLINICAL/'main/tables/09_total_clinical_effects.csv').set_index('dataset').loc['GSE246622']
    report.update(clinical_replayed=True,clinical_tables_compared=len(table_errors),clinical_max_numeric_error=max(table_errors.values()),
       whole_blood_scoring_API_max_error=api_error,whole_blood_profiles=262,whole_blood_scored=369,up=267,down=102,
       whole_blood_total={k:float(total[k]) for k in ['effect','p_value','primary_Holm_P','resampling_P','resampling_Holm_P']},
       S7E_max_numeric_error=alignment_error,S7E_whole_blood_rho=float(align.rho.iloc[0]),
       original_primary_family_size=2,bootstrap_draws=1000,Mild_splits=200,leave_one_evaluation_n=115,leave_one_training_n=114,leave_one_reference_n=261,
       GSE69606_recovery_pairs=15,table_differences=table_errors)
    return report
