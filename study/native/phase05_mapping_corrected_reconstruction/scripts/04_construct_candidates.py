#!/usr/bin/env python3
"""Build candidates from full corrected feature tables; never restrict to legacy 70."""
import argparse,csv,hashlib,json,math
from pathlib import Path

def read(path):
    with path.open(encoding='utf-8-sig',newline='') as h:return list(csv.DictReader(h))
def write(path,rows,columns):
    with path.open('w',encoding='utf-8-sig',newline='') as h:
        w=csv.DictWriter(h,fieldnames=columns);w.writeheader();w.writerows(rows)
def flag(x):return str(x).lower() in {'true','1'}
def num(x):
    try:return float(x)
    except (TypeError,ValueError):return math.nan
def significant(row,config):
    return num(row.get('FDR'))<config['primary_q_threshold'] and abs(num(row.get('logFC')))>config['absolute_log2FC_threshold']
def index(path):
    records=read(path);out={r['feature_id']:r for r in records}
    assert len(out)==len(records),path
    return out

def construct(root,project):
    out=root/'construction';out.mkdir(exist_ok=True)
    cfg=json.loads((root/'parameters/pre_analysis_v1.json').read_text())
    receipt=json.loads((root/'parameters/pre_analysis_declaration_receipt.json').read_text())
    assert hashlib.sha256((root/'parameters/pre_analysis_v1.json').read_bytes()).hexdigest()==receipt['config_sha256']
    hmap=index(root/'mapping/human_feature_map.csv')
    systems={
        'mouse':('M_lung_primary__infected_vs_control__edgeR_QLF.csv','mouse_feature_map.csv'),
        'cotton_rat':('R_D1_primary__age_average__edgeR_QLF.csv','cotton_rat_feature_map.csv'),
        'organoid_conditional':('O_conditional__strain_average__limma_logTPM.csv','organoid_feature_map.csv'),
    }
    by_system={}
    for name,(de,map_file) in systems.items():
        table=index(root/'de'/de)
        mapped={}
        for m in read(root/'mapping'/map_file):
            if name!='organoid_conditional' and not flag(m.get('eligible_mapping')):continue
            human=m.get('human_ensembl','');feature=m['feature_id']
            if human and feature in table:
                assert human not in mapped,(name,'Nonunique human target',human)
                mapped[human]=table[feature]
        by_system[name]=mapped
    summaries=[];freezes=[]
    branches={
      'D3_primary':'H_D3_primary__infected_change__edgeR_QLF.csv',
      'D7_exploratory':'H_D7_same_family__infected_change__edgeR_QLF.csv',
    }
    base_columns=['human_ensembl','display_symbol','human_log2FC','human_FDR','human_threshold_pass','frozen_direction','n_same_direction_supported_families','n_opposite_direction_supported_families','tier','cotton_rat_projection_caveat']
    extra_columns=[name+'_'+s for name in systems for s in ['mapped_and_tested','log2FC','FDR','same_direction_threshold_support','eligible_family_vote']]
    columns=base_columns+extra_columns
    full_by_branch={}
    for branch,filename in branches.items():
        h=index(root/'de'/filename);records=[]
        for ens,row in h.items():
            effect=num(row['logFC']);hp=significant(row,cfg);support=0;opposite=0
            record=dict(human_ensembl=ens,display_symbol=hmap.get(ens,{}).get('display_symbol',''),human_log2FC=row['logFC'],human_FDR=row['FDR'],human_threshold_pass=hp,
                frozen_direction='up' if effect>0 else 'down' if effect<0 else 'zero',
                cotton_rat_projection_caveat='Author-annotated mouse-protein bridge; not direct genome-wide one-to-one cotton-rat orthology')
            for name,lookup in by_system.items():
                v=lookup.get(ens);tested=v is not None
                same=tested and significant(v,cfg) and effect*num(v['logFC'])>0
                opp=tested and significant(v,cfg) and effect*num(v['logFC'])<0
                eligible=name!='organoid_conditional'
                if eligible:support+=int(same);opposite+=int(opp)
                record.update({name+'_mapped_and_tested':tested,name+'_log2FC':v['logFC'] if v else '',name+'_FDR':v['FDR'] if v else '',name+'_same_direction_threshold_support':same,name+'_eligible_family_vote':eligible})
            record['n_same_direction_supported_families']=support
            record['n_opposite_direction_supported_families']=opposite
            record['tier']='strict_cross_system_candidate' if hp and support>=2 else 'one_system_supported' if hp and support>=1 else 'human_response' if hp else 'not_selected'
            records.append(record)
        records.sort(key=lambda r:(not r['human_threshold_pass'],-r['n_same_direction_supported_families'],num(r['human_FDR']),r['human_ensembl']))
        full_by_branch[branch]=records
        write(out/f'{branch}_full_mapped_evidence.csv',records,columns)
        for label,minimum in [('human_response',0),('one_system_supported',1),('strict_cross_system_candidate',2)]:
            selected=[r for r in records if r['human_threshold_pass'] and r['n_same_direction_supported_families']>=minimum]
            write(out/f'{branch}_{label}.csv',selected,columns)
        core=[r for r in records if r['tier']=='strict_cross_system_candidate']
        summary=dict(branch=branch,human_features_tested=len(h),human_response_genes=sum(r['human_threshold_pass'] for r in records),one_system_supported_genes=sum(r['human_threshold_pass'] and r['n_same_direction_supported_families']>=1 for r in records),strict_candidate_genes=len(core),organoid_counts_as_vote=False,
                     status='No primary core under the declared rule' if branch=='D3_primary' and not core else 'Exploratory temporal candidate, not validated' if branch=='D7_exploratory' else 'Unvalidated primary candidate',
                     mapping_qualification='Cotton-rat support uses an explicitly qualified author-annotation bridge')
        summaries.append(summary)
        membership=out/f'{branch}_strict_cross_system_candidate.csv'
        freezes.append(dict(branch=branch,n_genes=len(core),membership_file=membership.name,sha256=hashlib.sha256(membership.read_bytes()).hexdigest(),scoring='Signed equal-weight within-dataset gene z-score mean, signs from this branch human response; no learned outcome weights',modules='Not defined; no forced five modules',status=summary['status']))
    # Rebuild the entire remaining pool for family omission. Independent family
    # DE fits are cached because they have no cross-family statistical inputs.
    loo=[]
    for branch,records in full_by_branch.items():
        original={r['human_ensembl'] for r in records if r['tier']=='strict_cross_system_candidate'}
        for omitted in ['mouse','cotton_rat']:
            remaining=[n for n in ['mouse','cotton_rat'] if n!=omitted]
            rebuilt={r['human_ensembl'] for r in records if r['human_threshold_pass'] and sum(bool(r[n+'_same_direction_threshold_support']) for n in remaining)>=2}
            union=original|rebuilt
            loo.append(dict(branch=branch,omitted_family=omitted,minimum_family_support_unchanged=2,eligible_families_remaining=len(remaining),original_size=len(original),rebuilt_size=len(rebuilt),jaccard=len(original&rebuilt)/len(union) if union else '',interpretation='Cannot meet unchanged two-family rule with one eligible family. Structural dependence, not evidence of robustness.',remaining_fit_policy='Unchanged independent within-family DE cached; full human pool reselected, never retention-only within old genes'))
        loo.append(dict(branch=branch,omitted_family='human_anchor',minimum_family_support_unchanged=2,eligible_families_remaining=2,original_size=len(original),rebuilt_size='',jaccard='',interpretation='Not estimable under mandatory-human-anchor rule',remaining_fit_policy='No anchor-free redefinition silently substituted'))
    write(root/'stability/family_omission_reconstruction.csv',loo,list(loo[0]))
    # Historical comparison is deliberately calculated after new membership.
    legacy=read(project/'results/tables/refined_crrm_definition/CRRM_core_gene_table.csv')
    comparison=[]
    for old in legacy:
        r={'legacy_gene':old['gene'],'legacy_module':old['module']}
        for branch,records in full_by_branch.items():
            hits=[x for x in records if old['gene'].upper() in x['display_symbol'].upper().split(';')]
            r[branch+'_testable_human_features']=len(hits)
            r[branch+'_human_response']=any(x['human_threshold_pass'] for x in hits)
            r[branch+'_strict_candidate']=any(x['tier']=='strict_cross_system_candidate' for x in hits)
        comparison.append(r)
    write(out/'legacy70_overlap_audit.csv',comparison,list(comparison[0]))
    (out/'candidate_build_summary.json').write_text(json.dumps(summaries,indent=2)+'\n')
    (out/'candidate_freeze_manifest.json').write_text(json.dumps(freezes,indent=2)+'\n')
    print(json.dumps(summaries,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',required=True);p.add_argument('--project',required=True);a=p.parse_args()
    construct(Path(a.root),Path(a.project))
