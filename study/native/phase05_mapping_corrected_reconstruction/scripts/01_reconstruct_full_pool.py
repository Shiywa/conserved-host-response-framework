#!/usr/bin/env python3
"""Reconstruct an explicitly exploratory cross-model consensus from full DE pools after the phase04 annotation correction."""
import argparse,csv,hashlib,json,math
from pathlib import Path
from importlib import import_module
u=import_module('04_construct_candidates')
read,write,index,num=u.read,u.write,u.index,u.num

def main(root,project):
    cfg=json.loads((root/'parameters/experimental_consensus_v2.json').read_text())
    receipt=json.loads((root/'parameters/mapping_correction_declaration_receipt.json').read_text())
    assert hashlib.sha256((root/'parameters/experimental_consensus_v2.json').read_bytes()).hexdigest()==receipt['config_sha256']
    output=root/'construction/experimental_consensus';output.mkdir(exist_ok=True)
    mapping={
      'A549':('A549_primary__infected_vs_mock__edgeR_QLF.csv','A549_feature_map.csv',None),
      'mouse':('M_lung_primary__infected_vs_control__edgeR_QLF.csv','mouse_feature_map.csv','eligible_mapping'),
      'cotton_rat':('R_D1_primary__age_average__edgeR_QLF.csv','cotton_rat_feature_map.csv','eligible_mapping'),
    }
    families={};symbols={}
    for family,(de,mapfile,eligibility) in mapping.items():
        stats=index(root/'de'/de);lookup={}
        for record in read(root/'mapping'/mapfile):
            if eligibility and not u.flag(record[eligibility]):continue
            ensembl=record.get('human_ensembl','');feature=record['feature_id']
            if ensembl and feature in stats:
                assert ensembl not in lookup,(family,ensembl)
                lookup[ensembl]=stats[feature]
                if family=='A549':symbols[ensembl]=record.get('display_symbol','')
        families[family]=lookup
    for r in read(root/'mapping/human_feature_map.csv'):
        if r['display_symbol']:symbols[r['human_ensembl']]=r['display_symbol']
    def passing(row):
        return row is not None and num(row.get('FDR'))<cfg['per_family_q_lt'] and abs(num(row.get('logFC')))>cfg['per_family_abs_log2FC_gt']
    universe=sorted(set().union(*(set(v) for v in families.values())))
    def build(eligible_families):
        results=[]
        for gene in universe:
            up=[];down=[];tested=[];row=dict(human_ensembl=gene,display_symbol=symbols.get(gene,''))
            for family in mapping:
                stat=families[family].get(gene);present=stat is not None
                if family in eligible_families and present:
                    tested.append(family)
                    if passing(stat):(up if num(stat['logFC'])>0 else down).append(family)
                row.update({family+'_mapped_and_tested':present,family+'_log2FC':stat['logFC'] if present else '',family+'_FDR':stat['FDR'] if present else ''})
            selected=(len(up)>=2 and not down) or (len(down)>=2 and not up)
            direction='up' if len(up)>=2 and not down else 'down' if len(down)>=2 and not up else ''
            row.update(eligible_families_tested=len(tested),n_up_supported=len(up),n_down_supported=len(down),up_support_families=';'.join(up),down_support_families=';'.join(down),candidate=selected,frozen_direction=direction,
                       interpretation='Exploratory experimental consensus; source-unit and cotton-rat bridge qualifications apply')
            results.append(row)
        return results
    full=build(list(mapping))
    full.sort(key=lambda r:(not r['candidate'],-max(r['n_up_supported'],r['n_down_supported']),r['human_ensembl']))
    core=[r for r in full if r['candidate']]
    all3=[r for r in core if max(r['n_up_supported'],r['n_down_supported'])==3]
    write(output/'all_mapped_gene_evidence.csv',full,list(full[0]))
    write(output/'experimental_consensus_candidate.csv',core,list(full[0]))
    write(output/'all_three_families_supported.csv',all3,list(full[0]))
    conflicts=[r for r in full if (r['n_up_supported']>=2 and r['n_down_supported']) or (r['n_down_supported']>=2 and r['n_up_supported'])]
    write(output/'direction_conflict_audit.csv',conflicts,list(full[0]))
    original={r['human_ensembl'] for r in core}
    loo=[]
    for excluded in mapping:
        rebuilt=build([f for f in mapping if f!=excluded])
        chosen=[r for r in rebuilt if r['candidate']]
        genes={r['human_ensembl'] for r in chosen};union=original|genes
        write(output/f'leave_family_out_{excluded}.csv',chosen,list(full[0]))
        loo.append(dict(omitted_family=excluded,required_support=2,original_size=len(original),rebuilt_size=len(genes),retained=len(original&genes),gained=len(genes-original),lost=len(original-genes),jaccard=len(original&genes)/len(union) if union else '',method='Full mapped universe reselected using independent cached family DE; unchanged two-family threshold'))
    write(root/'stability/experimental_family_omission.csv',loo,list(loo[0]))
    # Method sensitivity; membership stays fixed to the primary edgeR results.
    sensitivity=[]
    vfiles={'A549':'A549_primary__infected_vs_mock__limma_voom.csv','mouse':'M_lung_primary__infected_vs_control__limma_voom.csv','cotton_rat':'R_D1_primary__age_average__limma_voom.csv'}
    maps={family:read(root/'mapping'/spec[1]) for family,spec in mapping.items()}
    voom={family:index(root/'de'/file) for family,file in vfiles.items()}
    by_human={}
    for family in mapping:
        by_human[family]={r['human_ensembl']:r['feature_id'] for r in maps[family] if r.get('human_ensembl') and (mapping[family][2] is None or u.flag(r.get(mapping[family][2])))}
    for r in core:
        row={'human_ensembl':r['human_ensembl'],'display_symbol':r['display_symbol'],'frozen_direction':r['frozen_direction']}
        supporting=[]
        for family in mapping:
            feature=by_human[family].get(r['human_ensembl']);v=voom[family].get(feature)
            same=v is not None and passing(v) and ((num(v['logFC'])>0)==(r['frozen_direction']=='up'))
            row[family+'_voom_FDR']=v['FDR'] if v else '';row[family+'_voom_same_direction_support']=same
            if same:supporting.append(family)
        row['voom_supporting_families']=len(supporting);sensitivity.append(row)
    scol=['human_ensembl','display_symbol','frozen_direction']+[f+s for f in mapping for s in ['_voom_FDR','_voom_same_direction_support']]+['voom_supporting_families']
    write(output/'fixed_candidate_method_sensitivity.csv',sensitivity,scol)
    legacy=read(project/'results/tables/refined_crrm_definition/CRRM_core_gene_table.csv')
    old={r['gene'].upper() for r in legacy}
    overlap=[r for r in core if any(x in old for x in r['display_symbol'].upper().split(';'))]
    write(output/'legacy70_overlap.csv',overlap,list(full[0]))
    summary=dict(version='experimental_consensus_mapping_corrected_v2',n_mapped_gene_universe=len(universe),n_candidate_genes=len(core),n_up=sum(r['frozen_direction']=='up' for r in core),n_down=sum(r['frozen_direction']=='down' for r in core),n_all_three_supported=len(all3),n_conflicts_excluded=len(conflicts),n_with_at_least_two_voom_support=sum(r['voom_supporting_families']>=2 for r in sensitivity),legacy70_overlap=len(overlap),n_modules='not defined',experimental_family_omission=loo,
                 status='Exploratory candidate, frozen for subsequent checking; not an externally validated final CRRM',replication_limit='A549 has three deposited culture profiles per group; independent experiment provenance remains qualified. Cotton-rat mapping is an author-annotated bridge, not a proven genome-wide orthology map.',human_primary_negative_preserved=True,external_validation_used=False,single_cell_used=False)
    (output/'build_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    artifacts=['experimental_consensus_candidate.csv','all_three_families_supported.csv']
    freeze=dict(version=summary['version'],n_genes=len(core),scoring='Signed equal-weight mean gene z-score, common supported direction frozen; dataset-relative; missing and invariant genes reported',modules='Not defined; no imposed five-module split',config_sha256=receipt['config_sha256'],membership_files=[dict(file=f,sha256=hashlib.sha256((output/f).read_bytes()).hexdigest()) for f in artifacts],status=summary['status'],validation_reserved=True)
    (output/'candidate_freeze_manifest.json').write_text(json.dumps(freeze,indent=2)+'\n')
    print(json.dumps(summary,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',required=True);p.add_argument('--project',required=True);a=p.parse_args();main(Path(a.root),Path(a.project))
