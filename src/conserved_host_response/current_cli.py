"""CLI for the current CRRM388 paper. `chrf` remains the explicit legacy CLI."""
from __future__ import annotations
import argparse,json
from pathlib import Path
import pandas as pd
from .crrm388 import reconstruct,verify_definition
from .reference import FrozenReference,score
from .study import study_directory,prepare,clinical,figures,run_native,export_tables,study_input_path

def main(argv=None):
    parser=argparse.ArgumentParser(prog='crrm388',description='Frozen signed CRRM388 study workflows')
    parser.add_argument('--study-dir')
    sub=parser.add_subparsers(dest='command',required=True)
    p=sub.add_parser('definition');p.add_argument('--output-dir',required=True);p.add_argument('--evidence')
    sub.add_parser('check-definition')
    p=sub.add_parser('fit-reference');p.add_argument('--expression',required=True);p.add_argument('--reference-samples',required=True);p.add_argument('--output',required=True);p.add_argument('--units');p.add_argument('--sd-threshold',type=float,default=1e-8);p.add_argument('--nonzero-support',help='Initial-human raw-count support TSV:human_ensembl,nonzero;not inferred from logCPM')
    p=sub.add_parser('score');p.add_argument('--expression',required=True);p.add_argument('--reference',required=True);p.add_argument('--output-dir',required=True)
    p=sub.add_parser('prepare');p.add_argument('--workspace',required=True);p.add_argument('--data-root');p.add_argument('--resource-root');p.add_argument('--input-bundle')
    p=sub.add_parser('clinical');p.add_argument('--workspace',required=True);p.add_argument('--use-saved-splines',action='store_true');p.add_argument('--rscript',default='Rscript')
    p=sub.add_parser('figures');p.add_argument('--workspace',required=True)
    p=sub.add_parser('export-tables');p.add_argument('--workspace',required=True);p.add_argument('--output-dir',required=True)
    p=sub.add_parser('native');p.add_argument('--workspace',required=True);p.add_argument('--stage',required=True);p.add_argument('--script',required=True);p.add_argument('--rscript',default='Rscript');p.add_argument('--argument',action='append',default=[])
    p=sub.add_parser('verify');p.add_argument('--workspace');p.add_argument('--output',required=True)
    args=parser.parse_args(argv);bundle=study_directory(args.study_dir)
    genes=pd.read_csv(bundle/'frozen/crrm388_genes.tsv',sep='\t');sets=pd.read_csv(bundle/'frozen/crrm388_subsets.tsv',sep='\t')
    if args.command=='check-definition':
        print(json.dumps(verify_definition(genes,sets),indent=2))
    elif args.command=='definition':
        d=reconstruct(pd.read_csv(Path(args.evidence) if args.evidence else study_input_path('evidence/common_tested_discovery.tsv.gz'),sep='\t'))
        expected=genes.set_index('human_ensembl').frozen_direction;actual=d.set_index('human_ensembl').frozen_direction
        if not actual.sort_index().equals(expected.sort_index()):raise AssertionError('Reconstructed frozen definition differs')
        out=Path(args.output_dir);out.mkdir(parents=True,exist_ok=True);d.to_csv(out/'reconstructed_crrm388.tsv',sep='\t',index=False)
        (out/'definition_checks.json').write_text(json.dumps(verify_definition(genes,sets),indent=2)+'\n')
    elif args.command=='fit-reference':
        x=pd.read_csv(args.expression,index_col=0);ids=Path(args.reference_samples).read_text().splitlines();units=pd.read_csv(args.units,sep='\t').set_index('sample_id').unit if args.units else None
        support=pd.read_csv(args.nonzero_support,sep='\t').set_index('human_ensembl').nonzero.astype(str).str.lower().isin(['true','1']) if args.nonzero_support else None
        FrozenReference.fit(x,ids,units=units,sd_threshold=args.sd_threshold,nonzero_support=support).save(args.output)
    elif args.command=='score':
        x=pd.read_csv(args.expression,index_col=0);ref=FrozenReference.load(args.reference);subsets={n:d.human_ensembl.tolist() for n,d in sets.groupby('set_name',sort=False)}
        result,coverage,z=score(x,ref,genes.set_index('human_ensembl'),subsets)
        out=Path(args.output_dir);out.mkdir(parents=True,exist_ok=True);result.to_csv(out/'scores.tsv',sep='\t');coverage.to_csv(out/'coverage.tsv',sep='\t',index=False)
    elif args.command=='prepare':prepare(args.workspace,study_dir=bundle,data_root=args.data_root,resource_root=args.resource_root,input_bundle=args.input_bundle)
    elif args.command=='clinical':clinical(args.workspace,regenerate_splines=not args.use_saved_splines,rscript=args.rscript)
    elif args.command=='figures':figures(args.workspace)
    elif args.command=='export-tables':export_tables(args.workspace,args.output_dir,study_dir=bundle)
    elif args.command=='native':run_native(args.workspace,args.stage,args.script,rscript=args.rscript,arguments=args.argument)
    elif args.command=='verify':
        from .replay import verify
        result=verify(bundle,args.workspace);Path(args.output).parent.mkdir(parents=True,exist_ok=True);Path(args.output).write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

if __name__=='__main__':main()
