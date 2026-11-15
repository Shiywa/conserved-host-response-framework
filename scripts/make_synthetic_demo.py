"""Generate an explicitly synthetic CRRM388 CLI example; no biological data."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from conserved_host_response.study import study_directory


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output-dir',required=True)
    a=p.parse_args();out=Path(a.output_dir);out.mkdir(parents=True,exist_ok=True)
    genes=pd.read_csv(study_directory()/'frozen/crrm388_genes.tsv',sep='\t')
    values=np.random.default_rng(388).normal(size=(388,8))
    sign=genes.frozen_direction.map({'up':1.,'down':-1.}).to_numpy()
    values[:,4:]+=1.5*sign[:,None]
    ids=[f'synthetic_{i}' for i in range(8)]
    pd.DataFrame(values,index=genes.human_ensembl,columns=ids).to_csv(out/'synthetic_expression.csv')
    (out/'synthetic_reference_samples.txt').write_text('\n'.join(ids[:4])+'\n')
    (out/'README.txt').write_text('Synthetic software demonstration only; not experimental or clinical evidence.\n')


if __name__=='__main__':main()
