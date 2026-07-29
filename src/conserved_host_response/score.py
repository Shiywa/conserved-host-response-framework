from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd

from .io import load_yaml, read_numeric_table, read_table, write_json, write_table


def _safe_column(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", value).strip("_")


def score_framework(
    expression_path: str | Path,
    gene_manifest_path: str | Path,
    scoring_spec_path: str | Path,
    output_dir: str | Path,
    expression_gene_column: str = "gene",
):
    expression = read_numeric_table(expression_path)
    if expression_gene_column not in expression.columns:
        raise ValueError(f"Expression table lacks gene column: {expression_gene_column}")
    sample_columns = [column for column in expression.columns if column != expression_gene_column]
    if len(sample_columns) < 2:
        raise ValueError("At least two sample columns are required for within-dataset z scoring")
    expression[expression_gene_column] = expression[expression_gene_column].astype(str).str.strip().str.upper()
    for column in sample_columns:
        expression[column] = pd.to_numeric(expression[column], errors="coerce")
    expression = expression[expression[expression_gene_column] != ""]
    collapsed = expression.groupby(expression_gene_column, sort=False)[sample_columns].mean()

    manifest = read_table(gene_manifest_path)
    required = {"gene", "module", "direction"}
    if not required.issubset(manifest.columns):
        raise ValueError(f"Frozen manifest must include: {', '.join(sorted(required))}")
    manifest["gene"] = manifest["gene"].str.strip().str.upper()
    scoring = load_yaml(scoring_spec_path)
    invert_down = bool(scoring.get("invert_down_genes", False))
    ddof = int(scoring.get("zscore_ddof", 1))

    available = [gene for gene in manifest["gene"] if gene in collapsed.index]
    values = collapsed.loc[available].to_numpy(dtype=float)
    means = np.nanmean(values, axis=1, keepdims=True)
    sds = np.nanstd(values, axis=1, ddof=ddof, keepdims=True)
    variable = np.isfinite(sds[:, 0]) & (sds[:, 0] > 0)
    available_variable = [gene for gene, keep in zip(available, variable) if keep]
    z = (values[variable] - means[variable]) / sds[variable]
    z_frame = pd.DataFrame(z, index=available_variable, columns=sample_columns)
    if invert_down:
        directions = manifest.set_index("gene")["direction"]
        down = [gene for gene in available_variable if directions.get(gene) == "down"]
        z_frame.loc[down] = -z_frame.loc[down]

    scores = pd.DataFrame({"sample_id": sample_columns})
    scores["framework_total"] = z_frame.mean(axis=0, skipna=True).reindex(sample_columns).to_numpy()
    coverage_rows = [
        {
            "component": "framework_total",
            "genes_total": int(len(manifest)),
            "genes_detected": int(len(available)),
            "genes_scored": int(len(available_variable)),
            "coverage_fraction": len(available_variable) / len(manifest) if len(manifest) else np.nan,
        }
    ]
    for module, module_rows in manifest.groupby("module", sort=False):
        genes = module_rows["gene"].tolist()
        scored = [gene for gene in genes if gene in z_frame.index]
        column = f"module__{_safe_column(module)}"
        scores[column] = (
            z_frame.loc[scored].mean(axis=0, skipna=True).reindex(sample_columns).to_numpy()
            if scored
            else np.nan
        )
        coverage_rows.append(
            {
                "component": module,
                "genes_total": len(genes),
                "genes_detected": sum(gene in collapsed.index for gene in genes),
                "genes_scored": len(scored),
                "coverage_fraction": len(scored) / len(genes) if genes else np.nan,
            }
        )

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    write_table(scores, output_dir / "sample_framework_scores.tsv")
    write_table(pd.DataFrame(coverage_rows), output_dir / "framework_gene_coverage.tsv")
    write_table(
        z_frame.reset_index().rename(columns={"index": "gene"}),
        output_dir / "frozen_gene_zscores.tsv",
    )
    write_json(
        {
            "expression_gene_rows": int(len(collapsed)),
            "sample_count": len(sample_columns),
            "frozen_gene_count": int(len(manifest)),
            "detected_frozen_genes": len(available),
            "variable_scored_genes": len(available_variable),
            "invert_down_genes": invert_down,
            "zscore_ddof": ddof,
            "score_interpretation": "relative within-dataset score",
        },
        output_dir / "scoring_audit.json",
    )
    return scores, pd.DataFrame(coverage_rows), z_frame

