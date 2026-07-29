#!/usr/bin/env python3
"""Create a toolkit-local, auditable adapter for the original RSV evidence.

This script reads the original project's frozen/refined tables without changing
them. It emits the generic long evidence schema plus an explicit anchor-evidence
field required to replay the original formal 70-gene rule faithfully.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


SOURCES = (
    ("GSE155237", "human_anchor", "human_discovery", "human_evidence"),
    ("GSE166161", "human_anchor", "human_replication", "human_replication_evidence"),
    ("GSE197868", "mouse", "mouse", "mouse_evidence"),
    ("GSE297903", "cotton_rat", "cotton_rat", "cotton_rat_evidence"),
    ("GSE263272", "organoid", "organoid", "organoid_evidence"),
)


def as_flag(series: pd.Series) -> pd.Series:
    return series.fillna(False).astype(bool)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    project = Path(args.project_root).resolve()
    output = Path(args.output_dir).resolve()
    audit_path = project / "results/tables/frozen_crrm_definition/CRRM_gene_table_full_audit.csv"
    refined_path = project / "results/tables/refined_crrm_definition/CRRM_refined_gene_table_full.csv"
    for path in (audit_path, refined_path):
        if not path.is_file():
            raise FileNotFoundError(path)

    audit = pd.read_csv(audit_path, keep_default_na=False)
    refined = pd.read_csv(refined_path, keep_default_na=False)
    modules = refined[["gene", "module", "refined_confidence_level"]].copy()
    if modules["gene"].duplicated().any():
        raise ValueError("Original refined table has duplicated gene symbols")
    table = audit.merge(modules, on="gene", how="inner", validate="one_to_one")
    if len(table) != len(refined):
        raise ValueError("Audit/refined tables do not provide a one-to-one 571-gene universe")

    rows: list[pd.DataFrame] = []
    for dataset_id, system, prefix, evidence_column in SOURCES:
        required = [
            f"{prefix}_detected",
            f"{prefix}_n_primary_deg",
            f"{prefix}_direction",
            evidence_column,
        ]
        missing = [column for column in required if column not in table.columns]
        if missing:
            raise ValueError(f"Original table lacks fields for {dataset_id}: {', '.join(missing)}")
        rows.append(
            pd.DataFrame(
                {
                    "gene": table["gene"].astype(str).str.strip().str.upper(),
                    "dataset_id": dataset_id,
                    "system": system,
                    "direction": table[f"{prefix}_direction"].astype(str).str.strip().str.lower(),
                    "detected": as_flag(table[f"{prefix}_detected"]),
                    "primary_deg": pd.to_numeric(table[f"{prefix}_n_primary_deg"], errors="coerce").fillna(0).gt(0),
                    # Original 0/1/2 evidence records same-direction evidence.
                    "supported": pd.to_numeric(table[evidence_column], errors="coerce").fillna(0).gt(0),
                    "module": table["module"].astype(str).str.strip(),
                    "source_evidence_positive": pd.to_numeric(
                        table[evidence_column], errors="coerce"
                    ).fillna(0).gt(0),
                    "original_evidence_points": pd.to_numeric(
                        table[evidence_column], errors="coerce"
                    ).fillna(0).astype(int),
                    "original_refined_confidence_level": table["refined_confidence_level"].astype(str),
                }
            )
        )

    evidence = pd.concat(rows, ignore_index=True)
    if evidence.duplicated(["gene", "dataset_id"]).any():
        raise ValueError("Adapter produced duplicated gene/dataset records")
    output.mkdir(parents=True, exist_ok=True)
    evidence.to_csv(output / "rsv_server_evidence_generic.tsv", sep="\t", index=False)

    mapping = pd.DataFrame(
        [
            ("gene", "Original refined/audit gene symbol, uppercased."),
            ("dataset_id", "Original source field; GSE155237/GSE166161 remain separate fields within one declared study family."),
            ("system", "human_anchor, mouse, cotton_rat, or organoid."),
            ("direction", "Original source-specific collapsed direction."),
            ("detected", "Original source-specific detection field; retained without reinterpretation."),
            ("primary_deg", "Original source-specific number of primary DE contrasts > 0 (FDR < 0.05 and abs(log2FC) > 0.58)."),
            ("supported", "Original source-specific 0/1/2 evidence > 0; conservation use is additionally direction-checked by construct."),
            ("module", "Original pre-freeze refined module assignment; no GO reassignment is performed by this adapter."),
            ("source_evidence_positive", "Original source-specific 0/1/2 evidence > 0. Used only as the declared anchor-presence rule because original core eligibility required positive evidence in both human fields, not primary DEG status in both."),
        ],
        columns=["field", "mapping"],
    )
    mapping.to_csv(output / "rsv_server_evidence_adapter_mapping.csv", index=False)
    (output / "README.md").write_text(
        "# RSV server-evidence adapter\n\n"
        "This directory is generated inside the toolkit. It reads the original project tables as read-only. "
        "The adapter does not redefine a gene, module, direction, or threshold.\n\n"
        "`source_evidence_positive` is an explicit adapter field, not a replacement for `primary_deg`. "
        "It faithfully represents the original formal requirement that both related human fields have positive direction-aware evidence. "
        "GSE155237 and GSE166161 remain related fields in one declared study family and must not be described as independent cohorts.\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
