#!/usr/bin/env python3
"""Verify the Linux-generated RSV example against frozen manuscript source tables."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--example-dir", default="examples/rsv-crrm")
    args = parser.parse_args()
    root = Path(args.example_dir).resolve()
    outputs = root / "outputs"
    checks = []

    def record(check, passed, detail):
        checks.append({"check": check, "status": "PASS" if passed else "FAIL", "detail": detail})

    audit = pd.read_csv(outputs / "audit/manifest_audit.tsv", sep="\t")
    family = audit.loc[audit["check"] == "Construction study-family independence"].iloc[0]
    shared = audit.loc[audit["check"] == "Shared construction controls"].iloc[0]
    record("Related human study family detected", family["status"] == "WARN", family["detail"])
    record("Shared baseline detected", shared["status"] == "WARN", shared["detail"])

    metadata = json.loads((outputs / "frozen/framework_metadata.json").read_text(encoding="utf-8"))
    record("Frozen gene count", metadata["gene_count"] == 70, str(metadata["gene_count"]))
    record("Frozen module count", metadata["module_count"] == 5, str(metadata["module_count"]))
    record(
        "Frozen direction counts",
        metadata["direction_counts"].get("up") == 68 and metadata["direction_counts"].get("down") == 2,
        json.dumps(metadata["direction_counts"], sort_keys=True),
    )

    coverage = pd.read_csv(outputs / "GSE231666_scores/framework_gene_coverage.tsv", sep="\t")
    total_coverage = coverage.loc[coverage["component"] == "framework_total"].iloc[0]
    record(
        "GSE231666 frozen-gene scoring coverage",
        int(total_coverage["genes_scored"]) == 70,
        f"{int(total_coverage['genes_scored'])}/70",
    )

    scores = pd.read_csv(outputs / "GSE231666_scores/sample_framework_scores.tsv", sep="\t")
    reference = pd.read_csv(root / "inputs/GSE231666_metadata_reference_scores.csv")
    comparison = reference[["sample_id", "CRRM_score"]].merge(
        scores[["sample_id", "framework_total"]], on="sample_id", how="inner", validate="one_to_one"
    )
    max_error = float(np.max(np.abs(comparison["CRRM_score"] - comparison["framework_total"])))
    record("GSE231666 total-score numerical reproduction", max_error < 1e-10, f"max_abs_error={max_error:.3g}")

    validation = pd.read_csv(outputs / "GSE231666_validation/validation_effects.tsv", sep="\t")
    total = validation.loc[validation["score"] == "framework_total"].iloc[0]
    record(
        "GSE231666 exact-pair count",
        int(total["n_case"]) == 22 and int(total["n_control"]) == 22,
        f"{int(total['n_case'])} pairs",
    )
    record("GSE231666 paired total effect direction", float(total["effect"]) > 0, f"effect={float(total['effect']):.6g}")

    frame = pd.DataFrame(checks)
    frame.to_csv(outputs / "RSV_acceptance_checks.tsv", sep="\t", index=False)
    report = [
        "# RSV/CRRM Linux acceptance report",
        "",
        "| Check | Status | Detail |",
        "|---|---|---|",
    ]
    report.extend(
        f"| {row.check} | {row.status} | {str(row.detail).replace('|', '/')} |" for row in frame.itertuples()
    )
    report.extend(
        [
            "",
            "This acceptance test validates role auditing, freezing, dataset-relative scoring, and paired post-freeze reuse.",
            "It does not validate a complete 571-to-70 construction replay unless the original server-side pre-freeze evidence is adapted and tested separately.",
        ]
    )
    (outputs / "RSV_ACCEPTANCE_REPORT.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    if (frame["status"] == "FAIL").any():
        print(frame.to_string(index=False))
        raise SystemExit(1)
    print(frame.to_string(index=False))


if __name__ == "__main__":
    main()

