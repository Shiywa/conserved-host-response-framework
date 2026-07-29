from __future__ import annotations

from pathlib import Path

import pandas as pd

from .io import read_table, write_json, write_table


REQUIRED_DATASET_COLUMNS = {
    "dataset_id",
    "study_family",
    "role",
    "biological_role",
    "shared_control_group",
    "species",
    "tissue",
}


def _finding(check: str, status: str, detail: str) -> dict[str, str]:
    return {"check": check, "status": status, "detail": detail}


def audit_manifest(manifest_path: str | Path, output_dir: str | Path) -> list[dict[str, str]]:
    manifest = read_table(manifest_path)
    missing = sorted(REQUIRED_DATASET_COLUMNS - set(manifest.columns))
    if missing:
        raise ValueError(f"Dataset manifest is missing columns: {', '.join(missing)}")

    for column in REQUIRED_DATASET_COLUMNS:
        manifest[column] = manifest[column].astype(str).str.strip()
    manifest["role_norm"] = manifest["role"].str.lower()
    manifest["bio_norm"] = manifest["biological_role"].str.lower()
    construction = manifest[manifest["role_norm"] == "construction"].copy()
    nonconstruction = manifest[manifest["role_norm"] != "construction"].copy()
    findings: list[dict[str, str]] = []

    duplicate_ids = manifest.loc[manifest["dataset_id"].duplicated(False), "dataset_id"].unique().tolist()
    findings.append(
        _finding(
            "Unique dataset identifiers",
            "PASS" if not duplicate_ids else "FAIL",
            "All dataset_id values are unique." if not duplicate_ids else f"Duplicated: {', '.join(duplicate_ids)}",
        )
    )

    findings.append(
        _finding(
            "Construction datasets present",
            "PASS" if len(construction) else "FAIL",
            f"{len(construction)} construction datasets.",
        )
    )
    findings.append(
        _finding(
            "Post-freeze datasets present",
            "PASS" if len(nonconstruction) else "WARN",
            f"{len(nonconstruction)} non-construction datasets.",
        )
    )

    disallowed_terms = (
        "severity", "clinical", "outcome", "hospital", "single_cell", "single-cell",
        "cross_virus", "cross-virus", "recovery",
    )
    disallowed = construction[
        construction["bio_norm"].apply(lambda value: any(term in value for term in disallowed_terms))
    ]
    findings.append(
        _finding(
            "Outcome and cellular data excluded from construction",
            "PASS" if disallowed.empty else "FAIL",
            "No downstream outcome/cellular roles enter construction."
            if disallowed.empty
            else "Construction includes downstream roles: "
            + ", ".join(f"{row.dataset_id} ({row.biological_role})" for row in disallowed.itertuples()),
        )
    )

    family_groups = (
        construction[construction["study_family"] != ""]
        .groupby("study_family")["dataset_id"]
        .apply(list)
    )
    repeated_families = {family: ids for family, ids in family_groups.items() if len(ids) > 1}
    findings.append(
        _finding(
            "Construction study-family independence",
            "PASS" if not repeated_families else "WARN",
            "Each construction dataset belongs to a distinct study family."
            if not repeated_families
            else "; ".join(f"{family}: {', '.join(ids)}" for family, ids in repeated_families.items()),
        )
    )

    shared_groups = (
        construction[construction["shared_control_group"] != ""]
        .groupby("shared_control_group")["dataset_id"]
        .apply(list)
    )
    shared = {group: ids for group, ids in shared_groups.items() if len(ids) > 1}
    findings.append(
        _finding(
            "Shared construction controls",
            "PASS" if not shared else "WARN",
            "No shared control groups are declared among construction datasets."
            if not shared
            else "; ".join(f"{group}: {', '.join(ids)}" for group, ids in shared.items()),
        )
    )

    construction_ids = set(construction["dataset_id"])
    validation_ids = set(nonconstruction["dataset_id"])
    overlap = sorted(construction_ids & validation_ids)
    findings.append(
        _finding(
            "Construction/post-freeze dataset separation",
            "PASS" if not overlap else "FAIL",
            "No dataset identifier appears in both roles." if not overlap else f"Overlap: {', '.join(overlap)}",
        )
    )

    missing_sources = manifest[manifest.get("source_reference", pd.Series([""] * len(manifest))).astype(str).str.strip() == ""]
    findings.append(
        _finding(
            "Dataset source provenance",
            "PASS" if missing_sources.empty else "WARN",
            "Every dataset has a source reference or repository identifier."
            if missing_sources.empty
            else f"{len(missing_sources)} dataset rows lack source_reference.",
        )
    )

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    findings_frame = pd.DataFrame(findings)
    write_table(findings_frame, output_dir / "manifest_audit.tsv")
    write_json(
        {
            "n_datasets": int(len(manifest)),
            "n_construction_datasets": int(len(construction)),
            "n_construction_study_families": int(construction["study_family"].replace("", pd.NA).nunique()),
            "n_postfreeze_datasets": int(len(nonconstruction)),
            "status_counts": findings_frame["status"].value_counts().to_dict(),
        },
        output_dir / "manifest_audit_summary.json",
    )

    report = [
        "# Dataset-role and dependence audit",
        "",
        f"- Construction datasets: {len(construction)}",
        f"- Construction study families: {construction['study_family'].replace('', pd.NA).nunique()}",
        f"- Post-freeze datasets: {len(nonconstruction)}",
        "",
        "| Check | Status | Detail |",
        "|---|---|---|",
    ]
    report.extend(
        f"| {row['check']} | {row['status']} | {row['detail'].replace('|', '/')} |" for row in findings
    )
    report.extend(
        [
            "",
            "Warnings disclose dependence or incomplete provenance; they do not automatically invalidate a framework.",
            "Failures must be resolved before construction or post-freeze validation is interpreted.",
        ]
    )
    (output_dir / "manifest_audit.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    return findings
