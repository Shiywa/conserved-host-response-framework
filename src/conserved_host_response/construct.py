from __future__ import annotations

from pathlib import Path

import pandas as pd

from .io import as_bool, load_yaml, read_table, write_table


REQUIRED_EVIDENCE_COLUMNS = {
    "gene",
    "dataset_id",
    "system",
    "direction",
    "detected",
    "primary_deg",
    "supported",
    "module",
}


def construct_framework(evidence_path: str | Path, config_path: str | Path, output_dir: str | Path):
    evidence = read_table(evidence_path)
    missing = sorted(REQUIRED_EVIDENCE_COLUMNS - set(evidence.columns))
    if missing:
        raise ValueError(f"Evidence table is missing columns: {', '.join(missing)}")
    config = load_yaml(config_path)
    rules = config.get("construction", {})
    anchors = [str(x) for x in rules.get("anchor_datasets", [])]
    conservation_systems = [str(x) for x in rules.get("conservation_systems", [])]
    min_systems = int(rules.get("minimum_conservation_systems", 2))
    min_primary_anchors = int(rules.get("minimum_primary_anchor_datasets", 1))
    anchor_presence_column = str(rules.get("anchor_presence_column", "detected"))
    excluded_modules = {str(x).strip().lower() for x in rules.get("excluded_modules", [])}
    if not anchors:
        raise ValueError("construction.anchor_datasets must list at least one dataset_id")
    if not conservation_systems:
        raise ValueError("construction.conservation_systems must list at least one system")
    if anchor_presence_column not in evidence.columns:
        raise ValueError(
            "construction.anchor_presence_column is not present in the evidence table: "
            f"{anchor_presence_column}"
        )

    evidence["gene"] = evidence["gene"].str.strip().str.upper()
    evidence["direction"] = evidence["direction"].str.strip().str.lower()
    evidence["detected_bool"] = as_bool(evidence["detected"])
    evidence["primary_bool"] = as_bool(evidence["primary_deg"])
    evidence["supported_bool"] = as_bool(evidence["supported"])
    evidence["anchor_presence_bool"] = as_bool(evidence[anchor_presence_column])

    rows = []
    for gene, gene_rows in evidence.groupby("gene", sort=True):
        anchor_rows = gene_rows[gene_rows["dataset_id"].isin(anchors)]
        detected_anchors = set(anchor_rows.loc[anchor_rows["anchor_presence_bool"], "dataset_id"])
        anchor_directions = set(
            anchor_rows.loc[
                anchor_rows["anchor_presence_bool"] & anchor_rows["direction"].isin({"up", "down"}), "direction"
            ]
        )
        anchor_complete = all(anchor in detected_anchors for anchor in anchors)
        direction_consistent = anchor_complete and len(anchor_directions) == 1
        consensus_direction = next(iter(anchor_directions)) if direction_consistent else ""
        primary_anchor_count = int(
            anchor_rows.loc[anchor_rows["primary_bool"], "dataset_id"].nunique()
        )
        supported_systems = set(
            gene_rows.loc[
                gene_rows["system"].isin(conservation_systems)
                & gene_rows["supported_bool"]
                & (gene_rows["direction"] == consensus_direction),
                "system",
            ]
        )
        modules = [x for x in gene_rows["module"].str.strip().unique().tolist() if x]
        module = modules[0] if len(modules) == 1 else ""
        module_eligible = bool(module) and module.strip().lower() not in excluded_modules
        selected = (
            anchor_complete
            and direction_consistent
            and primary_anchor_count >= min_primary_anchors
            and len(supported_systems) >= min_systems
            and module_eligible
        )
        failure_reasons = []
        if not anchor_complete:
            failure_reasons.append("anchor_detection")
        if anchor_complete and not direction_consistent:
            failure_reasons.append("anchor_direction")
        if primary_anchor_count < min_primary_anchors:
            failure_reasons.append("anchor_primary_deg")
        if len(supported_systems) < min_systems:
            failure_reasons.append("cross_system_support")
        if not module_eligible:
            failure_reasons.append("module_assignment")
        rows.append(
            {
                "gene": gene,
                "module": module,
                "direction": consensus_direction,
                "anchor_complete": anchor_complete,
                "anchor_presence_column": anchor_presence_column,
                "module_eligible": module_eligible,
                "anchor_direction_consistent": direction_consistent,
                "primary_anchor_datasets": primary_anchor_count,
                "supporting_conservation_systems": len(supported_systems),
                "conservation_system_names": ";".join(sorted(supported_systems)),
                "selected": selected,
                "failure_reasons": ";".join(failure_reasons),
            }
        )

    audit = pd.DataFrame(rows)
    selected = audit[audit["selected"]][["gene", "module", "direction"]].copy()
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    write_table(audit, output_dir / "construction_gene_audit.tsv")
    write_table(selected, output_dir / "selected_gene_table.tsv")
    return audit, selected
