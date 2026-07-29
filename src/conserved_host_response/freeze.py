from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from .io import (
    canonical_json,
    load_yaml,
    read_table,
    sha256_file,
    sha256_text,
    write_json,
    write_table,
    write_yaml,
)


def freeze_framework(gene_table_path: str | Path, config_path: str | Path, output_dir: str | Path):
    source = read_table(gene_table_path)
    config = load_yaml(config_path)
    columns = config.get("columns", {})
    gene_column = columns.get("gene", "gene")
    module_column = columns.get("module", "module")
    direction_column = columns.get("direction", "direction")
    required = [gene_column, module_column, direction_column]
    missing = [column for column in required if column not in source.columns]
    if missing:
        raise ValueError(f"Gene table is missing columns: {', '.join(missing)}")

    frozen = pd.DataFrame(
        {
            "gene": source[gene_column].astype(str).str.strip().str.upper(),
            "module": source[module_column].astype(str).str.strip(),
            "direction": source[direction_column].astype(str).str.strip().str.lower(),
            "source_row": range(1, len(source) + 1),
        }
    )
    if (frozen["gene"] == "").any() or (frozen["module"] == "").any():
        raise ValueError("Frozen genes and modules must be non-empty")
    duplicates = frozen.loc[frozen["gene"].duplicated(False), "gene"].unique().tolist()
    if duplicates:
        raise ValueError(f"Frozen gene symbols must be unique: {', '.join(duplicates)}")
    invalid_directions = sorted(set(frozen["direction"]) - {"up", "down", "neutral"})
    if invalid_directions:
        raise ValueError(f"Invalid directions: {', '.join(invalid_directions)}")

    modules = (
        frozen.groupby("module", sort=False)
        .agg(genes_total=("gene", "size"), first_source_row=("source_row", "min"))
        .reset_index()
        .sort_values("first_source_row")
        .drop(columns="first_source_row")
    )
    scoring = config.get("scoring", {})
    framework = config.get("framework", {})
    canonical_definition = {
        "framework": framework,
        "scoring": scoring,
        "genes": frozen[["gene", "module", "direction"]].to_dict(orient="records"),
    }
    framework_hash = sha256_text(canonical_json(canonical_definition))
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    gene_path = write_table(frozen, output_dir / "frozen_gene_manifest.tsv")
    module_path = write_table(modules, output_dir / "frozen_module_manifest.tsv")
    score_path = write_yaml(scoring, output_dir / "frozen_scoring_spec.yaml")
    metadata = {
        "framework_name": framework.get("name", "unnamed-framework"),
        "disease": framework.get("disease", ""),
        "version": str(framework.get("version", "0.0.0")),
        # Frozen objects are content-addressed. A wall-clock timestamp would make
        # otherwise identical manifests and checksum files non-deterministic.
        "creation_policy": "deterministic_content_addressed",
        "gene_count": int(len(frozen)),
        "module_count": int(len(modules)),
        "direction_counts": frozen["direction"].value_counts().to_dict(),
        "input_gene_table_sha256": sha256_file(gene_table_path),
        "input_config_sha256": sha256_file(config_path),
        "framework_sha256": framework_hash,
    }
    metadata_path = write_json(metadata, output_dir / "framework_metadata.json")
    checksum_paths = [gene_path, module_path, score_path, metadata_path]
    checksum_text = "\n".join(f"{sha256_file(path)}  {path.name}" for path in checksum_paths) + "\n"
    (output_dir / "SHA256SUMS").write_text(checksum_text, encoding="utf-8")
    return frozen, modules, metadata
