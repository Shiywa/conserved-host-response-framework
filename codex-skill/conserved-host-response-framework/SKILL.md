---
name: conserved-host-response-framework
description: Build, audit, freeze, score, validate, and position conserved host-response frameworks from heterogeneous transcriptomic datasets. Use when designing a multi-dataset or cross-species disease signature, separating construction from validation, checking shared samples or study-family dependence, creating immutable gene/module/scoring manifests, applying a frozen framework to bulk or pseudobulk expression, interpreting severity or cellular localization after freezing, or preparing bounded manuscript claims for such a framework.
---

# Conserved Host-Response Framework

Use the bundled toolkit for deterministic calculations and this workflow for scientific decisions. Treat each disease-specific framework as a biological object; do not treat the workflow as a universally validated discovery algorithm.

## Workflow

1. Define the intended object and claim domain: comparative biology, diagnosis, prediction, prognosis, mechanism, or another use.
2. Build the dataset manifest before gene selection. Read `references/data-schema.md`.
3. Assign every dataset to construction or a named post-freeze role. Record study family, shared controls, biological replicate, species, tissue, endpoint, and source.
4. Run `audit`. Resolve every `FAIL`. Preserve every `WARN` in the scientific report; do not silently relabel related series as independent cohorts.
5. Build direction-aware construction evidence. Run `construct` only after the rules are declared in configuration.
6. Inspect the gene-level selection audit. Require interpretable module assignment and record all failures.
7. Run `freeze` before opening severity, recovery, cellular, cross-virus, or other downstream results. Preserve the generated manifests, metadata, and checksums.
8. Run `score` without changing genes, modules, directions, weights, or missing-gene rules. Report coverage and dataset-relative normalization.
9. Run `validate` using biological replicates. For single-cell resources, create donor/patient pseudobulks before formal inference.
10. Separate demonstrated findings, literature context, and interpretation. Read `references/evidence-levels.md` and `references/claim-guardrails.md`.

## Commands

Run the wrapper from this skill directory:

```bash
python scripts/run_toolkit.py audit --manifest dataset_manifest.tsv --output-dir audit
python scripts/run_toolkit.py construct --evidence evidence_long.tsv --config framework_config.yaml --output-dir construction
python scripts/run_toolkit.py freeze --gene-table selected_gene_table.tsv --config framework_config.yaml --output-dir frozen
python scripts/run_toolkit.py score --expression expression.csv --gene-manifest frozen/frozen_gene_manifest.tsv --scoring-spec frozen/frozen_scoring_spec.yaml --output-dir scores
python scripts/run_toolkit.py validate --scores scores/sample_framework_scores.tsv --metadata metadata.csv --group-column condition --case infected --control control --output-dir validation
```

Use `assets/dataset_manifest_template.tsv` and `assets/framework_config_template.yaml` as starting points.
For a prespecified subset, add `--filter-column COLUMN --filter-value VALUE`; never select a subset after inspecting effects.

## Decision gates

- Stop construction if roles overlap, outcome data enter selection without declaration, gene identities are duplicated, or direction rules are undefined.
- Continue with an explicit warning when construction datasets share a study family or controls; downgrade independence claims.
- Call the object frozen only when manifests and scoring rules are versioned and checksummed before downstream analysis.
- Call validation post-freeze only when no feature selection or coefficient learning uses the validation labels.
- Call clinical results associations unless prospective prediction, calibration, thresholds, and decision utility were evaluated.
- Call cross-species results directional or module-level conservation unless homologous magnitude and pathogenesis were directly established.
- Call cellular results localization unless paired compartments or perturbations establish signaling or mechanism.

## Outputs

Return:

- dataset-role and dependence audit;
- construction gene audit;
- frozen gene/module/scoring manifests and checksums;
- sample and module scores with coverage;
- post-freeze effect estimates;
- a claim-boundary summary stating what the analysis does and does not establish.

Keep negative, heterogeneous, low-coverage, and non-independent results visible.
