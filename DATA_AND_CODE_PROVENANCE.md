# Data and code provenance

## Public repository boundary

This repository publishes software, templates, tests, a frozen RSV/CRRM
example, and small derived validation artifacts. It does not publish raw
expression matrices, restricted clinical data, manuscript files, figures, or
the original analysis-server scripts.

The GSE231666 example contains a derived, within-dataset 70-gene z-score matrix
and derived score metadata required for numerical reproducibility. It is not a
raw sequencing or microarray matrix. Sample labels are repository accessions or
study-derived aliases; no direct personal identifiers are included.

## RSV example sources

| Role | Public accession or repository | Public artifact used here |
|---|---|---|
| Human construction field | GSE155237 | Aggregated direction-aware evidence |
| Related human construction field | GSE166161 | Aggregated direction-aware evidence |
| Mouse conservation | GSE197868 | Aggregated direction-aware evidence |
| Cotton rat conservation | GSE297903 | Aggregated direction-aware evidence |
| Human respiratory organoid conservation | GSE263272 | Aggregated direction-aware evidence |
| Post-freeze paired scoring example | GSE231666 | Derived 70-gene z-scores, metadata, and reference score |

The construction replay begins with aggregated evidence from the original
571-gene candidate resource. `scripts/adapt_rsv_server_evidence.py` maps those
pre-freeze summaries to the toolkit schema. It does not reconstruct the
original differential-expression analyses or signed rank integration from the
13 contrasts.

## What was validated

### Software validation

The unit tests cover data-role auditing, direction-aware construction, freezing,
coverage reporting, paired validation, and deterministic artifacts. A toy
non-RSV configuration demonstrates that the implementation itself is not
hard-coded to RSV.

### RSV biological-example validation

The bundled RSV acceptance test reproduces the predeclared frozen definition,
the GSE155237/GSE166161 dependence warnings, GSE231666 numerical score
agreement, and a post-freeze paired association. The construction replay
reproduces 70/70 frozen membership, directions, and modules from aggregated
evidence.

### Claims not established

This repository does not establish diagnostic accuracy, prediction,
calibration, clinical utility, causal mechanisms, RSV exclusivity, or
cross-disease biological generalizability. A second real disease framework is
required before making a cross-disease biological validation claim.
