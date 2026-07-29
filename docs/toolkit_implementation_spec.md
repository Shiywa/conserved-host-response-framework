# Toolkit implementation specification

## Scientific object

The toolkit distinguishes:

- a disease-specific frozen biological framework;
- a disease-agnostic implementation of construction, freezing, scoring, and validation;
- downstream association from diagnosis, prediction, and causality.

## Invariants

1. Dataset roles are declared before construction.
2. Shared samples, controls, and study families are reported.
3. Clinical severity, recovery, single-cell localization, and cross-virus context do not alter frozen membership unless explicitly declared as construction inputs.
4. Gene, module, direction, weight, normalization, missing-gene, and orthology rules are versioned.
5. Post-freeze analyses never modify the frozen manifest.
6. Dataset-relative scores are not presented as absolute cross-platform values.
7. Negative and heterogeneous validation results remain visible.

## Command contract

| Command | Input | Deterministic output |
|---|---|---|
| `audit` | Dataset manifest | Dependence/leakage findings and summary |
| `construct` | Direction-aware evidence table and configuration | Gene-level rule audit and selected genes |
| `freeze` | Selected/final gene table and configuration | Immutable manifests, metadata, and checksums |
| `score` | Expression matrix and frozen files | Sample/module scores, coverage, and scoring audit |
| `validate` | Scores and metadata | Effect estimates, FDR, and bounded interpretation report |

## Required validation

- Unit tests for dependence, shared controls, leakage, direction consistency, missing genes, duplicated symbols, and paired effects.
- RSV acceptance test for 70 genes, five modules, two down-direction genes, shared-baseline warning, score reproduction, and 22 exact pairs.
- Toy disease example for software-level disease independence.
- A second real disease is required before claiming cross-disease scientific validation.

## Construction adapter note

The required generic evidence columns retain their literal meanings: `detected`
records measurement/detection, `primary_deg` records the declared primary DEG
criterion, and `supported` records eligible same-direction conservation support.
Some legacy projects define anchor eligibility using a distinct direction-aware
evidence field (for example, a 0/1/2 evidence score). In that situation an
adapter may add an explicit Boolean field and set
`construction.anchor_presence_column` to that field. The adapter must preserve
the original fields, document the mapping, and must not relabel directional or
rank support as a primary DEG.

Projects that formally prohibit an unassigned or otherwise non-interpretable
module from entering a frozen core may declare those values through
`construction.excluded_modules`. This is an explicit project rule, not a
software-imposed functional annotation.

## Manuscript boundary

Version 0.1.0 may support a claim of a reusable, configuration-driven implementation. It does not support claims of universal performance, superiority, clinical prediction, or cross-disease validation.
