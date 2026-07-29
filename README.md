# Conserved Host-Response Framework Toolkit

`conserved-host-response-framework` is a configuration-driven Python toolkit
for constructing and applying a **frozen** transcriptomic host-response
framework across heterogeneous datasets. It was developed alongside a
respiratory syncytial virus (RSV) study, but the software is not limited to
RSV.

## What the toolkit does

The workflow has five explicit stages:

1. **audit** dataset roles, shared controls, related study families, and
   construction/post-freeze separation;
2. **construct** a direction-aware candidate framework from declared evidence;
3. **freeze** gene membership, modules, directions, and scoring rules into
   checksummed manifests;
4. **score** an unchanged framework in a new expression matrix; and
5. **validate** a post-freeze association using independent biological
   replicates or verified pairs.

The package does not establish diagnosis, prediction, calibration, clinical
utility, causality, or universal cross-disease biological validity.

## Installation

```bash
python -m pip install -e .
```

The package requires Python 3.10 or later. The command-line entry point is
`chrf`.

## Minimal workflow

```bash
chrf audit --manifest dataset_manifest.tsv --output-dir audit
chrf construct --evidence evidence_long.tsv --config framework_config.yaml --output-dir construction
chrf freeze --gene-table construction/selected_gene_table.tsv --config framework_config.yaml --output-dir frozen
chrf score --expression expression.csv --gene-manifest frozen/frozen_gene_manifest.tsv \
  --scoring-spec frozen/frozen_scoring_spec.yaml --output-dir scores
chrf validate --scores scores/sample_framework_scores.tsv --metadata metadata.csv \
  --group-column condition --case infected --control control --output-dir validation
```

For a prespecified paired subset, add `--pair-column`, `--filter-column`, and
`--filter-value`. The paired-validation output includes `pairing_audit.tsv` so
complete, incomplete, and duplicated pairs remain visible.

## RSV/CRRM example

`examples/rsv-crrm/` contains a frozen 70-gene, five-module Conserved RSV
Response Module (CRRM) example. It is a validated example of the toolkit, not
the only disease to which the software may be applied.

The bundled Linux acceptance test verifies:

- five construction datasets in four study families;
- declared shared study-family and baseline-control dependence between
  GSE155237 and GSE166161;
- 70 frozen genes in five modules, with 68 construction directions up and two
  down;
- the preserved unsigned, equal-weight sample scoring rule
  (`invert_down_genes: false`);
- 70/70 frozen-gene coverage and numerical reproduction of the GSE231666
  reference total score; and
- 22 exact, prespecified donor-age pairs with a positive paired total effect.

Run the complete bundled example from the repository root:

```bash
PYTHONPATH=src python -m unittest discover -s tests -v
python codex-skill/conserved-host-response-framework/scripts/run_toolkit.py audit \
  --manifest examples/rsv-crrm/dataset_manifest.tsv \
  --output-dir examples/rsv-crrm/outputs/audit
python codex-skill/conserved-host-response-framework/scripts/run_toolkit.py freeze \
  --gene-table examples/rsv-crrm/inputs/crrm_core70_source.csv \
  --config examples/rsv-crrm/framework_config.yaml \
  --output-dir examples/rsv-crrm/outputs/frozen
python codex-skill/conserved-host-response-framework/scripts/run_toolkit.py score \
  --expression examples/rsv-crrm/inputs/GSE231666_core70_gene_zscores.csv \
  --gene-manifest examples/rsv-crrm/outputs/frozen/frozen_gene_manifest.tsv \
  --scoring-spec examples/rsv-crrm/outputs/frozen/frozen_scoring_spec.yaml \
  --output-dir examples/rsv-crrm/outputs/GSE231666_scores
python codex-skill/conserved-host-response-framework/scripts/run_toolkit.py validate \
  --scores examples/rsv-crrm/outputs/GSE231666_scores/sample_framework_scores.tsv \
  --metadata examples/rsv-crrm/inputs/GSE231666_metadata_reference_scores.csv \
  --group-column condition --case RSV --control mock --pair-column pair_key \
  --filter-column primary_pair --filter-value true \
  --output-dir examples/rsv-crrm/outputs/GSE231666_validation
python scripts/verify_rsv_example.py --example-dir examples/rsv-crrm
```

## Construction replay

`outputs/construction_replay/` records a Linux replay of the original frozen
CRRM construction from the original project's **aggregated 571-gene evidence
tables**. With the documented adapter and formal selection rules, the replay
reproduced all 70 frozen genes, directions, and module assignments exactly:

| Measure | Result |
|---|---:|
| Shared membership | 70/70 |
| Missing genes | 0 |
| Additional genes | 0 |
| Direction mismatches | 0 |
| Module mismatches | 0 |

This replay does **not** rerun differential expression or rank integration
from the original 13 contrasts or from raw expression matrices. It validates
the implementation against the existing aggregated pre-freeze evidence.

## Reproducibility and boundaries

- Read `LINUX_VALIDATION_REPORT.md` for the acceptance and replay record.
- Read `DATA_AND_CODE_PROVENANCE.md` for public data provenance and the
  distinction between software validation and biological claims.
- Read `UNRESOLVED_LIMITATIONS.md` before interpreting results or extending a
  framework to another disease.
- The software-level toy test demonstrates disease-independent behavior; it is
  not cross-disease biological validation.

## Citation

Please cite the associated CRRM manuscript when available and this version of
the toolkit. Citation metadata are provided in `CITATION.cff`. A Zenodo DOI is
not yet assigned; versioned GitHub Releases should be used until archival DOI
registration is completed.

## License

MIT. See `LICENSE`.
