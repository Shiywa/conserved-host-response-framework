# Current manuscript workflow and input contract

This code release keeps matrices, complete result tables and workbooks outside
Git. `crrm388 prepare --workspace replay_run/study` creates the native-code and
configuration scaffold. Supply public data through `--data-root` and optional
curated resources through `--resource-root`, then run preparation/analysis stages.
The separate local replay bundle can be loaded using `--input-bundle` or
CRRM_INPUT_BUNDLE; it is not a publicly distributed archive.

`input_requirements.tsv` records local replay file paths/checksums.
`configuration_manifest.tsv` maps readable parameter files; the native-code
manifest records source/portable script checksums. Current clinical source
aliases refer to 262-profile results. Historical statistical helper functions
have no retired clinical driver entry point.

## Stage order and commands

The code/scaffold entry points are in the root README; data-dependent replay requires separate inputs.
Individual native scripts run with:

```bash
crrm388 native --workspace replay_run/study --stage STAGE --script SCRIPT
```

Add `--argument VALUE` for each positional argument expected by the original
script, and `--rscript PATH` when Rscript is outside PATH. Raw stages require
the external inputs and original preprocessing/mapping records described below.
Do not execute a new protocol-declaration script over an already executed
analysis; the recorded frozen parameter files are supplied for replay.

| Stage / manuscript output | Actual scripts under `study/native/` | Inputs and replay status |
|---|---|---|
| Three-model discovery and mapping | `phase02_crrm_reconstruction/scripts/01_corrected_bulk_de.R`, `02_export_orthology.R`, `03_build_feature_mapping.R`, `05_A549_discovery_fit.R`; mapped-pool functions in `phase05_mapping_corrected_reconstruction` | Raw counts, sample manifests, frozen mapping tables and configurations required for raw DE. Evidence reconstruction is executed by `crrm388 definition`. Earlier broad-pool/anchor branches in historical native sources are not the current definition entry. The current strict three-model selector is `crrm388.py`. |
| Functional annotation and comparison | `phase10_functional_comparison_20261112/scripts/00_retrieve_public_sources.py`, `02_functional_atlas.R`, `03_Wang2026_comparison.py`, `04_IFN_reference_decomposition.R` | Source GO/Reactome/Hallmark libraries, mapping and comparator lists. Retrieval code is included; cached numerical outputs remain in the local analysis workspace. Native inference is not silently replaced by a generic enrichment tool. |
| Primary-human culture and other culture evaluation | `phase05_mapping_corrected_reconstruction/validation/GSE231666/scripts/02_frozen_human_validation.R`; `phase06_external_validation_figure1/scripts/01_manifest_and_download_counts.py` through `05_independent_analysis_checks.py` | Raw culture counts for full preprocessing. Separately held processed GSE231666 388-by-44 logCPM gives executed scoring/exact-inference regression. Historical comparison slots retained for original correction families are not a new 70-gene definition. |
| External systems / Fig. 2 | `phase18_Result2_external_response_20261112/scripts/01_prepare_external_inputs.R`, `01b_numerical_variance_guard.R`, `02_export_prepared_inputs.R`, `03_analyze_frozen_external_response.py`, `04_independent_analysis_checks.R` | Public arrays/counts, platform annotations and challenge pairing records. Current prepared results remain external prerequisites for rendering. GSE155237/GSE166161 are one external study family. |
| Clinical severity/recovery / Fig. 3, S3 | `clinical_otherwisehealthy_262_r1/scripts/01_generate_spline_designs.R`, `02_rebuild_clinical.py` | Executed locally from separately held processed expression and the fixed 262-profile inclusion list. RMA and author-corrected views have separate 262-profile references; other studies retain their native scopes. Raw preparation: `phase20_Result3_clinical_analysis_20261112/scripts/02_prepare_clinical_inputs.R`, using the current inclusion ledger. |
| Clinical incremental/conditional analysis / S4 | `clinical_otherwisehealthy_262_r1/scripts/04_rebuild_incremental.py`; upstream comparator mapping in `phase22_Result3_incremental_information_20261112/scripts/02_prepare_comparison_inputs.R` | Executed processed-input replay. Includes conditional scores, measured CyTOF composition, expression proxies, member omission and internal holdout. Reference mappings and public processed composition tables are separate input prerequisites. |
| Cellular context / Fig. 4, S5 | `Result4_analysis/scripts/01_audit_descriptor_atlas.py`, `02_freeze_atlas_analysis.py`, `03_extract_native_expression_profiles.py`, `04_score_and_test_profiles.py`, `05_cellular_interpretation.py`, `07_measurement_robustness.py`, `10_author_scope_sensitivity.py`, `15_correct_DC_time_matching.py` | Extraction needs curated h5ad, feature/availability metadata and biological-unit annotation. Large objects and row-membership metadata are external. The repository includes actual extraction/scoring code and readable frozen configurations; result tables and the full Data Descriptor resource remain external. Primary inference: infant PBMCs (19 RSV/17 healthy), airway cultures (6 donors), DC cultures (3 donors). Broad resource summaries are descriptive. |
| Cross-virus / Fig. 5, S6 | `Result5_analysis/scripts/00_audit_and_freeze.py`, `01_prepare_arrays.R`, `02_score_and_compare.py`, `03_supplementary_checks.py`, `06_extended_measurement_checks.py`, `04_independent_statistics.R` | Dallas GSE38900 GPL6884 and UK GSE42026 arrays; PBMC/DC processed cellular inputs. Uses Result-5-specific expanded-virus references. Current result/source tables remain in the local analysis workspace; raw arrays/large cells remain external. |
| Host perturbation / Fig. 6, S7 | `Result6_analysis/scripts/00_fetch_metadata.py`, `01_audit_and_freeze.py`, `02_extract_mouse_profiles.py`, `03_score_and_test.py`, `05_mouse_identity_sensitivity.py`, `04_independent_statistics.R`; current bridge `clinical_otherwisehealthy_262_r1/scripts/05_update_clinical_bridge.py` | GSE179353 author abundance matrix/culture design; GSE291797 curated mouse object/profiles. Current perturbation results and culture arrays remain in the local analysis workspace. S7E replay code uses the current raw clinical gene contrasts and all 70 whole-culture label allocations. Full raw host extraction is acquisition-dependent. |

The first three discovery R scripts take original positional roots: bulk DE
`DATA_ROOT PHASE1_MANIFEST_ROOT OUTPUT_ROOT`, orthology `OUTPUT_ROOT`, and
mapping `DATA_ROOT OUTPUT_ROOT`; the A549 script takes `DATA_ROOT OUTPUT_ROOT`.
Inspect their argument headers before raw runs. Large raw-stage prerequisites
are recorded input contracts, not files guaranteed to exist in a prepared
processed-data workspace. Archived native sources may contain earlier
exploratory branches; only the frozen current definition/configuration and
current result blocks are manuscript inputs.

## Figure entry points and source data

`crrm388 figures` invokes these seven scripts sequentially; all 13 outputs
were replayed in the recorded Linux environment. No model is refitted here.

| Figures | Stage / script | Required machine-readable inputs |
|---|---|---|
| 1, S1 | `phase16_identifier12pt_20261112 / 01_build_figures.py` | `tables/` plus primary/external score tables under `phase07_A4_7pt_figures` and `phase06_external_validation_figure1` |
| 2, S2 | `phase18_Result2_external_response_20261112 / 05_plot_figures.py` | Stage `tables/` and frozen `parameters/` |
| 3, S3 | `clinical_otherwisehealthy_262_r1/main / 07_plot_clinical_figures.py` | `main/tables/`; clinical `source_data/Fig3_SourceData.csv`, `FigS3_SourceData.csv` |
| S4 | `clinical_otherwisehealthy_262_r1/incremental / 07_plot_exploratory_figure.py` | `incremental/tables/`; clinical `source_data/FigS4_SourceData.csv` |
| 4, S5 | `Result4_analysis / 20_plot_cellular_focus_r1.py` | Stage `results/`, `tables/`, `parameters/` |
| 5, S6 | `Result5_analysis / 15_plot_crossvirus_focus_r2.py` | Stage `results/`, including `Figure5S6_source_data_crossvirus_focus_r2.csv` |
| 6, S7 | `clinical_otherwisehealthy_262_r1/host / 17_plot_component_focus_r3.py` | `host/results/`; clinical `source_data/Fig6_FigS7_SourceData.csv` |

Current clinical `source_data/CRRM_clinical_panel_results.tsv` is the complete
234-row ledger for the 24 affected clinical/host panels. Other figures have
machine-readable source-table requirements listed in `input_requirements.tsv`.

## Tables and independent checks

`crrm388 export-tables` exports current S1–S7 XLSX workbooks and S8 TSV/title,
with SHA256 values in `EXPORT_MANIFEST.json`. These are traceable frozen table
exports. Native table-building scripts are also included under their stages;
rerunning them requires their recorded prerequisite tables and R `openxlsx`.
S8 has 46 statistical design rules. Its first line is the column header;
the title is separate. Historical correction slots remain where required by
the original family, without reviving a historical definition or population.

After clinical replay and figures, optional independent checks are:

```bash
crrm388 native --workspace replay_run/study --stage clinical_otherwisehealthy_262_r1 \
  --script 10_prepare_verification.py
crrm388 native --workspace replay_run/study --stage clinical_otherwisehealthy_262_r1 \
  --script 11_independent_R_checks.R
```

The independent-check preparation step writes local verification-only inputs and relocates
obsolete copied declaration files out of active parameter folders. Statistical
settings and result tables are unchanged. The executed independent R checks
cover 21 families/131 P positions, 31 HC3 fits, 11 Welch contrasts, and 70
culture allocations. Their summaries are in `study/validation/` (summaries only).

## Public acquisition and external layout

GEO sources are accessible through
`https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=ACCESSION`.
For a specific published file, use the small acquisition helper:

```bash
python scripts/fetch_geo.py --accession GSE246622 --category matrix \
  --filename GSE246622_series_matrix.txt.gz --output-dir external_data/GSE246622
```

Use the exact deposited filename shown by GEO; an accession need not have the
illustrated single-platform matrix filename. Optional `--sha256` verifies a
recorded source hash. Native culture and functional scripts also contain
source-specific public retrieval code; `Result6_analysis/00_fetch_metadata.py`
fetches GEO/ENA/article design records. No credential is required for those
public endpoints. This helper does not invent author processed objects or
perform normalization.

Discovery raw inputs include the GSE231788 total-RNA counts/Ensembl annotation,
GSE197868 mouse count metadata, and GSE297903 author feature-count annotation.
Full clinical preparation requires deposited expression and platform mapping
plus original author-derived clinical/cytometry objects. The resource-based
cellular stages expect `DATA/expression/<source>.h5ad`, curated metadata,
feature tables, and availability masks under `--resource-root`. These are
upstream processed inputs; this repository does not rebuild the entire Data
Descriptor resource. Source-specific paths/checksums remain visible in the
normalized native parameter files. A missing external dependency must be
acquired from its source rather than substituted by zeros or guessed metadata.
