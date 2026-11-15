# CRRM388 analysis and scoring code

Code for the current **conserved RSV response module (CRRM)** manuscript.
Version 0.2.0 uses **388 unique genes, 280 up and 108 down**, and seven
fixed, overlapping annotation subsets. The current CLI is `crrm388`.

This release contains analysis/plotting code, frozen gene definitions,
readable configurations, public data/sample manifests, and validation
summaries. Current expression matrices, single-cell objects, complete
per-sample result tables, and supplementary workbooks are held separately
from this code repository.

## Installation and a runnable code example

Use the source checkout, which includes `study/frozen/` and native scripts:

```bash
python -m pip install -r requirements-study.txt
python -m pip install -e '.[study]'
crrm388 check-definition
python -m unittest discover -s tests -v

# Generate an explicitly synthetic example; no biological data are downloaded.
python scripts/make_synthetic_demo.py --output-dir replay_run/synthetic
crrm388 fit-reference --expression replay_run/synthetic/synthetic_expression.csv \
  --reference-samples replay_run/synthetic/synthetic_reference_samples.txt \
  --output replay_run/synthetic/reference.json
crrm388 score --expression replay_run/synthetic/synthetic_expression.csv \
  --reference replay_run/synthetic/reference.json --output-dir replay_run/synthetic/scores
```

The synthetic example checks software operation, not biological validation.
For real expression data, replace those paths with the study-specific matrix
and reference sample list. Bulk matrices have genes in rows and observations
in columns. `--units units.tsv` selects cellular weighting (`sample_id`, `unit`).

## Frozen definition

Discovery uses A549 total RNA at 48 h (GSE231788), mouse lung at 12 h
(GSE197868), and cotton-rat lung day-1 RSV versus day-0 PBS (GSE297903), with
equal weighting of adult and aged strata. Filtering and native DE testing
precede mapping. The common mapped/tested background is **8,607 genes**;
members have concordant direction, BH-adjusted P < 0.05, and absolute log2
fold change > 0.58 in all three models.

`study/frozen/crrm388_genes.tsv` has one row per gene;
`crrm388_subsets.tsv` records overlapping membership. The seven subsets are
IFN-alpha (47), IFN-gamma (68), TNF/NF-κB (53), TLR signaling (18),
degranulation annotation (25), RTK signaling (13), and microtubule annotation
(10). IFN-list overlap has 79 members and Outside IFN lists has 309.
The total counts each gene once. Annotation means are not added to form it.

To reconstruct membership from separately supplied mapped/tested evidence:

```bash
crrm388 definition --evidence local_inputs/common_tested_discovery.tsv.gz \
  --output-dir replay_run/definition
```

`check-definition` validates the supplied frozen lists; it does not rerun DE.
GSE155237/GSE166161 are one external challenge study family, not discovery
anchors. Unique human target projection is not uniformly strict one-to-one
orthology across species.

## Signed scoring and references

For eligible measured members of A:

`S_Aib = mean_g [s_g * (x_gib - mu_gRb) / sigma_gRb]`.

The discovery sign is fixed; expression scale and reference are study-specific.
Fit once, save gene means/SDs/eligibility, and project new observations without
refitting. Bulk uses equal profile weights and sample SD. Cellular references
use equal total biological-unit weights, equal profiles within unit, and
weighted population SD. Scores are not directly comparable across studies.

Initial primary-human scoring uses SD > 0 and raw nonzero-count support.
For a new initial-human input supply `--sd-threshold 0 --nonzero-support FILE`;
nonzero support cannot be inferred from logCPM. Later native workflows retain
their recorded numerical guards; generic CLI default is 1e-8.

Missing/invariant members are omitted, never zero-imputed. Totals require
≥70% overall and ≥50% per direction; subsets require ≥5 members and ≥50%
coverage. Failure gives missing scores and `UNESTIMABLE_COVERAGE`. Up and
reverse-down contributions use the actual total denominator and sum to the
total. Omission/common-member analyses inherit the parent reference.

## Analysis and figure workflows

```bash
# Initialize native code and readable configurations, without data matrices.
crrm388 prepare --workspace replay_run/study --data-root external_data \
  --resource-root external_resource
```

Obtain public source inputs and prepare the required stage files before running
analysis scripts. [study/WORKFLOWS.md](study/WORKFLOWS.md) maps all six Results,
Figures 1–6/S1–S7, and table-generation/export scripts to their prerequisites.
[study/input_requirements.tsv](study/input_requirements.tsv) records local
replay input/reference paths and checksums; it is not a hosted data archive.
`study/configuration_manifest.tsv` maps readable frozen configurations.

The existing local processed-input bundle can be used for author-side replay:

```bash
# This directory is supplied separately and is not downloaded from this repo.
export CRRM_INPUT_BUNDLE=../local_validation_inputs
crrm388 prepare --workspace replay_run/full
crrm388 clinical --workspace replay_run/full
crrm388 verify --workspace replay_run/full --output replay_run/verification.json
crrm388 figures --workspace replay_run/full
crrm388 export-tables --workspace replay_run/full --output-dir replay_run/tables
```

Rscript is required to regenerate age splines; `--use-saved-splines` uses the
separate saved R bases. Figure rendering requires the recorded plotting stack
and an installed Arial font. No font is distributed. Table export copies
existing local workbooks; native table builders are also included.

## Current clinical design and statistics

GSE246622 uses the author's otherwise-healthy population: Healthy 56, Mild 79,
Moderate 91, Severe 36 (262 profiles; 206 acute). A thin public inclusion list
is under `study/selection/`; matrices/covariate records remain external.
The author-corrected and RMA views each use their own 262-profile reference,
with 369 scorable members (267 up/102 down). Severe–Mild has 115 profiles.
Severity age knots are 55/177 days, boundaries 10/356; four-group gene models
use 60/186 and 10/363. Bootstrap has 1,000 grouped-profile draws; 200 Mild
splits use 39/40 profiles and exchange roles.

Internal holdout uses 114 training cases and a 261-profile gene reference
per fold. Predictor scaling is training-only; author preprocessing and the
115-profile age basis remain fixed. This is exploratory internal evaluation.
Unknown patient identities remain missing. The retired 313-profile analysis
and separate r3 restriction block are not current result blocks.

Native analyses retain mean differences, Welch/paired t, HC3, exact sign flips
and label resampling as applicable. Model/resampling P channels and original
Holm/BH families remain separate, including unavailable positions. Pointwise
intervals, bootstrap intervals, influence ranges and descriptive prediction
losses retain their distinct meanings. See the 46 rules in
[Table S8](study/tables/TableS8_statistical_design.tsv).

## Validation and historical version

Recorded local checks rebuilt 388 members from saved evidence, replayed real
44-profile primary-human scoring and the current clinical/S7E analyses, and
regenerated all 13 figures. They did not freshly rerun all raw DE/arrays or
the multi-million-cell resource. Details and tolerances are in
[LINUX_VALIDATION_REPORT.md](LINUX_VALIDATION_REPORT.md).

Default source-only tests run core/definition and legacy checks; real-data
regressions explicitly skip when separate inputs are unavailable. Setting
`CRRM_INPUT_BUNDLE` enables those regressions. CI does not claim a real-data
replay using absent inputs.

The [v0.1.0 release](https://github.com/Shiywa/conserved-host-response-framework/releases/tag/v0.1.0)
retains the historical 70-gene/five-module unsigned implementation. `chrf`,
`examples/rsv-crrm/`, original templates and `tests/test_legacy_v010.py` remain
explicit legacy material. Their Wilcoxon/Hedges' g workflow is not the current
paper's statistical analysis. See [docs/LEGACY.md](docs/LEGACY.md).

## Citation and license

Software authors remain Liu Jing and Wang Shiyou; paper authorship is separate.
Use `CITATION.cff` and cite source studies separately. No manuscript/archive
DOI is asserted. The existing MIT license/copyright are preserved. Source
data retain their original terms; see [provenance](DATA_AND_CODE_PROVENANCE.md).
