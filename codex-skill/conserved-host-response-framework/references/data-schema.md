# Data schema

## Dataset manifest

Required columns:

- `dataset_id`: stable resource identifier.
- `study_family`: publication, challenge, trial, or experiment family.
- `role`: `construction` or a named post-freeze role.
- `biological_role`: infection, severity, recovery, cellular localization, cross-virus context, or another declared use.
- `shared_control_group`: identifier shared by datasets reusing controls; blank otherwise.
- `species`, `tissue`: biological context.
- `source_reference`: DOI, accession, or repository record.

Recommended columns:

- `endpoint`, `replicate_unit`, `paired_design`, `platform`, `orthology_rule`, `freeze_relation`.

## Construction evidence

One row per gene/dataset evidence field:

- `gene`, `dataset_id`, `system`, `direction`;
- `detected`, `primary_deg`, `supported`;
- `module`.

Use `up`, `down`, or `neutral` for direction and explicit booleans for evidence flags.

## Expression matrix

Use one row per gene and one column per biological sample. The first column must be the gene identifier. Duplicate symbols are collapsed by mean. Scoring uses within-dataset gene-wise z scores and reports unmeasured or invariant frozen genes.

