# Linux validation report: Conserved Host-Response Framework Toolkit

## Scope and boundary

Validation was performed only within this toolkit directory. Original RSV
`data/`, `results/`, figure directories, source tables, manuscript files, and
R analysis scripts were read-only inputs. No scientific result, threshold,
frozen CRRM gene, module, or manuscript file was modified.

## Phase 1: environment and structural audit

- Environment: recorded in `environment_linux.txt`.
- Toolkit structure: all generated artifacts are under `examples/rsv-crrm/outputs/`
  or `outputs/construction_replay/`; no validation command writes to the
  original project result tree.
- Codex skill validation: the locally available `quick_validate.py` returned
  `Skill is valid!`.
- Unit tests: 8/8 passed with
  `PYTHONPATH=src python3 -m unittest discover -s tests -v`.

Toolkit-only corrections made during Linux validation:

1. Added `construction.anchor_presence_column` so an adapter can faithfully
   represent an explicit legacy anchor-evidence rule without relabeling rank
   support as a primary DEG.
2. Added optional `construction.excluded_modules`, used in the RSV replay to
   enforce the original formal exclusion of `Unassigned` genes.
3. Made freeze manifests content-addressed rather than wall-clock dependent;
   repeated freezes now have identical checksum manifests.
4. Added paired-validation auditing. Incomplete pairs are reported in
   `pairing_audit.tsv`; duplicated case/control observations within a pair stop
   validation with an error.
5. Expanded tests for duplicate genes, invariant/missing genes, shared
   controls, related series, clinical-construction leakage, incomplete and
   duplicated pairs, unsigned down-gene scoring, deterministic freezing, and
   a non-RSV toy framework.

## Phase 2: bundled RSV acceptance test

The documented bundled commands were executed with
`--filter-column primary_pair --filter-value true` for validation. Outputs are
in `examples/rsv-crrm/outputs/`.

| Acceptance item | Result |
|---|---|
| Construction datasets / study families | PASS: 5 / 4 |
| GSE155237/GSE166161 related-family warning | PASS |
| GSE155237/GSE166161 shared-baseline warning | PASS |
| Frozen genes / modules | PASS: 70 / 5 |
| Frozen directions | PASS: 68 up, 2 down |
| `invert_down_genes` | PASS: false |
| GSE231666 frozen-gene coverage | PASS: 70/70 scored |
| Reference total-score agreement | PASS: max absolute error 8.88e-16 |
| Exact donor-age pairs after prespecified filter | PASS: 22 complete pairs |
| Paired total effect direction | PASS: positive paired dz = 1.37015 |

The paired audit records 22 complete pairs and no silently excluded incomplete
pairs in the accepted primary population.

## Phase 3: full Linux construction replay

Read-only source inputs were present and audited:

- `results/tables/frozen_crrm_definition/CRRM_allowed_discovery_evidence_long.csv`
- `results/tables/frozen_crrm_definition/CRRM_gene_table_full_audit.csv`
- `results/tables/refined_crrm_definition/CRRM_refined_gene_table_full.csv`
- `results/tables/refined_crrm_definition/CRRM_core_gene_table.csv`
- `results/tables/refined_crrm_definition/CRRM_refined_definition.md`

The original Figure 1 construction logic was inspected in
`scripts/38_define_frozen_crrm_signature.R`,
`scripts/39_refine_crrm_modules.R`, and the Figure 1B/1C scripts. The
toolkit-local adapter `scripts/adapt_rsv_server_evidence.py` maps the existing
server audit into the generic long schema while retaining original `detected`,
`primary_deg`, direction, and module fields. It adds the documented
`source_evidence_positive` field because formal CRRM-core eligibility requires
positive direction-aware evidence in both human fields, not a primary DEG in
both fields.

GSE155237 and GSE166161 are retained as two source fields but remain one
declared study family with shared baseline controls. They were not treated as
independent cohorts.

Replay results (`outputs/construction_replay/replay_comparison_summary.csv`):

| Metric | Value |
|---|---:|
| Original frozen core | 70 |
| Generic replay selected genes | 70 |
| Exact shared membership | 70 |
| Missing from replay | 0 |
| Additional in replay | 0 |
| Direction mismatches | 0 |
| Module mismatches | 0 |

The generic logic reproduced the original 70-gene membership without changing
scientific selection thresholds. The run uses the original pre-freeze module
assignments and excludes only the original `Unassigned` module category.

Two independent `freeze` executions from the replay selected table produced
identical `SHA256SUMS` files.

## Phase 4: safety and robustness coverage

The expanded unit tests demonstrate toolkit-level safeguards for:

- duplicate gene symbols;
- invariant and missing frozen genes;
- explicit low coverage reporting;
- shared controls and related construction series;
- clinical/outcome data incorrectly assigned to construction;
- complete, incomplete, and duplicated paired designs;
- unsigned equal-weight scoring with `invert_down_genes: false`;
- deterministic frozen-manifest checksums; and
- toy non-RSV software behavior.

These are software tests, not cross-disease biological validation.

## Recommendation

**Ready as an internal reproducibility tool.**

All stated Linux acceptance and replay checks passed. Public archival release
still requires repository-level review and a second real-disease biological
demonstration before any cross-disease scientific generalizability claim. The
toolkit must not be described as clinically predictive, causal, or universally
validated.
