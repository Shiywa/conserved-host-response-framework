# Changelog

## v0.1.0 - 2026-07-29

- Initial public release of the configuration-driven audit, construct, freeze,
  score, and validate workflow.
- Added a bundled RSV/CRRM acceptance example with frozen 70-gene manifests.
- Added a Linux construction replay from aggregated 571-gene evidence that
  exactly reproduces the frozen CRRM membership, directions, and modules.
- Added safety tests for duplicate genes, missing/invariant genes, study-family
  dependence, shared controls, construction leakage, paired designs, unsigned
  scoring, and deterministic freeze artifacts.
