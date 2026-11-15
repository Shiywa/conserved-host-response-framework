# Changelog

## v0.2.0 - 2026-11-15

- Code-only current-study distribution: matrices/complete results stay local.

- Make the frozen 388-gene conserved RSV response module (280 up, 108 down)
  the default study entry point, with seven overlapping annotation subsets.
- Add signed reference fitting/projection, saved gene eligibility, bulk and
  equal-biological-unit cellular scaling, and current coverage gates.
- Include actual R/Python analysis and plotting workflows for Results 1–6,
  Figures 1–6/S1–S7, and Tables S1–S8 generation/export code.
- Replay the author-defined otherwise-healthy GSE246622 collection (262
  profiles; 115 strict severity profiles), incremental analyses, and S7E.
- Add local-input regression checks and release automation;
  distinguish evidence reconstruction, processed-data replay, and raw reruns.
- Preserve v0.1.0 and unsigned `chrf` behavior. Its eight compatibility tests
  and bundled example are explicitly legacy material.

## v0.1.0 - 2026-07-29

- Initial public release of the configuration-driven audit, construct, freeze,
  score, and validate workflow.
- Added a bundled RSV/CRRM acceptance example with frozen 70-gene manifests.
- Added a Linux construction replay from aggregated 571-gene evidence that
  exactly reproduces the frozen CRRM membership, directions, and modules.
- Added safety tests for duplicate genes, missing/invariant genes, study-family
  dependence, shared controls, construction leakage, paired designs, unsigned
  scoring, and deterministic freeze artifacts.
