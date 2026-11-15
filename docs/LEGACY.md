# Historical v0.1.0 material

The historical release is
[v0.1.0](https://github.com/Shiywa/conserved-host-response-framework/releases/tag/v0.1.0),
original commit `a1025aa8c671f6f14004ba164a6c7cc1d3797d3f`.
It used 70 genes, five modules, and an explicitly unsigned example specification.

The `chrf` CLI, generic `framework_config_template.yaml`, old dataset/contrast
templates, `examples/rsv-crrm`, `outputs/construction_replay`, the original
toolkit specification, and bundled Codex skill retain this context. The
current CLI is `crrm388`; its reference and statistical modules are separate.
Legacy score/validate semantics were not silently changed. Legacy Wilcoxon,
Hedges' g, and 571-to-70 evidence replay do not validate the current manuscript.

The eight historical tests are in `tests/test_legacy_v010.py`. Current
scientific migration tests are in `tests/test_crrm388.py`. Run either using
unittest discovery with the corresponding filename pattern, or run both.

Historical reports remain available at the v0.1.0 tag. Current root reports
describe v0.2.0. No old release, tag, or copyright is removed or redirected.
