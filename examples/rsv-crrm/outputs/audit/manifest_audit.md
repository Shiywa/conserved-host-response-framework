# Dataset-role and dependence audit

- Construction datasets: 5
- Construction study families: 4
- Post-freeze datasets: 15

| Check | Status | Detail |
|---|---|---|
| Unique dataset identifiers | PASS | All dataset_id values are unique. |
| Construction datasets present | PASS | 5 construction datasets. |
| Post-freeze datasets present | PASS | 15 non-construction datasets. |
| Outcome and cellular data excluded from construction | PASS | No downstream outcome/cellular roles enter construction. |
| Construction study-family independence | WARN | RSV_human_challenge_2020: GSE155237, GSE166161 |
| Shared construction controls | WARN | RSV_challenge_shared_baseline: GSE155237, GSE166161 |
| Construction/post-freeze dataset separation | PASS | No dataset identifier appears in both roles. |
| Dataset source provenance | PASS | Every dataset has a source reference or repository identifier. |

Warnings disclose dependence or incomplete provenance; they do not automatically invalidate a framework.
Failures must be resolved before construction or post-freeze validation is interpreted.
