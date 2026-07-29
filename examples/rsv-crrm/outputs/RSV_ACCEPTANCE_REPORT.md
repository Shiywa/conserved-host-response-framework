# RSV/CRRM Linux acceptance report

| Check | Status | Detail |
|---|---|---|
| Related human study family detected | PASS | RSV_human_challenge_2020: GSE155237, GSE166161 |
| Shared baseline detected | PASS | RSV_challenge_shared_baseline: GSE155237, GSE166161 |
| Frozen gene count | PASS | 70 |
| Frozen module count | PASS | 5 |
| Frozen direction counts | PASS | {"down": 2, "up": 68} |
| GSE231666 frozen-gene scoring coverage | PASS | 70/70 |
| GSE231666 total-score numerical reproduction | PASS | max_abs_error=8.88e-16 |
| GSE231666 exact-pair count | PASS | 22 pairs |
| GSE231666 paired total effect direction | PASS | effect=1.37015 |

This acceptance test validates role auditing, freezing, dataset-relative scoring, and paired post-freeze reuse.
It does not validate a complete 571-to-70 construction replay unless the original server-side pre-freeze evidence is adapted and tested separately.
