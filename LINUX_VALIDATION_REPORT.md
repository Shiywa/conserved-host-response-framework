# CRRM388 Linux replay record — 2026-11-15

## Executed scope

The current release was checked against independently saved scientific outputs,
not values filled into a fresh expected table. Discovery reconstruction starts
from the full common mapped/tested DE evidence; primary-human and whole-blood
replay starts from processed expression. All 13 figures were regenerated from
current result tables. Raw DE and large-cell extraction were not freshly rerun.

Commands (from the source checkout; observed environment in environment_linux.txt):

```bash
export CRRM_INPUT_BUNDLE=../local_validation_inputs
python -m unittest discover -s tests -v
crrm388 definition --output-dir replay_run/definition
crrm388 prepare --workspace replay_run/study
crrm388 clinical --workspace replay_run/study
crrm388 figures --workspace replay_run/study
crrm388 verify --workspace replay_run/study --output replay_run/verification.json
crrm388 export-tables --workspace replay_run/study --output-dir replay_run/tables
crrm388 native --workspace replay_run/study --stage clinical_otherwisehealthy_262_r1 --script 10_prepare_verification.py
crrm388 native --workspace replay_run/study --stage clinical_otherwisehealthy_262_r1 --script 11_independent_R_checks.R
```

The full local 28-test suite contains 20 current checks and eight explicit v0.1.0
compatibility checks. Current checks cover discovery signs/uniqueness, down
reversal, overlapping subsets, immutable fitted references, actual missing and
invariant handling, coverage gates, unit weighting/SD, weighted contributions,
current statistical kernels/family sizes, and independently saved clinical
results. The native processed clinical pipeline is executed when CRRM_INPUT_BUNDLE supplies the separate local inputs. The code-only default runs 21 tests and skips one regression class containing seven data-dependent tests. No matrices or complete result tables are distributed in this release.

## Recorded numerical checks

| Check | Result |
|---|---|
| Common discovery background | 8,607 |
| Frozen definition | 388 unique; 280 up / 108 down; exact match |
| Primary-human profiles | 44; max score error 1.2212453270876722e-15 |
| Primary-human exact inference | max error 6.661338147750939e-16 |
| Clinical numeric tables compared | 54; max numeric error 0.0 |
| Generic signed API vs native whole-blood score | max error 1.6653345369377348e-16 |
| Whole-blood population/reference | 262; 369 measured members, 267 up / 102 down |
| Whole-blood total Severe–Mild | 0.095171119063625 |
| Welch raw / Holm2 P | 0.0233212655353988 / 0.0466425310707977 |
| Monte Carlo raw / separate Holm2 P | 0.0393 / 0.0786 |
| Whole-blood degranulation effect | 0.420539 (full precision in source table) |
| S7E whole-blood rho | 0.4460869565217391 |
| S7E permutation / BH3 P | 18/70 = 0.2571428571428571 / 0.3857142857142857 |
| Recovery | GSE69606 original 15 pairs unchanged |
| Internal evaluation | 115 profiles; training 114; gene reference 261 per fold |
| Figures | All 13 regenerated PNGs byte-identical to current author figure files |

The 54-table comparison includes retained other-study tables as well as
recomputed clinical outputs; it is not a claim that every row was refitted.
The declared Python comparison tolerance is 1e-10. The independent R checks
matched 131 adjustment positions (max error 1.11e-16), 31 HC3 fits (max error
7.60e-12), 11 Welch contrasts (max error 1.11e-16), and the 70 S7E allocations.
R HC3/Welch verification uses 1e-8; the larger HC3 tolerance accommodates
independent linear-algebra implementations and is not used to change results.

## Reproduction tiers

1. **Saved evidence reconstruction:** three-model mapped/tested effects and FDR
   -> frozen members and signs; no fresh raw DE claim.
2. **Processed-expression replay:** real primary-human scores/exact inference,
   current 262-profile clinical models, 1,000 grouped profile-bootstrap draws,
   200 swapped Mild 39/40 splits, internal holdout, S7E and separate P families.
   The full local run regenerated both spline bases with R. Locally enabled regression tests use saved
   R bases and do not claim to rerun R. Public CI runs source-only tests and a synthetic CLI example.
3. **Rendering and export:** all 13 figure scripts; S1–S8 traceable table export.
   PNG byte identity is environment-specific. Figure output generation is not
   raw expression reprocessing.
4. **Raw/acquisition-dependent stages:** actual source code is included with
   original scientific parameters; raw counts/arrays and large h5ad resources
   were not all rerun or bundled.

Validation summaries/hashes are in `study/validation/`; the underlying data and complete comparison tables remain local. The superseded 313-profile
clinical results and r3 restriction block are not current regression targets.
Historical eight-test/70-gene acceptance records remain at the v0.1.0 tag.
