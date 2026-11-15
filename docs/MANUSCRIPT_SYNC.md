# Author synchronization checklist

This repository update does not edit or publish manuscript text/Word files.
The manuscript should cite the actual published software version/tag/commit
and retain the distinction between historical v0.1.0 and current CRRM388.
Use the final release URL and commit recorded after publication; do not cite
an unpushed local identifier.

- Definition: conserved RSV response module; 388 unique members, 280 up and
  108 down; seven overlapping annotation subsets; IFN overlap 79/outside 309.
- GSE246622: author-defined otherwise-healthy 262 (56/79/91/36), strict 115;
  369 scored members (267/102), two separate 262-profile reference views.
- Severity and four-group spline populations/knots, 1,000 profile bootstrap
  draws, 200 Mild 39/40 splits, 114-case training and 261-profile reference
  per internal holdout must match the code. Unknown patient IDs remain missing.
- Main whole-blood Welch Holm2 P is 0.0466425310707977; separately adjusted
  label-resampling P is 0.0786. The unchanged PBMC total comparison does not
  become a replicated whole-score severity marker. Table S8 retains families.
- S7E uses raw clinical gene contrasts, 70 culture-label allocations, and
  whole-blood rho 0.4460869565217391; BH3 P 0.3857142857142857. No clinical
  treatment-efficacy claim follows from these exploratory correspondences.
- Current observed replay: Python 3.12.13, numpy 2.2.6, pandas 2.3.3,
  scipy 1.16.3, matplotlib 3.10.9, PyYAML 6.0.3; R 4.2.0. Historical package
  versions belong to their recorded stages, not a fabricated uniform runtime.
- State reproduction scope as saved-evidence construction, processed-input
  scoring/clinical replay, and figure/table generation/export. Do not describe
  this release as a new raw-data or full multi-million-cell integration run.
- Software authors remain Liu Jing and Wang Shiyou; manuscript authorship and
  eventual publication metadata are separate. No manuscript DOI is asserted.

- The current release distributes code/configurations/definitions and validation summaries. It does not distribute the processed-input bundle, full source tables, or workbooks. Real-data replay was executed locally; public CI uses core tests and an explicitly synthetic example.
