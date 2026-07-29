# Unresolved limitations

1. The toolkit replay verifies that the original server-side, already
   aggregated evidence tables can reproduce the frozen 70-gene CRRM. It does
   not rerun raw-expression differential expression or rank integration from
   the 13 source contrasts.
2. GSE155237 and GSE166161 share a study family and baseline controls. The
   toolkit preserves and reports this dependence; they must not be described
   as two independent human cohorts.
3. The adapter relies on the original project's audited 0/1/2 evidence fields
   and pre-freeze module assignments. It does not independently validate the
   biological annotation or ortholog mapping underlying those fields.
4. The GSE231666 acceptance example demonstrates numerical score reproduction
   and a paired post-freeze association, not diagnosis, prediction,
   calibration, causal validity, or clinical utility.
5. The toy non-RSV unit test demonstrates disease-independent software
   behavior only. A second real disease framework would be required for a
   cross-disease scientific generalizability claim.
6. The included LICENSE explicitly does not grant public redistribution or
   public-use rights. Author-approved licensing, repository metadata, and
   archival review remain required before public release.
