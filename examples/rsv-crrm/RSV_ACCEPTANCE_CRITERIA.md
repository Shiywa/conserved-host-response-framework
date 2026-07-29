# RSV/CRRM acceptance criteria

- Dataset audit reports five construction datasets and four construction study families.
- GSE155237 and GSE166161 trigger study-family and shared-baseline warnings.
- Severity, recovery, single-cell, and cross-virus resources do not enter construction.
- Freeze output contains 70 unique genes, five modules, 68 up directions, and two down directions.
- Sample scoring detects and scores all 70 genes in the GSE231666 precomputed-z input.
- Recomputed total scores match the manuscript source table within numerical tolerance.
- Paired validation retains 22 donor-age pairs and a positive total CRRM effect.
- Construction command is tested on a toy dataset, but the complete 571-to-70 RSV reconstruction is not claimed without the server-side pre-freeze evidence tables.

