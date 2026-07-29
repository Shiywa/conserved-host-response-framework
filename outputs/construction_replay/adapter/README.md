# RSV server-evidence adapter

This directory is generated inside the toolkit. It reads the original project tables as read-only. The adapter does not redefine a gene, module, direction, or threshold.

`source_evidence_positive` is an explicit adapter field, not a replacement for `primary_deg`. It faithfully represents the original formal requirement that both related human fields have positive direction-aware evidence. GSE155237 and GSE166161 remain related fields in one declared study family and must not be described as independent cohorts.
