# Contributing

Contributions are welcome through issues and pull requests.

## Principles

1. Do not alter a frozen framework because of downstream validation results.
2. Declare dataset roles, related study families, and shared controls before
   construction.
3. Preserve negative, heterogeneous, and low-coverage results in outputs.
4. Do not add restricted data, credentials, participant-identifiable metadata,
   manuscript files, or unpublished figures to the repository.
5. Add focused tests for every behavior change.

## Development checks

```bash
python -m pip install -e .
PYTHONPATH=src python -m unittest discover -s tests -v
```

Changes to `examples/rsv-crrm/` must preserve the frozen 70-gene definition,
the declared GSE155237/GSE166161 dependence, and the unsigned equal-weight
scoring rule unless a new, explicitly versioned example is added.
