# API overview

## Public reference surface

The documentation treats the following objects as the effective public API:

- `Microenvironment`
- `GeneConfig`
- `ModelResult`
- `CancerRLTOModel`
- `ModelDiagnostics`
- `example_usage`

## Reference sections

- [Core model API](model.md)
- [Diagnostics API](diagnostics.md)

## Design notes

This is a curated API reference, not a raw symbol dump.

The focus is on the parts that are useful for:

- constructing models,
- evaluating and optimizing fitness,
- simulating perturbations,
- generating diagnostics.

Low-value internals and helper details are intentionally not elevated beyond what `mkdocstrings` already renders for the public objects.

## Import convention

```python
from rlto_model import (
    Microenvironment,
    GeneConfig,
    ModelResult,
    CancerRLTOModel,
    ModelDiagnostics,
    example_usage,
)
```