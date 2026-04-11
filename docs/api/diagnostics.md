# Diagnostics API

## `ModelDiagnostics`

::: rlto_model.ModelDiagnostics

## `example_usage`

::: rlto_model.example_usage

## Diagnostics guide

The plotting methods are intended to answer different questions:

| Method                           | Intended diagnostic question                                                      |
|----------------------------------|-----------------------------------------------------------------------------------|
| `plot_component_decomposition()` | do robustness, burden, toxicity, and net fitness produce a well-shaped objective? |
| `plot_noise_stability()`         | how sensitive is the local objective region around the optimum?                   |
| `plot_convergence_check()`       | does optimization converge from varied starting points?                           |
| `plot_bivariate_landscape()`     | does the joint landscape show compensation ridges or a clear optimum?             |
| `plot_sensitivity_profile()`     | does the optimum move smoothly under parameter perturbation?                      |

## Saving plots

```python
from rlto_model import ModelDiagnostics

diag = ModelDiagnostics(model)
fig = diag.plot_bivariate_landscape("MYC", "KRAS")
fig.savefig("bivariate_MYC_KRAS.png")
```

## Operational note

The example script closes figures after saving them. That pattern is worth preserving in batch runs to avoid figure
accumulation.
