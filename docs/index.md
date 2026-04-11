# Cancer RLTO Model

A cancer-adapted RLTO-inspired Python model for deterministic analysis of robustness, burden, toxicity, optimization, and stress-test diagnostics.

## What this project does

The repository provides a single Python module, `rlto_model.py`, that models tumor-cell fitness as a trade-off between:

- robustness against abundance fluctuations,
- linear burden from producing protein,
- superlinear toxicity from overabundance,
- and a global resource constraint at the tumor level.

The implementation supports:

- per-gene fitness evaluation,
- optimization of gene abundance,
- global optimization across all genes,
- simulation and optimization of gene inhibition,
- pathway-level summaries,
- synthetic stress-test data,
- and diagnostic plots for landscape shape, convergence, sensitivity, and bivariate interactions.

## Project status

This documentation is generated from the current repository state and is intentionally aligned to the code rather than to narrative claims in older prose.

!!! note "Current repository shape"
    The project currently looks like a script-oriented module rather than a packaged library. The documentation therefore emphasizes direct script execution and module import patterns.

## Typical workflow

1. Build a dataframe of gene parameters or use the built-in demo dataset.
2. Construct a `Microenvironment`.
3. Build a `CancerRLTOModel` from the dataframe.
4. Evaluate fitness or run optimization routines.
5. Generate diagnostic plots if needed.
6. Interpret outputs in terms of robustness, overabundance, burden, toxicity, and intervention response.

## Quick example

```python
from rlto_model import CancerRLTOModel, Microenvironment

env = Microenvironment(hypoxia=0.6, nutrient_limitation=0.4)
df = CancerRLTOModel.demo_dataset()
model = CancerRLTOModel.from_dataframe(df, microenvironment=env)

summary = model.evaluate()
optimized = model.optimize_all()
score = model.tumor_fitness_score()

print(summary.head())
print(optimized[["gene", "optimal_abundance", "optimal_net_fitness"]])
print(score)
```

## Output conventions

Running the module as a script:

```bash
python rlto_model.py
```

will:

* print optimization summaries to stdout,
* optimize the built-in synthetic “Rogue’s Gallery” stress-test dataset,
* and write diagnostic figures into `diagnostic_plots/`.

See [Outputs](outputs.md) for the expected file structure. 

## Feature summary

=== "Core model"

    - `Microenvironment` to represent environmental stress factors.
    - `GeneConfig` for per-gene parameterization.
    - `CancerRLTOModel` for deterministic analytical evaluation and optimization.
    - `ModelResult` for structured result output.

=== "Analysis workflow"

    - baseline abundance estimation,
    - threshold abundance calculation,
    - gamma-distributed abundance modeling,
    - robustness probability from the analytical gamma survival function,
    - burden and toxicity costs,
    - net gene fitness,
    - tumor-level optimization with a quadratic resource penalty.

=== "Diagnostics"

    - component decomposition plots,
    - local noise-stability inspection,
    - convergence checks across start points,
    - bivariate landscape plots,
    - sensitivity profiling.