# Usage

## Entry points

The repository currently exposes two practical usage modes:

1. direct script execution via `python rlto_model.py`,
2. importing objects and functions from `rlto_model.py`.

There is no dedicated CLI parser in the current codebase.

## Script mode

```bash
python rlto_model.py
```

This triggers `example_usage()` and performs a full demonstration workflow on the built-in synthetic dataset.

### What the script does

1. constructs a moderately stressful `Microenvironment`,
2. loads the synthetic demo dataframe,
3. creates a `CancerRLTOModel`,
4. runs gene-level optimization across all genes,
5. performs tumor-level global optimization,
6. optimizes inhibition for `TOXIC_DRIVER`,
7. writes diagnostic plots into `diagnostic_plots/`,
8. prints a JSON summary.

## Import mode

### Build from the demo dataset

```python
from rlto_model import CancerRLTOModel, Microenvironment

env = Microenvironment(hypoxia=0.6, nutrient_limitation=0.4)
model = CancerRLTOModel.from_dataframe(
    CancerRLTOModel.demo_dataset(),
    microenvironment=env,
)
```

### Build from your own dataframe

```python
import pandas as pd
from rlto_model import CancerRLTOModel, Microenvironment

df = pd.read_csv("genes.csv")
env = Microenvironment(hypoxia=0.2, drug_pressure=0.5)

model = CancerRLTOModel.from_dataframe(df, microenvironment=env)
```

## Common operations

### Evaluate the current state

```python
df = model.evaluate()
```

### Summarize by pathway

```python
pathways = model.pathway_summary()
```

### Compute a scalar tumor score

```python
score = model.tumor_fitness_score()
```

### Optimize one gene

```python
gene = model.genes[0]
best = model.optimize_gene_abundance(gene)
```

### Optimize all genes independently

```python
opt_df = model.optimize_all()
```

### Optimize the whole tumor state jointly

```python
x_opt, f_opt = model.optimize_global(global_resource_weight=0.5)
```

### Simulate inhibition

```python
post = model.simulate_intervention("MYC", inhibition=0.3)
```

### Find the best inhibition level

```python
u_opt, f_residual = model.optimize_inhibition("MYC")
```

## Recommended analysis pattern

```python
results = model.evaluate()
opt_gene = model.optimize_all()
pathways = model.pathway_summary()
score = model.tumor_fitness_score()
```

Then move to diagnostics only when you need to inspect model shape or optimization behavior.

## Reproducibility notes

* The main analytical fitness calculation is deterministic.
* The model still initializes an RNG and exposes `sample_abundance()`.
* The diagnostics API contains labels such as “noise stability”, but the current analytical implementation does not rely on Monte Carlo sampling for `net_gene_fitness()`.

!!! note
    The `n_samples` parameter is preserved in several method signatures for compatibility, but the analytical path means it is largely not used by the core fitness calculation.
