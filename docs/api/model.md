# Core model API

## Data containers

### `Microenvironment`

::: rlto_model.Microenvironment

### `GeneConfig`

::: rlto_model.GeneConfig

### `ModelResult`

::: rlto_model.ModelResult

## Main model

### `CancerRLTOModel`

::: rlto_model.CancerRLTOModel

## Practical usage notes

### Recommended constructor path

For most workflows, prefer dataframe ingestion:

```python
from rlto_model import CancerRLTOModel, Microenvironment

env = Microenvironment(hypoxia=0.3, drug_pressure=0.2)
df = CancerRLTOModel.demo_dataset()
model = CancerRLTOModel.from_dataframe(df, microenvironment=env)
```

### Methods you will use most often

| Method                      | Purpose                                    |
|-----------------------------|--------------------------------------------|
| `evaluate()`                | rank the current gene set                  |
| `optimize_gene_abundance()` | optimize one gene                          |
| `optimize_all()`            | optimize all genes independently           |
| `optimize_global()`         | jointly optimize the tumor state           |
| `simulate_intervention()`   | apply an inhibition fraction to one target |
| `optimize_inhibition()`     | search for the best inhibition level       |
| `pathway_summary()`         | group outputs by pathway                   |
| `tumor_fitness_score()`     | derive a scalar summary score              |

### Behavior caveat

Several methods retain `n_samples` arguments for compatibility, but the main fitness calculation is analytical.
