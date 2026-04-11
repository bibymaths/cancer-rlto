# Outputs

## Overview

The module returns either structured Python objects, pandas dataframes, scalar optimization results, or generated plots.

## Structured result object

`net_gene_fitness()` returns a `ModelResult` dataclass with:

- `gene`
- `baseline_abundance`
- `threshold_abundance`
- `overabundance`
- `abundance_mean`
- `abundance_std`
- `robustness_probability`
- `burden_cost`
- `toxicity_cost`
- `net_fitness`
- `pathway`

## Tabular outputs

### `evaluate()`

Returns a dataframe sorted by descending `net_fitness` and then `overabundance`, with an added `rank` column.

### `optimize_all()`

Returns a dataframe with:

- `gene`
- `baseline_abundance`
- `optimal_abundance`
- `critical_threshold`
- `optimal_overabundance`
- `optimal_net_fitness`

### `pathway_summary()`

Returns grouped statistics by pathway:

- number of genes,
- mean and median overabundance,
- mean net fitness,
- mean robustness.

### `simulate_intervention()`

Returns a dataframe of post-inhibition per-gene fitness results.

## Scalar outputs

### `tumor_fitness_score()`

Returns a scalar tumor score derived from the clipped mean gene fitness.

### `optimize_global()`

Returns a tuple:

```python
(optimal_abundance_vector: np.ndarray, best_global_fitness: float)
````

### `optimize_inhibition()`

Returns a tuple:

```python
(best_inhibition: float, residual_fitness: float)
```

## Script-generated files

When running:

```bash
python rlto_model.py
```

the code creates:

```text
diagnostic_plots/
```

and writes figures such as:

```text
diagnostic_plots/MYC_component.png
diagnostic_plots/MYC_noise.png
diagnostic_plots/MYC_convergence.png
diagnostic_plots/MYC_sensitivity.png
diagnostic_plots/bivariate_MYC_KRAS.png
diagnostic_plots/bivariate_CHAOS_HOUSEKEEPER.png
diagnostic_plots/bivariate_TOXIC_STRUCTURAL.png
...
...
...
```

## Example result interpretation

| Field                         | Interpretation                                               |
|-------------------------------|--------------------------------------------------------------|
| `overabundance > 1`           | mean abundance exceeds critical threshold                    |
| high `robustness_probability` | abundance distribution remains above threshold more reliably |
| high `burden_cost`            | production cost is materially penalizing                     |
| high `toxicity_cost`          | superlinear overexpression penalty is material               |
| low or negative `net_fitness` | the gene state is unfavorable under current settings         |

## JSON-style terminal summary

The script entrypoint ends by printing a JSON payload similar to:

```json
{
  "Best_global_fitness": 5.8508,
  "Optimal_Inhibition_Target": 0.000006,
  "All_Genes_Optimization": [
    {
      "gene": "MYC",
      "baseline_abundance": 2301.83,
      "optimal_abundance": 70.10,
      "optimal_net_fitness": 7.92
    }
  ]
}
```

!!! note
Exact values depend on the model configuration and dataset. The numeric example above reflects the current README
narrative and should be treated as an expected pattern rather than a guaranteed regression fixture.
