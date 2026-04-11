# Quickstart

## Fastest path

Run the model with the bundled synthetic stress-test dataset:

```bash
python rlto_model.py
```

This executes `example_usage()` and will:

* evaluate and optimize all demo genes,
* compute a global tumor-level optimum,
* optimize inhibition against `TOXIC_DRIVER`,
* generate diagnostic plots under `diagnostic_plots/`,
* print a final JSON summary.

## Minimal import example

```python
from rlto_model import CancerRLTOModel, Microenvironment

env = Microenvironment(
    hypoxia=0.6,
    nutrient_limitation=0.4,
    immune_pressure=0.0,
    drug_pressure=0.0,
    oxidative_stress=0.0,
)

df = CancerRLTOModel.demo_dataset()
model = CancerRLTOModel.from_dataframe(df, microenvironment=env)

results = model.evaluate()
print(results[["gene", "net_fitness", "robustness_probability"]])
```

## Optimize all genes

```python
opt_df = model.optimize_all()
print(opt_df[["gene", "baseline_abundance", "optimal_abundance", "optimal_net_fitness"]])
```

## Optimize the global tumor state

```python
opt_vector, best_global_fitness = model.optimize_global()
print(best_global_fitness)
```

## Simulate and optimize an intervention

```python
intervention_df = model.simulate_intervention("TOXIC_DRIVER", inhibition=0.4)
best_inhibition, residual_fitness = model.optimize_inhibition("TOXIC_DRIVER")

print(best_inhibition)
print(residual_fitness)
```

## Generate diagnostics

```python
from rlto_model import ModelDiagnostics

diag = ModelDiagnostics(model)
fig = diag.plot_component_decomposition("MYC")
fig.savefig("MYC_component.png")
```

## Interpreting the first results

Focus on these fields first:

* `net_fitness`: the composite objective value per gene,
* `robustness_probability`: analytical survival probability above threshold,
* `overabundance`: mean abundance divided by threshold abundance,
* `optimal_abundance`: abundance maximizing gene-level net fitness.

See [Methods and assumptions](methods.md) for the model interpretation.