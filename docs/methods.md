# Methods and assumptions

## Model philosophy

The implementation encodes a robustness–load trade-off for cancer systems biology.

At a high level:

- higher abundance can increase the probability of remaining above a critical threshold,
- but abundance also incurs burden and toxicity penalties,
- and the tumor-level optimizer applies a global resource constraint.

## Deterministic analytical core

The current code uses an analytical gamma survival function to compute the probability that abundance exceeds the
critical threshold.

That matters because it makes the main fitness calculation deterministic.

```python
robustness_prob = gammaincc(shape, threshold / scale)
```

## Abundance estimation

If `baseline_abundance` is not explicitly supplied, it is inferred from:

* transcription rate,
* translation efficiency,
* copy number,
* oncogenic boost,
* clone fraction,
* regulation strength,
* stress penalty,
* mRNA and protein decay rates.

Conceptually:

```text
abundance ∝ synthesis × regulation_gain × stress_penalty / (mrna_decay × protein_decay)
```

## Threshold abundance

If `threshold_abundance` is not explicitly provided, it is inferred from:

* square-root dependence on baseline abundance,
* stress sensitivity and environmental stress,
* essentiality,
* regulation-based discounting.

## Noise model

Protein abundance is parameterized as a gamma distribution.

The code derives:

* `shape = 1 / cv²`
* `scale = mean / shape`

with `cv²` influenced by:

* burstiness,
* transcription rate,
* copy number,
* regulation strength,
* microenvironment stress.

## Fitness components

Per-gene net fitness is modeled as:

```text
net = robustness_term - burden_cost - toxicity_cost
```

### Robustness term

```text
robustness_term = robustness_scale × expected_growth × essentiality_weight
```

where `expected_growth` is approximated from the analytical survival probability.

### Burden cost

Linear in abundance after scaling:

```text
burden_cost ∝ burden_weight × (abundance / 1000)
```

### Toxicity cost

Superlinear in abundance:

```text
toxicity_cost ∝ toxicity_weight × (abundance / 1000)^1.25 / (1 + regulation_strength)
```

## Optimization layers

### Layer 1: gene-level optimization

`optimize_gene_abundance()` uses bounded scalar minimization to maximize per-gene net fitness.

### Layer 2: global tumor optimization

`optimize_global()` uses SLSQP over all genes jointly and subtracts a quadratic penalty on total abundance normalized by
baseline total abundance.

### Layer 3: intervention optimization

`simulate_intervention()` scales the target gene baseline by `(1 - inhibition)`.

`optimize_inhibition()` searches `u ∈ [0, 1]` to minimize mean tumor fitness after intervention.

## Important implementation note

Several methods still accept an `n_samples` parameter, but `net_gene_fitness()` is now analytical and no longer depends
on Monte Carlo sampling for the main calculation.

## Scope and limits

This implementation should be read as a practical engineering model, not a full formal recovery of any paper’s
supplementary derivation.
