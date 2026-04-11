# Inputs

## Overview

The main input abstraction is `GeneConfig`, either provided directly or constructed from a pandas dataframe using `CancerRLTOModel.from_dataframe()`.

A secondary input is `Microenvironment`, which encodes environmental stress conditions.

## Gene-level inputs

### Required field

Only one dataframe column is strictly required beyond identifiers:

- `gene`
- `transcription_rate`

### Optional fields

The code supports the following optional columns:

| Column | Type | Default | Meaning |
|---|---:|---:|---|
| `translation_efficiency` | float | `1.0` | Translation gain factor |
| `mrna_half_life` | float | `6.0` | mRNA half-life |
| `protein_half_life` | float | `12.0` | Protein half-life |
| `copy_number` | float | `2.0` | Effective gene dosage |
| `clone_fraction` | float | `1.0` | Fraction of clones carrying the gene state |
| `essentiality_weight` | float | `1.0` | Weight on robustness benefit |
| `threshold_abundance` | float or null | inferred | Explicit critical abundance threshold |
| `baseline_abundance` | float or null | inferred | Explicit baseline abundance override |
| `burden_weight` | float | `0.02` | Linear burden weight |
| `toxicity_weight` | float | `0.0` | Superlinear toxicity weight |
| `regulation_strength` | float | `0.0` | Stabilizing regulation term |
| `burstiness` | float | `1.0` | Noise control used in gamma parameterization |
| `stress_sensitivity` | float | `0.2` | Environmental stress sensitivity |
| `oncogenic_boost` | float | `0.0` | Abundance boost from oncogenic pressure |
| `pathway` | str | `generic` | Pathway/grouping label |

## Microenvironment inputs

`Microenvironment` defines five bounded stress axes:

- `hypoxia`
- `nutrient_limitation`
- `immune_pressure`
- `drug_pressure`
- `oxidative_stress`

Each value is clipped to `[0, 1]`.

## Validation rules

`GeneConfig.validate()` enforces:

- `transcription_rate > 0`
- `translation_efficiency > 0`
- `copy_number > 0`
- `clone_fraction` in `[0, 1]`
- `essentiality_weight >= 0`
- `burden_weight >= 0`
- `toxicity_weight >= 0`
- `burstiness > 0`

Invalid values raise `ValueError`.

## Example input table

```csv
gene,transcription_rate,essentiality_weight,burden_weight,toxicity_weight,pathway
MYC,8.0,0.8,0.015,0.020,proliferation
KRAS,5.5,0.9,0.012,0.015,MAPK
PASSENGER,2.0,0.0,0.01,0.01,none
```

## Demo dataset

The module ships with `CancerRLTOModel.demo_dataset()`, a synthetic stress-test dataset containing:

* typical anchor genes,
* zero-essentiality and zero-penalty edge cases,
* noise extremes,
* burden and toxicity extremes,
* microenvironment-sensitive genes.

!!! warning
    If your real project expects omics-derived inputs, sample metadata, or file-based ingestion beyond pandas dataframe construction, that behavior is not present in the supplied module and should be documented only after implementation.
