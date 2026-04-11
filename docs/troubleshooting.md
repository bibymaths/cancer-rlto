# Troubleshooting

## `ImportError: No module named scipy` or similar

Install runtime dependencies first:

```bash
uv sync
```

## Plots are not being written

Check that the script can create the `diagnostic_plots/` directory in the current working directory.

Also verify that your environment has a working Matplotlib backend for non-interactive figure saving.

## Optimization returns unexpected values

Review the inputs first:

* very high `toxicity_weight` will strongly suppress abundance,
* `essentiality_weight = 0` removes the robustness benefit,
* very low `clone_fraction` can collapse baseline abundance,
* high `burstiness` may push abundance optima upward.

## A gene raises `ValueError`

`GeneConfig.validate()` rejects invalid parameter ranges. Check:

* positive transcription rate,
* positive translation efficiency,
* positive copy number,
* clone fraction within `[0, 1]`,
* non-negative weights,
* positive burstiness.