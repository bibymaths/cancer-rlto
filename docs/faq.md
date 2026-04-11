# FAQ

## Is this a packaged Python library?

Not yet. The supplied repository state looks like a single-module project rather than a fully packaged distribution.

## Is the model stochastic?

The abundance distribution is stochastic in concept, but the main fitness calculation is currently analytical and deterministic.

## Why do some methods still accept `n_samples`?

Backwards compatibility. The signatures were preserved even though `net_gene_fitness()` no longer relies on Monte Carlo estimation.

## Is there a CLI?

No dedicated CLI parser exists in the provided code. Current usage is either script execution or importing the module in Python.

## Can I use my own dataset?

Yes, if you can convert it into a dataframe with the required columns, especially `gene` and `transcription_rate`.

## Where do the plots go?

The example workflow writes figures into `diagnostic_plots/`.

## Why are the docs not embedding the figures from the README?

Because no image files were supplied. The docs describe expected outputs and where they should appear, but they do not pretend those assets exist in the repository.

## Does the project include tests or packaging metadata?

Not in the provided files.

## Is the model tied exactly to a published RLTO paper?

The code frames itself as an RLTO-inspired engineering implementation adapted to cancer systems biology. It should not be documented as a strict formal reconstruction unless supporting derivations are added.
