# Plot select Lee's L curves

Plot Lee's L curves from select gene pairs, coloring by biological
conditions

## Usage

``` r
plotLeeSelect(df, lmm_res, pairs_use, title = NULL)
```

## Arguments

- df:

  Data frame with Lee's L results from all samples, with columns pair,
  side, lee, sample, and group (for biological conditions).

- lmm_res:

  Linear mixed model results data frame from [`runBinLMM`](runBinLMM.md)

- pairs_use:

  Which pairs to plot

- title:

  Title of the plot

## Value

A `ggplot2` oboject
