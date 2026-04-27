# Plot Lee's L curves

The curves show how Lee's L changes over spatial scale for each pair of
genes, as indicated by bin size used to aggregate the data.

## Usage

``` r
plotLeeCurves(df, facet_by = NULL, show_median = FALSE, sample_n = NULL)
```

## Arguments

- df:

  Data frame output from [`clusterLeeCurves`](clusterLeeCurves.md), or
  any data frame with a column named "pair" for gene pairs, "side" for
  bin size as in side length, and "lee" for Lee's L values, and
  optionally other categorical columns indicating cluster membership of
  genes.

- facet_by:

  Name of a categorical or integer column in `df` to facet the plot, not
  tidyeval.

- show_median:

  Logical, whether to plot the median of clusters, only used when
  `facet_by` is specified.

- sample_n:

  A smaller number of gene pairs to plot. Often the number of gene pairs
  is large so downsampling makes the plot more readable.

## Value

A `ggplot2` object.
