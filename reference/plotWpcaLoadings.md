# Plot Wpca loadings

The loadings are for each combination of bin size and group. This helps
interpreting variations in multi-scale behaviors in spatial metrics.

## Usage

``` r
plotWpcaLoadings(
  wpca_out,
  dims = 1:4,
  nfeatures = 10,
  balanced = TRUE,
  ncol = NULL
)
```

## Arguments

- wpca_out:

  Output from [`runWpca`](runWpca.md).

- dims:

  Which PCs to plot, as ingeter indices

- nfeatures:

  Total number of features whose loadings are to be plotted

- balanced:

  Logical, whether to plot equal number of features with the most
  positive and most negative loadings.

- ncol:

  Passed to
  [`facet_wrap`](https://ggplot2.tidyverse.org/reference/facet_wrap.html),
  for the number of columns in the facets when plotting multiple PCs.

## Value

A `ggplot2` object
