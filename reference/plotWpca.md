# Plot projections of genes/gene pairs in Wpca space

It's similar to `scater::plotPCA`, but each point is a gene or gene pair
instead of a cell.

## Usage

``` r
plotWpca(wpca_out, dims = 1:2, color_by = NULL, size = 0.5, label = TRUE)
```

## Arguments

- wpca_out:

  Output of [`runWpca`](runWpca.md).

- dims:

  Which PCs to plot, as ingeter indices

- color_by:

  Vector to color the points, must be in the same order as the row names
  of `wpca_out$x`

- size:

  Point size

## Value

A `ggplot2` object. `fact_matrix` will be used when plotting more than 2
PCs.
