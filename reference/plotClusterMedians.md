# Plot cluster median curves

To better juxtapose different cluster patterns in the same plot.

## Usage

``` r
plotClusterMedians(df, cluster_col)
```

## Arguments

- df:

  Data frame output from either
  [`clusterLeeCurves`](clusterLeeCurves.md) or
  [`clusterMoranCurves`](clusterMoranCurves.md). There must be a column
  named either "moran" or "lee" and a column called "side" in addition
  to the column specified in `cluster_col`.

- cluster_col:

  Name of a categorical column in `df` to use as cluster labels.

## Value

A `ggplot2` object.
