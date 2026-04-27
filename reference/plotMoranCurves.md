# Plot Moran's I curves

The curves show how Moran's I changes over spatial scale, as indicated
by bin size used to aggregate the data. Despite the name, this function
can be applied to other univariate spatial results such as Geary's C.

## Usage

``` r
plotMoranCurves(
  df,
  name = "moran",
  color_by = NULL,
  facet_by = NULL,
  show_null = FALSE,
  mean_vars = NULL,
  color_name = NULL,
  title_name = "Moran's I"
)
```

## Arguments

- df:

  Data frame output from [`clusterMoranCurves`](clusterMoranCurves.md),
  or any data frame with a column named "gene" for genes, "side" for bin
  size as in side length, and "moran" for Moran's I values, and
  optionally other categorical columns indicating cluster membership of
  genes.

- name:

  Column name in `df` with the spatial results.

- color_by:

  A data frame with column "gene" for gene symbols or IDs and one other
  column for color. Alternatively, name of a column in `df` to use for
  coloring.

- facet_by:

  Name of a categorical or integer column in `df` to facet the plot, not
  tidyeval.

- show_null:

  Logical, whether to show expected values of Moran's I under null
  hypothesis (values are randomly permuted in space) and the interval
  within which the null is not to be rejected (p \>= 0.05 after
  Bonferroni correction accounting for the number of genes and number of
  sides).

- mean_vars:

  Data frame for means and variances of Moran's I under null hypothesis
  from [`getMoranMeanVar`](getMoranMeanVar.md), required for
  `show_null = TRUE`.

- color_name:

  Name to show for colors in the legend.

- title_name:

  Name of the metric to put on the plot title.

## Value

A `ggplot` object.
