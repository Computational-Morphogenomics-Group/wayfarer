# Plot gene loadings of a given principal component across scales

This can visualize how the meanings of principal components change with
spatial scales.

## Usage

``` r
plotLoadingCurves(
  sfes,
  pc = 1,
  reduction = "PCA",
  field = "rotation",
  color_by = c("loading", "moran"),
  color_ind = 1,
  flip = FALSE,
  ...
)
```

## Arguments

- sfes:

  A list of SFE objects with names as bin sizes

- pc:

  Index of principal component or eigenvector or loading vector to plot.
  Facetting is used when length of this argument is greater than 1.

- reduction:

  Name of dimension reduction to plot

- field:

  Name of the attribute of the dimension reduction that holds the gene
  loadings

- color_by:

  Whether to color the curves by Moran's I or loading values. Can only
  color by Moran's I when multiple PCs are plotted because each PC has
  different gene loadings.

- color_ind:

  Since one curve can only have one color, choose the index of the SFE
  object in `sfes` whose values are used for coloring.

- flip:

  Logical, whether to flip the signs of the loadings.

- ...:

  Arguments passed to
  [`facet_wrap`](https://ggplot2.tidyverse.org/reference/facet_wrap.html).

## Value

A `ggplot` object
