# Plot pairs of genes in space from multiple SFE objects

This function is useful to visualize pairs of genes whose Lee's L is of
interest.

## Usage

``` r
plotSFEsBiscale(
  sfes,
  feature1,
  feature2,
  df_lee,
  colGeometryName = 1L,
  bbox = NULL,
  palette = "BlueGold",
  exprs_value = "logcounts",
  swap_rownames = NULL,
  title = NULL,
  crop = TRUE,
  show_sizes = TRUE,
  sep = "_",
  ncol = NULL,
  nrow = NULL,
  side_use = NULL,
  widths = NULL,
  heights = NULL,
  design = NULL,
  legend_index = NULL
)
```

## Arguments

- sfes:

  A list of SFE objects, whose names must be the bin sizes.

- feature1:

  First feature to plot

- feature2:

  Second feature to plot

- df_lee:

  Data frame from [`clusterLeeCurves`](clusterLeeCurves.md), where Lee's
  L values are extracted to annotate the plot.

- colGeometryName:

  Name of a `colGeometry` `sf` data frame whose numeric columns of
  interest are to be used to compute the metric. Use `colGeometryNames`
  to look up names of the `sf` data frames associated with cells/spots.

- bbox:

  Named numeric vector specifying a bounding box, either to zoom into a
  smaller area or to show the box itself. The names should be "xmin",
  "xmax", "ymin", and "ymax" in any order. The same bbox will be used
  for all SFE objects in `sfes`.

- palette:

  Name of bivariate palette; see
  [`bi_pal`](https://chris-prener.github.io/biscale/reference/bi_pal.html)
  for complete list of built-in palettes in the \`biscale\` package.

- exprs_value:

  Which assay whose data should be plotted

- swap_rownames:

  Column name of `rowData(object)` to be used to identify features
  instead of `rownames(object)` when labeling plot elements. If not
  found in `rowData`, then rownames of the gene count matrix will be
  used.

- title:

  Title of the entire multi-panel plot

- crop:

  Logical, whether to crop the samples by `bbox`. If `FALSE`, then the
  `bbox` will be shown as a box on the plot to indicate its location.

- show_sizes:

  Logical, whether to show bin sizes in the titles of individual panels.
  If TRUE, then `sfes` must have names that are the bin sizes

- sep:

  Separator between the two gene symbols in `df_lee`

- ncol, nrow:

  The dimensions of the grid to create - if both are `NULL` it will use
  the same logic as
  [facet_wrap()](https://ggplot2.tidyverse.org/reference/facet_wrap.html)
  to set the dimensions

- side_use:

  Which side length whose Lee's L is to be plotted. It must be specified
  when plotting the same bin size from multiple samples.

- widths, heights:

  The relative widths and heights of each column and row in the grid.
  Will get repeated to match the dimensions of the grid. The special
  value of `NA`/`-1null` will behave as `1null` unless a fixed aspect
  plot is inserted in which case it will allow the dimension to expand
  or contract to match the aspect ratio of the content

- design:

  Specification of the location of areas in the layout. Can either be
  specified as a text string or by concatenating calls to
  [`area()`](https://patchwork.data-imaginist.com/reference/area.html)
  together. See the examples for further information on use.

- legend_index:

  Index of position to place the legend. If `NULL`, then the legend will
  not be shown. Unlike in univariate plots, the legend itself is a
  `ggplot2` object

## Value

A `patchwork` object.
