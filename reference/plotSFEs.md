# Plot gene expression from multiple SFE objects

This is useful to show what this gene looks like across different bin
sizes. When there are more than 200,000 cells or bins and no bbox is
specified, then scattermore is used to speed up plotting.

## Usage

``` r
plotSFEs(
  sfes,
  feature,
  bbox = NULL,
  title = NULL,
  show_sizes = TRUE,
  ncol = NULL,
  nrow = NULL,
  widths = NULL,
  heights = NULL,
  design = NULL,
  crop = TRUE,
  ...
)
```

## Arguments

- sfes:

  A list of SFE objects.

- feature:

  Gene to plot, only one gene

- bbox:

  Bounding box to plot a subarea

- title:

  Title of the entire multi-panel plot

- show_sizes:

  Logical, whether to show bin sizes in the titles of individual panels.
  If TRUE, then `sfes` must have names that are the bin sizes

- ncol, nrow:

  The dimensions of the grid to create - if both are `NULL` it will use
  the same logic as
  [facet_wrap()](https://ggplot2.tidyverse.org/reference/facet_wrap.html)
  to set the dimensions

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

- crop:

  Logical, whether to only plot a cropped area when `bbox` is specified.
  If FALSE, then the bounding box will be shown as a box on the plot
  while the entire tissue section will be shown. This is used to show
  where a bounding box is to provide context along side a plot cropping
  by the bounding box.

- ...:

  Other arguments to pass to
  [`plotSpatialFeature`](https://rdrr.io/pkg/Voyager/man/plotSpatialFeature.html)

## Value

A `patchwork` object
