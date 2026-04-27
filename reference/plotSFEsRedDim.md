# Plot dimension reduction in space for multiple SFE objects

This function not only plots dimension reductions in space, but also
aligns the components for a consistent look because in PCA, the PC's can
be flipped in sign.

## Usage

``` r
plotSFEsRedDim(
  sfes,
  pc = 1,
  reduction = "PCA",
  field = "rotation",
  divergent = TRUE,
  diverge_center = 0,
  bbox = NULL
)
```

## Arguments

- sfes:

  A list of SFE objects, for which the dimension reduction of interest
  has been computed.

- pc:

  Index of principal component or eigenvector or loading vector to plot.
  Facetting is used when length of this argument is greater than 1.

- reduction:

  Name of dimension reduction to plot

- field:

  Field in the attribute of the dimension reduction where variance
  explained or eigenvalues are stored. If this field is not found, then
  the components will not be flipped and flipping is not necessary for
  dimension reductions such as NMF which only produce non-negative
  results.

- divergent:

  Logical, whether a divergent palette should be used.

- diverge_center:

  If `divergent = TRUE`, the center from which the palette should
  diverge. If `NULL`, then not centering.

- bbox:

  A bounding box to specify a smaller region to plot, useful when the
  dataset is large. Can be a named numeric vector with names "xmin",
  "xmax", "ymin", and "ymax", in any order. If plotting multiple
  samples, it should be a matrix with sample IDs as column names and
  "xmin", "ymin", "xmax", and "ymax" as row names. If multiple samples
  are plotted but `bbox` is a vector rather than a matrix, then the same
  bounding box will be used for all samples. You may see points at the
  edge of the geometries if the intersection between the bounding box
  and a geometry happens to be a point there. If `NULL`, then the entire
  tissue is plotted.

## Value

A `patchwork` object
