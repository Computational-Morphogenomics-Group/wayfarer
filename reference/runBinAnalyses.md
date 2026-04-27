# Run basic analyses on aggregated SFE objects of all bin sizes

This function runs log normalization, bin-level QC, PCA, adjacency
graph, Moran's I, and Lee's L for all binned SFE objects in a directory,
given that bins overlapping too little with tissue have been removed.

## Usage

``` r
runBinAnalyses(
  dir,
  out_path,
  tissue_geometry,
  min_props = 0.9,
  quantiles = NULL,
  ncomponents = 30,
  queen = FALSE,
  zero.policy = TRUE,
  p.adjust.method = "BH",
  neg_regex = c("^NegPrb", "^Blank", "^BLANK", "^NegControl", "Unassigned"),
  BPPARAM = SerialParam(),
  ...
)
```

## Arguments

- dir:

  Directory where the binned outputs are located. Output of each bin
  size must be in a directory named "binx", such as "bin12" for 12
  micron bins.

- out_path:

  Directory to write the output

- tissue_geometry:

  Either tissue boundary (with holes if present in tissue) or cell
  segmentation polygons. Area of overlap of each bin with this geometry
  will be computed.

- min_props:

  Minimum proportion of each bin overlapping tissue; bins that don't
  overlap enough will be removed to mitigate edge effect. This should be
  a numeric vector of the same length as the number of bin sizes and can
  differ for different bin sizes. If length 1, then the same value will
  be applied to all bin sizes. Otherwise, it must have names
  corresponding to bin sizes.

- quantiles:

  Numeric vector of quantiles of area of bin overlapping tissue, can
  differ for different bin sizes. If not NULL, this will supercede
  `min_props`.

- ncomponents:

  Number of components to compute for PCA.

- queen:

  if TRUE, a single shared boundary point meets the contiguity
  condition, if FALSE, more than one shared point is required; note that
  more than one shared boundary point does not necessarily mean a shared
  boundary line

- zero.policy:

  default `attr(listw, "zero.policy")` as set when `listw` was created,
  if attribute not set, use global option value; if TRUE assign zero to
  the lagged value of zones without neighbours, if FALSE assign NA

- p.adjust.method:

  Method to correct for multiple testing, passed to
  [`p.adjustSP`](https://r-spatial.github.io/spdep/reference/p.adjustSP.html).
  Methods allowed are in
  [`p.adjust.methods`](https://rdrr.io/r/stats/p.adjust.html).

- neg_regex:

  Character vector of regex patterns indicating that a feature is for
  negative control. Features matching any of these patterns will be
  removed prior to analyses. Counts of these features are kept in the
  aggregated input data in `dir`.

- BPPARAM:

  A `BiocParallelParam` object specifying whether and how computing the
  metric for numerous genes shall be parallelized.

- ...:

  Other arguments passed to S4 method (for convenience wrappers like
  `calculateMoransI`) or method used to compute metrics as specified by
  the argument `type` (as in more general functions like
  `calculateUnivariate`). See documentation of functions with the same
  name as specified in `type` in the `spdep` package for the method
  specific arguments. For variograms, see
  [`.variogram`](https://rdrr.io/pkg/Voyager/man/variogram-internal.html).

## Value

Invisibly `out_path`; the log normalization, QC, PCA, adjacency graph,
and Moran's I will be stored in the SFE object; the SFE object with the
results will be written to `out_path` with directory names "binx_esda".
A data frame for Moran's I and Lee's L across bin sizes will be written
to `out_path` as CSV files. The "esda" means that exploratory spatial
data analysis (ESDA) has been performed.
