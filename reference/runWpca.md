# Run weighted PCA on fitted curves per stage

To further interpret variabilities in the multi-scale behaviors. This is
to be done on fitted curves from [`getPredCurves`](getPredCurves.md) to
better elucidate trends when there is more variability among individual
samples within the same group. This function is a thin wrapper around
[`Wpca`](https://rdrr.io/pkg/WCluster/man/Wpca.html), doing some data
wrangling to convert the input data frame into a matrix as input to that
function.

## Usage

``` r
runWpca(df, weights_col, ...)
```

## Arguments

- df:

  A data frame, output from [`getPredCurves`](getPredCurves.md) or
  `getEMM`.

- ...:

  More arguments to pass to
  [`Wpca`](https://rdrr.io/pkg/WCluster/man/Wpca.html), except for `x`
  and `wcol`.

- df_mean_var:

  A dara frame with a column `var` for the variance of the spatial
  metric for each bin size and each sample, and a column `group` for the
  groups being compared in the linear mixed models. For example, output
  from [`getMoranMeanVar`](getMoranMeanVar.md), but with groups added.
  See vignette. This data frame is used to compute the weights so the
  larger bin sizes and samples with fewer cells get less weight. The
  weight is 1/var. The maximum weight among all samples in a group is
  used.

## Value

Same list output as
[`Wpca`](https://rdrr.io/pkg/WCluster/man/Wpca.html), except that row
names have been added to the item "rotation" to facilitate further
visualization.
