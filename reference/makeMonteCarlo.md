# Create Monte Carlo simulations with no spatial structure

The purpose of this function is to check if edge effect is causing
spurious changes in Moran's I and Lee's L across scales.This function
uses geometries from existing aggregations to create Monte Carlo
simulations of gene expression with no spatial structure. Poisson
samples with mean informed by the total counts of a given gene are taken
on small bins, which are then aggregated. Because given complete spatial
randomness (CSR), the transcript spots should be a homogeneous Poisson
point process within the tissue boundary. The number of counts within
any given area follows a Poisson distribution whose mean is the mean of
the Poisson point process multiplied by the size of that area.
Therefore, to make computation faster, we can take Poisson samples with
means scaled by the area of bins in tissue.

## Usage

``` r
makeMonteCarlo(
  in_path,
  out_path,
  tissue_boundary,
  quantiles = c(0.25, 0.5, 0.75),
  nsim = 100,
  BPPARAM = SerialParam(),
  queen = FALSE
)
```

## Arguments

- in_path:

  Path with SFE objects written to subdirectories with
  [`runBinAnalyses`](runBinAnalyses.md).

- out_path:

  Output directory

- tissue_boundary:

  `sf` or `sfc` for tissue boundary. See
  [`getTissueBoundaryConcave`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/getTissueBoundaryConcave.html)
  and
  [`getTissueBoundaryImg`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/getTissueBoundaryImg.html)
  on obtaining the tissue boundary.

- quantiles:

  Vector of quantiles of total counts of genes to use for the intensity
  of the Poisson point process, because edge effect affects genes with
  different expression levels differently.

## Value

`out_path`; the output is written to disk there.

## Details

This function creates uncorrelated non-spatial patterns. They are stored
as SFE objects and written to disk. Then spatial metrics can be computed
on those patterns. This can be used to check if edge effect is causing
artifacts in the multi-scale analysis, and to compute Monte Carlo
p-values for the linear mixed models used to compare across biological
conditions where the log likelihood ratio test p-values are unreliable
due to correlation of spatial metrics between adjacent bin sizes and
heteroscedasticity, i.e. variance is higher for larger bins because
there are fewer bins.

Because edge effect artifacts affect genes with different expression
levels differently, `nsim` CSR patterns are generated for each quantile
of expression levels.
