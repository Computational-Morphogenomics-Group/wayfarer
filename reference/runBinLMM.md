# Use linear mixed model to compare among groups

This function uses a linear mixed model with a spline term to model the
Moran's I or Lee's L curves, so each group has its own intercept and
slope, distinguishing between variability within and between groups.
Then the entire random effects term or the random slope is dropped, and
a likelihood ratio test is used to compared the full model and the
reduced one to see if the random effects or random slope is significant.
If it is significant, then it indicates that the curves differ among
groups, though it does not indicate which groups are different.

## Usage

``` r
runBinLMM(df_res, degree = 2, BPPARAM = SerialParam())
```

## Arguments

- df_res:

  A data frame with Moran's I or Lee's L across bin sizes and samples,
  such as from [`readMoranSamples`](readMoranSamples.md) or
  [`readLeeSamples`](readLeeSamples.md). The columns should be renamed
  so that the gene or gene pair column is named "feature", the Moran's I
  or Lee's L values column is named "value", and the group column from
  additional sample info is named "group".

- degree:

  degree of the piecewise polynomial—default is `3` for cubic splines.

- BPPARAM:

  A `bpparam` object to parallelize computation over features.

## Value

A data frame with p-values and adjusted p-values for each feature
