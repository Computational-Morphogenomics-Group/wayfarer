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
runBinLMM(df_res, degree = 2, BPPARAM = SerialParam(), save_models = FALSE)
```

## Arguments

- df_res:

  A data frame with Moran's I or Lee's L across bin sizes and samples,
  such as from [`readMoranSamples`](readMoranSamples.md) or
  [`readLeeSamples`](readLeeSamples.md). The columns should be renamed
  so that the gene or gene pair column is named "feature", the Moran's I
  or Lee's L values column is named "value", and the group column from
  additional sample info is named "group". There should also be a column
  "weights" due to heteroscadiscity.

- degree:

  degree of the piecewise polynomial—default is `3` for cubic splines.

- BPPARAM:

  A `bpparam` object to parallelize computation over features.

- save_models:

  Logical, whether to save the fitted models.

## Value

A data frame with p-values, adjusted p-values, and log likelihood ratios
for each feature. If `save_models = TRUE`, then there will also be a
list column for the fitted models and a list column `data` for the data
used to fit the models, which can be used to get the predicted curves.
