# Get fitted curves per group

After running [`runBinLMM`](runBinLMM.md) with `save_models = TRUE`, use
this function to get the predicted curves for each group.

## Usage

``` r
getPredCurves(df_res)
```

## Arguments

- df_res:

  Data frame output from [`runBinLMM`](runBinLMM.md), with list columns
  "model" and "data".

## Value

A data frame with column "pred" for the predicted curves.
