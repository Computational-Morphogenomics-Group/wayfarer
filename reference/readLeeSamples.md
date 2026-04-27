# Read multi-sample Lee's L results from multiple samples

[`runBinAnalyses`](runBinAnalyses.md) should write Lee's L for all genes
in all bin sizes to a file `df_lee.csv`. This function reads this file
from multiple samples and concatenates them.

## Usage

``` r
readLeeSamples(dirs, sample_info = NULL, cutoff = 0L)
```

## Arguments

- dirs:

  A character vector of directories, one for each sample. There must be
  a `df_lee.csv` file in each directory.

- sample_info:

  Optional data frame with more info about each sample. There must be a
  column called "sample". Such info is helpful when plotting.

- cutoff:

  Gene pairs whose absolute values of Lee's L is below this cutoff in
  all bin sizes and all samples will not be included in the output. This
  will remove genes that are not spatially correlated at any length
  scale and any sample.

## Value

A data frame with columns moran, side, sample, and optionally other
columns in `sample_info`.
