# Read multi-scale Moran's I results from multiple samples

[`runBinAnalyses`](runBinAnalyses.md) should write Moran's I for all
genes in all bin sizes to a file `df_moran.csv`. This function reads
this file from multiple samples and concatenates them.

## Usage

``` r
readMoranSamples(dirs, sample_info = NULL)
```

## Arguments

- dirs:

  A character vector of directories, one for each sample. There must be
  a `df_moran.csv` file in each directory.

- sample_info:

  Optional data frame with more info about each sample. There must be a
  column called "sample". Such info is helpful when plotting.

## Value

A data frame with columns moran, side, sample, and optionally other
columns in `sample_info`.
