# Read one bin size from multiple samples

This function reads from the output of
[`runBinAnalyses`](runBinAnalyses.md). Results from one bin size are
read from multiple samples.

## Usage

``` r
readBinSamples(dirs, side)
```

## Arguments

- dirs:

  A character vector of directories, one for each sample. The SFE object
  for each bin size must be in a directory whose name begins with
  "binx_esda" where x is the bin size, such as "bin12_esda".

- side:

  One bin size to read

## Value

A list of SFE objects whose names are base names in `dirs`.
