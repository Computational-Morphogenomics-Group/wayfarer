# Read aggregated data from multiple bin sizes for the same sample

This function reads from the output of
[`runBinAnalyses`](runBinAnalyses.md). SFE objects of the selected bin
sizes will be loaded.

## Usage

``` r
readBins(dir, sides)
```

## Arguments

- dir:

  Directory where the results are stored. The SFE object for each bin
  size must be in a directory whose name begins with "binx_esda" where x
  is the bin size, such as "bin12_esda".

- sides:

  Numeric vector of bin sizes to read

## Value

A list of SFE objects whose names are the bin sizes
