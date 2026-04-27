# Get the mean and variance of Moran's I

This function computes the mean and variance of Moran's I under the null
hypothesis of no spatial autocorrelation, in the context of plotting
Moran's I curves over bin sizes in
[`plotMoranCurves`](plotMoranCurves.md).

## Usage

``` r
getMoranMeanVar(sfes, name = c("moran", "geary"))
```

## Arguments

- sfes:

  A list of SFE objects of different bin sizes from the same sample (can
  be read with [`readBins`](readBins.md)), required to compute and plot
  the mean and variance of Moran's I (or Geary's C). The list must have
  names corresponding to bin sizes.

- name:

  Whether mean and variance are to be computed for Moran's I (moran) or
  Geary's C (geary).

## Value

A data frame with columns side, mean, var, th1, and th2. Column th1 is
the 2.5 th2 is the 97.5 genes and number of bin sizes so the interval is
rather conservative.

## See also

\[spdep::moran.test()\]
