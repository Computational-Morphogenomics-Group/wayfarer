# Plot Moran's I of bin projection in PCA space

Basically see how spatially structured each principal component is and
how it may relate to bin size.

## Usage

``` r
plotMoranPCs(sfes, npcs = 20, reduction = "PCA")
```

## Arguments

- sfes:

  A list of SFE objects, with Moran's I computed for bin projections in
  the dimension reduction of interest.

- npcs:

  Number of components whose variance explained or eigenvalues are to be
  plotted.

- reduction:

  Name of the dimension reduction of interest

## Value

A `ggplot` object
