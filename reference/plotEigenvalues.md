# Plot variance explained across bin sizes

This function plots the iconic PCA elbow plot for variance explained or
eigenvalue for multiple SCE objects, distinguished by color.

## Usage

``` r
plotEigenvalues(
  sfes,
  npcs = 20,
  reduction = "PCA",
  field = "varExplained",
  xlab = "Principal component",
  ylab = "Variance explained"
)
```

## Arguments

- sfes:

  A list of SFE (SCE and SPE are fine) objects with the dimension. Names
  of the list are bin sizes. reduction of interest computed and stored
  in `reducedDims`.

- npcs:

  Number of components whose variance explained or eigenvalues are to be
  plotted.

- reduction:

  Name of the dimension reduction of interest

- field:

  Field in the attribute of the dimension reduction where variance
  explained or eigenvalues are stored.

- xlab:

  X-axis label

- ylab:

  Y-axis label

## Value

A `ggplot` object
