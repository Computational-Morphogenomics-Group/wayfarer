# Cluster Moran's I curves

Moran's I has been computed for the same genes across different spatial
scales. This function gets the Moran's I values for all genes and
resolutions and clusters the patterns with which Moran's I's vary
through scales, as a way to cluster genes. Leiden and hierarchical
clustering are performed on the Moran's I values themselves and on the
diff between adjacent scales. The
[`approxSilhouette`](https://rdrr.io/pkg/bluster/man/approxSilhouette.html)
function can be used to assess cluster quality.

## Usage

``` r
clusterMoranCurves(
  df,
  hclust_params = list(),
  leiden_params = list(resolution = 0.8, objective_function = "modularity"),
  mat = NULL
)
```

## Arguments

- df:

  Data frame from [`runBinAnalyses`](runBinAnalyses.md) for Moran's I,
  with columns moran, gene, and side.

- hclust_params:

  Hierarchical clustering parameters, passed to
  [`HclustParam`](https://rdrr.io/pkg/bluster/man/HclustParam-class.html).

- leiden_params:

  Leiden clustering parameters, passed to
  [`cluster_leiden`](https://r.igraph.org/reference/cluster_leiden.html).

- mat:

  Matrix with column as bin side lengths and rows as genes. If NULL,
  then it will be made from `df`.

## Value

The same data frame in the input but with cluster assignment of each
gene added.
