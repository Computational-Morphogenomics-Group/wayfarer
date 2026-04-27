# Cluster Lee's L curves

Lee's L has been computed for pairs of genes across spatial scales. This
function reads the results and clusters the curves of Lee's L of each
pair of genes across scales. The
[`approxSilhouette`](https://rdrr.io/pkg/bluster/man/approxSilhouette.html)
function can be used to assess cluster quality.

## Usage

``` r
clusterLeeCurves(
  df,
  hclust_params = list(),
  leiden_params = list(resolution = 0.8, objective_function = "modularity"),
  mat = NULL
)
```

## Arguments

- df:

  Data frame with Lee's L results from
  [`runBinAnalyses`](runBinAnalyses.md), which should have columns pair,
  side, and lee. The data frame can be read into R with
  [`readLeeSamples`](readLeeSamples.md), where a cutoff can be set to
  remove gene pairs with low Lee's L in all bin sizes.

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
gene pair added.
