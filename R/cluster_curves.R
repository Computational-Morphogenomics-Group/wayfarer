.cluster_mat <- function(mat, hclust_params, leiden_params) {
    tibble(gene = rownames(mat),
           hclust = clusterRows(mat, BLUSPARAM = do.call(HclustParam,
                                                         hclust_params)),
           leiden = clusterRows(mat,
                      BLUSPARAM = NNGraphParam(cluster.fun = "leiden",
                                               cluster.args = leiden_params)))
}



#' Cluster Moran's I curves
#'
#' Moran's I has been computed for the same genes across different spatial
#' scales. This function gets the Moran's I values for all genes and resolutions
#' and clusters the patterns with which Moran's I's vary through scales, as a
#' way to cluster genes. Leiden and hierarchical clustering are performed on the
#' Moran's I values themselves and on the diff between adjacent scales. The
#' \code{\link[bluster]{approxSilhouette}} function can be used to assess
#' cluster quality.
#'
#' @param df Data frame from \code{\link{runBinAnalyses}} for Moran's I, with
#'   columns moran, gene, and side.
#' @param hclust_params Hierarchical clustering parameters, passed to
#'   \code{\link[bluster]{HclustParam}}.
#' @param leiden_params Leiden clustering parameters, passed to
#'   \code{\link[igraph]{cluster_leiden}}.
#' @param mat Matrix with column as bin side lengths and rows as genes. If NULL,
#'   then it will be made from \code{df}.
#' @return The same data frame in the input but with cluster assignment of each
#' gene added.
#' @importFrom tibble tibble rownames_to_column column_to_rownames
#' @importFrom dplyr select rename left_join mutate filter if_any group_by
#'   summarize
#' @importFrom tidyr unnest pivot_wider pivot_longer unite
#' @importFrom SingleCellExperiment rowData
#' @importFrom bluster clusterRows HclustParam NNGraphParam
#' @export
clusterMoranCurves <- function(df, hclust_params = list(),
                               leiden_params = list(resolution = 0.8,
                                                    objective_function = "modularity"),
                               mat = NULL) {
    # df is in long form, mat is in wide form
    if (is.null(mat)) {
        mat <- df |>
            select(side, gene, moran) |>
            pivot_wider(names_from = side, values_from = moran) |>
            column_to_rownames("gene") |>
            as.matrix()
    }
    diffs <- apply(mat, 1, diff) |> t()
    df_clust <- .cluster_mat(mat, hclust_params, leiden_params)
    df_clust_diff <- .cluster_mat(diffs, hclust_params, leiden_params) |>
        dplyr::select(-gene, hclust_diffs = hclust, leiden_diffs = leiden)
    df_clust <- cbind(df_clust, df_clust_diff)
    df |> left_join(df_clust, by = "gene")
}

#' Cluster Lee's L curves
#'
#' Lee's L has been computed for pairs of genes across spatial scales. This
#' function reads the results and clusters the curves of Lee's L of each pair of
#' genes across scales. The \code{\link[bluster]{approxSilhouette}} function can
#' be used to assess cluster quality.
#'
#' @inheritParams clusterMoranCurves
#' @param df Data frame with Lee's L results from \code{\link{runBinAnalyses}},
#'   which should have columns pair, side, and lee. The data frame can be read
#'   into R with \code{\link{readLeeSamples}}, where a cutoff can be set to
#'   remove gene pairs with low Lee's L in all bin sizes.
#' @return The same data frame in the input but with cluster assignment of each
#'   gene pair added.
#' @export
clusterLeeCurves <- function(df, hclust_params = list(),
                             leiden_params = list(resolution = 0.8,
                                                  objective_function = "modularity"),
                             mat = NULL) {
    clusterMoranCurves(df_lee |> rename(gene = pair, moran = lee), hclust_params, leiden_params,
                    mat = mat) |>
        rename(lee = moran, pair = gene)
}
