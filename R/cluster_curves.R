.cluster_mat <- function(mat, hclust_params, leiden_params) {
    tibble(gene = rownames(mat),
           hclust = clusterRows(mat, BLUSPARAM = do.call(HclustParam,
                                                         hclust_params)),
           leiden = clusterRows(mat,
                      BLUSPARAM = NNGraphParam(cluster.fun = "leiden",
                                               cluster.args = leiden_params)))
}

.cluster_curves <- function(df, hclust_params, leiden_params, mat = NULL) {
    # df is in long form, mat is in wide form
    if (is.null(mat)) {
        mat <- df |>
            select(side, gene, value) |>
            pivot_wider(names_from = side, values_from = value) |>
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
#' @param sfes A list of \code{SpatialFeatureExperiment} objects which have the
#'   same genes, same sample IDs, and with Moran's I computed for the genes. The
#'   names of the list must be the bin size.
#' @param hclust_params Hierarchical clustering parameters, passed to
#'   \code{\link[bluster]{HclustParam}}.
#' @param leiden_params Leiden clustering parameters, passed to
#'   \code{\link[igraph]{cluster_leiden}}.
#' @return A data frame.
#' @importFrom tibble tibble rownames_to_column column_to_rownames
#' @importFrom dplyr select rename left_join mutate filter if_any group_by
#'   summarize
#' @importFrom tidyr unnest pivot_wider pivot_longer unite
#' @importFrom SingleCellExperiment rowData
#' @importFrom bluster clusterRows HclustParam NNGraphParam
#' @export
clusterMoranCurves <- function(sfes, hclust_params = list(),
                               leiden_params = list(resolution = 0.8,
                                        objective_function = "modularity")) {
    # Add data validation later
    sides <- as.integer(names(sfes))
    df_moran <- tibble(side = sides,
                       morans = lapply(sfes, function(x) {
                           rowData(x) |> as.data.frame() |>
                               rownames_to_column() |>
                               dplyr::select(gene = rowname,
                                             value = moran_sample01)
                           # OK, what if the sample_id is something else?
                       })) |>
        unnest(cols = morans)
    .cluster_curves(df_moran, hclust_params, leiden_params) |>
        rename(moran = value)
}

.get_pairs_df <- function(lees) {
    nr <- nrow(lees[[1]])
    inds_df <- tibble(j = unlist(lapply(seq_len(nr-1), function(x) rep(x+1, times = x))),
                      i = unlist(lapply(seq_len(nr-1), seq_len)))
    inds_df <- inds_df |>
        mutate(gene1 = rownames(lees[[1]])[i],
               gene2 = rownames(lees[[1]])[j]) |>
        unite("pair", gene1, gene2, sep = "_")
    inds <- upper.tri(lees[[1]]) # indices
    uts <- lapply(lees, function(l) l[inds])
    names(uts) <- names(lees)
    uts <- as.data.frame(uts, optional = TRUE)
    uts$pair <- inds_df$pair
    uts
}
#' Cluster Lee's L curves
#'
#' Lee's L has been computed for pairs of genes across spatial scales. This
#' function reads the results and clusters the curves of Lee's L of each pair of
#' genes across scales. The \code{\link[bluster]{approxSilhouette}} function can
#' be used to assess cluster quality.
#'
#' @inheritParams clusterMoranCurves
#' @param dir Directory where Lee's L results are stored. File names must have
#'   the pattern "bin<x>_lee.rds", such as "bin12_lee.rds".
#' @param cutoff Only gene pairs with Lee's L greater than the cutoff at at
#'   least one bin size are used for clustering, due to the potentially large
#'   number of pairs.
#' @param sides Numeric vector of bin sizes whose results are to be read.
#' @return A data frame.
#' @export
clusterLeeCurves <- function(dir, sides, cutoff = 0.2, hclust_params = list(),
                             leiden_params = list(resolution = 0.8,
                                        objective_function = "modularity")) {
    fns <- paste0("bin", sides, "_lee.rds")
    lees <- lapply(fns, readRDS)
    names(lees) <- sides
    lees <- lapply(lees, as.matrix)
    rns <- rownames(lees[[1]])
    lees <- lapply(lees, function(x) x[rns, rns])
    df_lee <- .get_pairs_df(lees)
    df_lee <- df_lee |>
        mutate(any_corr = if_any(-pair, ~ abs(.x) > cutoff)) |>
        filter(any_corr) |>
        dplyr::select(-any_corr)
    lee_mat <- df_lee |>
        column_to_rownames(var = "pair") |>
        as.matrix()
    df_lee <- df_lee |>
        pivot_longer(-pair, names_to = "side", values_to = "value") |>
        mutate(side = as.integer(side))
    .cluster_curves(df_lee |> rename(gene = pair), hclust_params, leiden_params,
                    mat = lee_mat) |>
        rename(lee = value, pair = gene)
}

# Actually, shall I write a function to run all the basic analyses for a list of SFEs at once?
# The manual part really is choosing a threshold of proportion of bin in cells
# The report is mostly just copy and paste. OK, so functions to:
# 1. Make the bin aggregated SFEs, save to disk with alabaster.sfe
# 2. Do whatever is in the notebooks, from filtering bins to running ESDA; I have
# rules of thumbs for the slightly manual parts.
# 3. Generate the Rmd reports from the analyses; I ran those notebooks manually because
# I wanted to see and comment on the plots. Oh right, the bboxes must be chosen
# manually for now. There might be a computational way to choose them based on
# concordex results.

