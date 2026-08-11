#' Run weighted PCA on fitted curves per stage
#'
#' To further interpret variabilities in the multi-scale behaviors. This is to
#' be done on fitted curves from \code{\link{getPredCurves}} to better elucidate
#' trends when there is more variability among individual samples within the
#' same group. This function is a thin wrapper around
#' \code{\link[WCluster]{Wpca}}, doing some data wrangling to convert the input
#' data frame into a matrix as input to that function.
#'
#' @param df A data frame, output from \code{\link{getPredCurves}} or
#' \code{\link{getEMM}}.
#' @param df_mean_var A dara frame with a column \code{var} for the variance of
#'   the spatial metric for each bin size and each sample, and a column
#'   \code{group} for the groups being compared in the linear mixed models. For
#'   example, output from \code{\link{getMoranMeanVar}}, but with groups added.
#'   See vignette. This data frame is used to compute the weights so the larger
#'   bin sizes and samples with fewer cells get less weight. The weight is
#'   1/var. The maximum weight among all samples in a group is used.
#' @param ... More arguments to pass to \code{\link[WCluster]{Wpca}}, except for
#'   \code{x} and \code{wcol}.
#' @return Same list output as \code{\link[WCluster]{Wpca}}, except that row
#'   names have been added to the item "rotation" to facilitate further
#'   visualization.
#' @importFrom rlang check_installed
#' @export
runWpca <- function(df, weights_col, ...) {
    check_installed("WCluster")
    mat <- df |>
        arrange(group, side) |>
        unite("id", side, group) |>
        pivot_wider(id_cols = feature, names_from = id, values_from = pred) |>
        column_to_rownames("feature") |>
        as.matrix()
    wts <- df[[weights_col]]
    wts <- wts/max(wts)
    out <- WCluster::Wpca(mat, wcol = wts, ...)
    rownames(out$rotation) <- colnames(mat)
    out
}

#' Plot Wpca loadings
#'
#' The loadings are for each combination of bin size and group. This helps
#' interpreting variations in multi-scale behaviors in spatial metrics.
#'
#' @param wpca_out Output from \code{\link{runWpca}}.
#' @param dims Which PCs to plot, as ingeter indices
#' @param nfeatures Total number of features whose loadings are to be plotted
#' @param balanced Logical, whether to plot equal number of features with the
#'   most positive and most negative loadings.
#' @param ncol Passed to \code{\link[ggplot2]{facet_wrap}}, for the number of
#'   columns in the facets when plotting multiple PCs.
#' @return A \code{ggplot2} object
#' @importFrom dplyr any_of slice_max slice_min
#' @importFrom tibble rownames_to_column
#' @export
plotWpcaLoadings <- function(wpca_out, dims = 1:4, nfeatures = 10,
                             balanced = TRUE, ncol = NULL) {
    df_loadings_pred <- as.data.frame(wpca_out$rotation)
    names(df_loadings_pred) <- paste0("PC", seq_len(ncol(df_loadings_pred)))
    df_loadings_pred <- rownames_to_column(df_loadings_pred, "id")
    pcs <- paste0("PC", dims)
    df_plt <- df_loadings_pred |>
        select(id, any_of(pcs)) |>
        pivot_longer(-id, names_to = "PC", values_to = "value") |>
        unite(col = "lab", id, PC, sep = "___", remove = FALSE)
    if (balanced) {
        nn <- floor(nfeatures/2)
        df_plt1 <- df_plt |>
            group_by(PC) |>
            slice_max(value, n = nn, with_ties = FALSE)
        df_plt2 <- df_plt|>
            group_by(PC) |>
            slice_min(value, n = nn, with_ties = FALSE)
        df_plt_use <- bind_rows(df_plt1, df_plt2)
    } else {
        df_plt_use <- df_plt |>
            group_by(PC) |>
            slice_max(abs(value), n = nfeatures)
    }
    df_plt_use <- df_plt_use |>
        mutate(lab = reorder(lab, value))
    reg <- "___.+$"
    ggplot(df_plt_use, aes(value, lab)) +
        geom_segment(aes(yend = lab), xend = 0, show.legend = FALSE) +
        geom_vline(xintercept = 0, linetype = 2) +
        geom_point(color = "blue") +
        scale_y_discrete(labels = function(x) gsub(reg, "", x)) +
        labs(x = "Loading", y = "Feature") +
        facet_wrap(~ PC, scales = "free_y", ncol = ncol)
}

#' Plot projections of genes/gene pairs in Wpca space
#'
#' It's similar to \code{scater::plotPCA}, but each point is a gene or gene pair
#' instead of a cell.
#'
#' @inheritParams plotWpcaLoadings dims
#' @param wpca_out Output of \code{\link{runWpca}}.
#' @param color_by Vector to color the points, must be in the same order as the
#' row names of \code{wpca_out$x}
#' @param size Point size
#' @return A \code{ggplot2} object. \code{\link[ggforce]{fact_matrix}} will be
#'   used when plotting more than 2 PCs.
#' @export
plotWpca <- function(wpca_out, dims = 1:2, color_by = NULL, size = 0.5,
                     label = TRUE) {
    if (length(dims > 2L)) check_installed("ggforce")
    if (label) check_installed("ggrepel")
    df_wpca_pred <- as.data.frame(wpca_out$x)
    names(df_wpca_pred) <- paste0("PC", seq_len(ncol(df_wpca_pred)))
    df_wpca_pred <- rownames_to_column(df_wpca_pred, "feature")
    if (!is.null(color_by))
        df_wpca_pred$color_by <- color_by
    if (!is.null(color_by))
        p <- ggplot(df_wpca_pred, aes(color = color_by))
    else p <- ggplot(df_wpca_pred)
    pcs <- paste0("PC", dims)
    if (length(dims) > 2L) {
        p <- p +
            geom_point(aes(x = .panel_x, y = .panel_y), size = size)
        if (label)
            p <- p + ggrepel::geom_label_repel(aes(x = .panel_x, y = .panel_y, label = feature))
        if (!is.null(color_by))
            p <- p + ggforce::geom_autodensity(fill = NA, linewidth = 0.5, position = "identity")
        else
            p <- p + ggforce::geom_autodensity(fill = NA, linewidth = 0.5, color = "black")
        p <- p +
            ggforce::facet_matrix(vars(all_of(pcs)), layer.diag = if (label) 3 else 2)
    } else {
        p <- p + geom_point(aes(.data[[pcs[1]]], .data[[pcs[2]]]))
        if (label)
            p <- p + ggrepel::geom_label_repel(aes(.data[[pcs[1]]], .data[[pcs[2]]],
                                          label = feature))
    }
    p
}
