#' Plot Moran's I curves
#'
#' The curves show how Moran's I changes over spatial scale, as indicated by bin
#' size used to aggregate the data.
#'
#' @param df Data frame output from \code{\link{clusterMoranCurves}}, or any
#'   data frame with a column named "gene" for genes, "side" for bin size as in
#'   side length, and "moran" for Moran's I values, and optionally other
#'   categorical columns indicating cluster membership of genes.
#' @param color_by A data frame with column "gene" for gene symbols or IDs and
#'   one other column for color. Alternatively, name of a column in \code{df} to
#'   use for coloring.
#' @param facet_by Name of a categorical or integer column in \code{df} to facet
#'   the plot, not tidyeval.
#' @param show_null Logical, whether to show expected values of Moran's I under
#'   null hypothesis (values are randomly permuted in space) and the interval
#'   within which the null is not to be rejected (p >= 0.05 after Bonferroni
#'   correction accounting for the number of genes and number of sides).
#' @param sfes The list of SFE objects from which \code{df} was made. Only
#'   required when \code{show_null = TRUE}, to extract the number of bins in
#'   order to compute the mean and variance of Moran's I for each bin size.
#' @param color_name Name to show for colors in the legend.
#' @importFrom ggplot2 ggplot aes geom_line scale_x_continuous
#'   scale_color_viridis_c facet_wrap geom_ribbon labs scale_color_manual
#' @importFrom rlang !!! sym
#' @importFrom spdep moran.test
#' @importFrom SpatialFeatureExperiment colGraph
#' @importFrom SingleCellExperiment logcounts
#' @importFrom scales breaks_log
#' @importFrom stats qnorm
#' @return A \code{ggplot} object.
#' @export
plotMoranCurves <- function(df, color_by = NULL, facet_by = NULL,
                            show_null = FALSE, sfes = NULL,
                            color_name = NULL) {
    # I'll do data validation later
    if (!is.null(color_by)) {
        if (is.data.frame(color_by)) {
            names(color_df)[names(color_df) != "gene"] <- "color"
            df <- df |>
                left_join(color_df, by = "gene")
        } else if (is.character(color_by) && color_by %in% names(df)) {
            df <- df |>
                mutate(color = !!!sym(color_by))
        }
        p <- ggplot(df, aes(side, moran, group = gene, color = color))
        if (is.numeric(color_df$color))
            p <- p + scale_color_viridis_c(option = "E")
        else
            p <- p + scale_color_manual(values = Voyager::ditto_colors)
    } else
        p <- ggplot(df, aes(side, moran, group = gene))
    if (show_null) {
        if (is.null(sfes)) {
            stop("sfes must be supplied when show_null = TRUE")
        }
        moran_mean_vars <- tibble(side = names(sfes),
                                  mvs = lapply(sfes, function(x) {
                                      foo <- moran.test(logcounts(x)[1,],
                                                        listw = colGraph(x))
                                      foo$estimate[2:3]
                                  }),
                                  mean = vapply(mvs, function(x) x[1],
                                                FUN.VALUE = numeric(1)),
                                  var = vapply(mvs, function(x) x[2],
                                               FUN.VALUE = numeric(2))) |>
            dplyr::select(-mvs)
        moran_mean_vars <- moran_mean_vars |>
            mutate(th1 = qnorm(0.025/nrow(sfes[[1]])/length(sfes),
                               mean = mean, sd = sqrt(var),
                               lower.tail = TRUE),
                   th2 = qnorm(0.025/nrow(sfes[[1]])/length(sfes),
                               mean = mean, sd = sqrt(var),
                               lower.tail = FALSE))
    }
    p <- p +
        geom_line(alpha = 0.5) +
        scale_x_continuous(transform = "log2", breaks = breaks_log(10, 2))
    if (!is.null(facet_by)) {
        p <- p + facet_wrap(~ !!!sym(facet_by))
    }
    if (show_null) {
        p <- p +
            geom_ribbon(data = moran_mean_vars, aes(side, ymin = th1, ymax = th2),
                        fill = "darkslategray1", color = "cyan", alpha = 0.5)
    }
    p + labs(title = "Moran's I across scales", x = "Bin size (μm)",
             y = "Moran's I", color = color_name)
}

#' Plot gene expression from multiple SFE objects
#'
#' This is useful to show what this gene looks like across different bin sizes.
#'
#' @param sfes A list of SFE objects.
#' @param feature Gene to plot
#' @param bbox Bounding box to plot a subarea
#' @param ... Other arguments to pass to
#'   \code{\link[Voyager]{plotSpatialFeature}}
#' @importFrom Voyager plotSpatialFeature
#' @importFrom ggplot2 ggtitle
#' @importFrom patchwork wrap_plots
#' @return A \code{patchwork} object
#' @export
plotSFEs <- function(sfes, feature, bbox = NULL, ...) {
    ps <- lapply(sfes, function(x) {
        p <- if (ncol(x) > 2e5 && is.null(bbox))
            plotSpatialFeature(x, feature, colGeometryName = "centroids",
                               scattermore = TRUE)
        else plotSpatialFeature(x, feature, bbox = bbox)
        p + ggtitle(paste0("I = ", format(rowData(x)[feature, "moran_sample01"],
                                          digits = 3)))
    })
    wrap_plots(ps)
}
