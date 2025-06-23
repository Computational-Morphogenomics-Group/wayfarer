ditto_colors <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00",
                  "#CC79A7", "#666666", "#AD7700", "#1C91D4", "#007756", "#D5C711",
                  "#005685", "#A04700", "#B14380", "#4D4D4D", "#FFBE2D", "#80C7EF",
                  "#00F6B3", "#F4EB71", "#06A5FF", "#FF8320", "#D99BBD", "#8C8C8C",
                  "#FFCB57", "#9AD2F2", "#2CFFC6", "#F6EF8E", "#38B7FF", "#FF9B4D",
                  "#E0AFCA", "#A3A3A3", "#8A5F00", "#1674A9", "#005F45", "#AA9F0D",
                  "#00446B", "#803800", "#8D3666", "#3D3D3D")

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
#' @importFrom ggplot2 ggplot aes geom_line scale_x_continuous ggtitle
#'   scale_color_viridis_c facet_wrap geom_ribbon labs scale_color_manual
#' @importFrom rlang .data
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
    if (!is.null(facet_by)) {
        df <- df |>
            mutate(facet = .data[[facet_by]])
    }
    if (!is.null(color_by)) {
        if (is.data.frame(color_by)) {
            names(color_by)[names(color_by) != "gene"] <- "color"
            df <- df |>
                left_join(color_by, by = "gene")
        } else if (is.character(color_by) && color_by %in% names(df)) {
            df <- df |>
                mutate(color = .data[[color_by]])
        }
        p <- ggplot(df) +
            geom_line(aes(side, moran, group = gene, color = color), alpha = 0.5)
        if (is.numeric(df$color))
            p <- p + scale_color_viridis_c(option = "E")
        else
            p <- p + scale_color_manual(values = ditto_colors, na.value = "gray80")
    } else
        p <- ggplot(df) +
            geom_line(aes(side, moran, group = gene), alpha = 0.5)
    if (show_null) {
        if (is.null(sfes)) {
            stop("sfes must be supplied when show_null = TRUE")
        }
        moran_mean_vars <- tibble(side = as.integer(names(sfes)),
                                  mvs = lapply(sfes, function(x) {
                                      foo <- moran.test(logcounts(x)[1,],
                                                        listw = colGraph(x))
                                      foo$estimate[2:3]
                                  }),
                                  mean = vapply(mvs, function(x) x[1],
                                                FUN.VALUE = numeric(1)),
                                  var = vapply(mvs, function(x) x[2],
                                               FUN.VALUE = numeric(1))) |>
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
        scale_x_continuous(transform = "log2", breaks = breaks_log(10, 2))
    if (!is.null(facet_by)) {
        p <- p + facet_wrap(~ facet)
    }
    if (show_null) {
        p <- p +
            geom_ribbon(data = moran_mean_vars, aes(side, ymin = th1, ymax = th2),
                        fill = "darkslategray1", color = "cyan", alpha = 0.5) +
            geom_line(data = moran_mean_vars, aes(side, mean), color = "blue")
    }
    p + labs(title = "Moran's I across scales", x = "Bin size (μm)",
             y = "Moran's I", color = color_name)
}

#' Plot gene expression from multiple SFE objects
#'
#' This is useful to show what this gene looks like across different bin sizes.
#' When there are more than 200,000 cells or bins and no bbox is specified, then
#' scattermore is used to speed up plotting.
#'
#' @param sfes A list of SFE objects.
#' @param feature Gene to plot, only one gene
#' @param bbox Bounding box to plot a subarea
#' @param ... Other arguments to pass to
#'   \code{\link[Voyager]{plotSpatialFeature}}
#' @importFrom Voyager plotSpatialFeature
#' @importFrom ggplot2 ggtitle
#' @importFrom patchwork wrap_plots
#' @return A \code{patchwork} object
#' @export
plotSFEs <- function(sfes, feature, bbox = NULL, show_sizes = TRUE,
                     ncol = NULL, nrow = NULL, widths = NULL,
                     heights = NULL, design = NULL, ...) {
    args <- list(...)
    swap_rownames <- args[["swap_rownames"]]
    ps <- lapply(seq_along(sfes), function(i) {
        if (!is.null(swap_rownames))
            feature1 <- rownames(sfes[[i]])[rowData(sfes[[i]])[[swap_rownames]] == feature]
        else feature1 <- feature
        p <- if (ncol(sfes[[i]]) > 2e5 && is.null(bbox))
            plotSpatialFeature(sfes[[i]], feature, colGeometryName = "centroids",
                               scattermore = TRUE, ...)
        else plotSpatialFeature(sfes[[i]], feature, bbox = bbox, ...)
        if (show_sizes) {
            tt <- paste0(names(sfes)[[i]], " μm, I = ",
                         format(rowData(sfes[[i]])[feature1, "moran_sample01"],
                                digits = 3))
        } else {
            tt <- paste0(names(sfes)[[i]], ", I = ",
                         format(rowData(sfes[[i]])[feature1, "moran_sample01"],
                                digits = 3))
        }
        p + ggtitle(tt)
    })
    wrap_plots(ps, ncol = ncol, nrow = nrow, widths = widths,
               heights = heights, design = design)
}

#' Plot Lee's L curves
#'
#' The curves show how Lee's L changes over spatial scale for each pair of
#' genes, as indicated by bin size used to aggregate the data.
#'
#' @inheritParams plotMoranCurves
#' @param df Data frame output from \code{\link{clusterLeeCurves}}, or any data
#'   frame with a column named "pair" for gene pairs, "side" for bin size as in
#'   side length, and "lee" for Lee's L values, and optionally other categorical
#'   columns indicating cluster membership of genes.
#' @param show_median Logical, whether to plot the median of clusters, only used
#'   when \code{facet_by} is specified.
#' @param sample_n A smaller number of gene pairs to plot. Often the number of
#'   gene pairs is large so downsampling makes the plot more readable.
#' @return A \code{ggplot2} object.
#' @export
plotLeeCurves <- function(df, facet_by = NULL, show_median = FALSE,
                          sample_n = NULL) {
    if (!is.null(facet_by)) {
        df <- df |>
            mutate(facet = .data[[facet_by]])
        if (show_median) {
            df_med <- df |>
                group_by(facet, side) |>
                summarize(median = median(lee))
        }
    }
    if (!is.null(sample_n)) {
        df <- df |>
            filter(pair %in% sample(unique(df$pair), sample_n))
    }
    alpha <- if (is.null(sample_n)) 0.05 else 0.1
    p <- ggplot(df) +
        geom_line(aes(side, lee, group = pair), alpha = alpha) +
        scale_x_continuous(transform = "log2", breaks = breaks_log(10, 2))
    if (!is.null(facet_by)) {
        p <- p +
            facet_wrap(~ facet)
        if (show_median) {
            p <- p +
                geom_line(data = df_med, aes(side, median), color = "magenta", linewidth = 1)
        }
    }
    p + labs(title = "Lee's L across scales", x = "Bin size (μm)", y = "Lee's L")
}

#' Plot cluster median curves
#'
#' To better juxtapose different cluster patterns in the same plot.
#'
#' @param df Data frame output from either \code{\link{clusterLeeCurves}} or
#'   \code{\link{clusterMoranCurves}}. There must be a column named either
#'   "moran" or "lee" and a column called "side" in addition to the column
#'   specified in \code{cluster_col}.
#' @param cluster_col Name of a categorical column in \code{df} to use as
#'   cluster labels.
#' @return A \code{ggplot2} object.
#' @export
plotClusterMedians <- function(df, cluster_col) {
    col_use <- names(df)[names(df) %in% c("moran", "lee")]
    df_med <- df |>
        group_by(.data[[cluster_col]], side) |>
        summarize(median = median(.data[[col_use]]))
    ggplot(df_med, aes(side, median, color = .data[[cluster_col]])) +
        geom_line() +
        scale_x_continuous(transform = "log2", breaks = breaks_log(10, 2)) +
        scale_color_manual(values = ditto_colors) +
        labs(title = "Lee's L across scales", x = "Bin size (μm)", y = "Lee's L",
             color = "Cluster")
}

#' Plot pairs of genes in space from multiple SFE objects
#'
#' This function is useful to visualize pairs of genes whose Lee's L is of
#' interest.
#'
#' @inheritParams plotSFEs
#' @param sfes A list of SFE objects, whose names must be the bin sizes.
#' @param pair_use Gene pair to plot; the two gene names are separated by "_" as
#' in \code{\link{clusterLeeCurves}} output.
#' @param inds Indices of which SFE objects to plot
#' @param df_lee Data frame from \code{\link{clusterLeeCurves}}, where Lee's L
#' values are extracted to annotate the plot.
#' @return A \code{patchwork} object.
#' @importFrom stringr str_split_fixed
#' @export
plotSFEsLee <- function(sfes, pair_use, inds, df_lee, bbox = NULL,
                        show_sizes = TRUE, annotations = NULL, ...) {
    pair1 <- str_split_fixed(pair_use, "_", n = 2)[1,]
    df <- df_lee |> filter(pair == pair_use)
    if ("sample" %in% names(df_lee)) {
        l <- df$lee[match(names(sfes), df$sample)]
    } else l <- df$lee
    ps <- lapply(seq_along(inds), function(j) {
        i <- inds[[j]]
        p <- if (ncol(sfes[[i]]) > 2e5 && is.null(bbox))
            plotSpatialFeature(sfes[[i]], pair1, colGeometryName = "centroids",
                               scattermore = TRUE, ...)
        else plotSpatialFeature(sfes[[i]], pair1, bbox = bbox, ...)
        if (show_sizes)
            title <- paste0(names(sfes)[[i]], " μm, L = ", format(l[i], digits = 3))
        else
            title <- paste0(names(sfes)[[i]], ", L = ", format(l[i], digits = 3))
        if (!is.null(annotations)) title <- paste(title, annotations[j], sep = ", ")
        p[[1]] <- p[[1]] + ggtitle(title)
        p
    })
    wrap_plots(ps, tag_level = "keep")
}
