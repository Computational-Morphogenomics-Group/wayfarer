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
#' @importFrom spdep moran.test geary.test
#' @importFrom SpatialFeatureExperiment colGraph
#' @importFrom SingleCellExperiment logcounts
#' @importFrom scales breaks_log
#' @importFrom stats qnorm
#' @return A \code{ggplot} object.
#' @export
plotUnivariateCurves <- function(df, name = "moran", color_by = NULL, facet_by = NULL,
                            show_null = FALSE, sfes = NULL,
                            color_name = NULL, title_name = "Moran's I") {
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
            geom_line(aes(side, .data[[name]], group = gene, color = color), alpha = 0.5)
        if (is.numeric(df$color))
            p <- p + scale_color_viridis_c(option = "E")
        else
            p <- p + scale_color_manual(values = ditto_colors, na.value = "gray80")
    } else
        p <- ggplot(df) +
            geom_line(aes(side, .data[[name]], group = gene), alpha = 0.5)
    if (show_null) {
        if (is.null(sfes)) {
            stop("sfes must be supplied when show_null = TRUE")
        }
        test_fun <- switch (name,
            moran = moran.test,
            geary = geary.test
        )
        mean_vars <- tibble(side = as.integer(names(sfes)),
                                  mvs = lapply(sfes, function(x) {
                                      foo <- test_fun(logcounts(x)[1,],
                                                      listw = colGraph(x))
                                      foo$estimate[2:3]
                                  }),
                                  mean = vapply(mvs, function(x) x[1],
                                                FUN.VALUE = numeric(1)),
                                  var = vapply(mvs, function(x) x[2],
                                               FUN.VALUE = numeric(1))) |>
            dplyr::select(-mvs)
        mean_vars <- mean_vars |>
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
            geom_ribbon(data = mean_vars, aes(side, ymin = th1, ymax = th2),
                        fill = "darkslategray1", color = "cyan", alpha = 0.5) +
            geom_line(data = mean_vars, aes(side, mean), color = "blue")
    }
    p + labs(title = paste(title_name, "across scales"), x = "Bin size (μm)",
             y = title_name, color = color_name)
}

#' @rdname plotUnivariateCurves
#' @export
plotMoranCurves <- function(df, color_by = NULL, facet_by = NULL,
                            show_null = FALSE, sfes = NULL,
                            color_name = NULL) {
    plotUnivariateCurves(df, name = "moran", color_by = color_by,
                         facet_by = facet_by, show_null = show_null,
                         sfes = sfes, color_name = color_name,
                         title_name = "Moran's I")
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
plotSFEs <- function(sfes, feature, bbox = NULL, title = NULL,
                     show_sizes = TRUE,
                     ncol = NULL, nrow = NULL, widths = NULL,
                     heights = NULL, design = NULL, crop = TRUE, ...) {
    args <- list(...)
    swap_rownames <- args[["swap_rownames"]]
    ps <- lapply(seq_along(sfes), function(i) {
        if (!is.null(swap_rownames))
            feature1 <- rownames(sfes[[i]])[rowData(sfes[[i]])[[swap_rownames]] == feature]
        else feature1 <- feature
        p <- if (ncol(sfes[[i]]) > 2e5 && is.null(bbox))
            plotSpatialFeature(sfes[[i]], feature, colGeometryName = "centroids",
                               scattermore = TRUE, ...)
        else if (is.matrix(bbox)) {
            # Each column of the bbox matrix for each sfe in the list
            if (crop) {
                plotSpatialFeature(sfes[[i]], feature, bbox = bbox[,i], ...)
            } else {
                plotSpatialFeature(sfes[[i]], feature, bbox = NULL, ...) +
                    geom_sf(data = st_as_sfc(st_bbox(bbox[,i])), linewidth = 1, fill = NA)
            }
        } else {
            if (crop || is.null(bbox)) {
                plotSpatialFeature(sfes[[i]], feature, bbox = bbox, ...)
            } else {
                plotSpatialFeature(sfes[[i]], feature, bbox = NULL, ...) +
                    geom_sf(data = st_as_sfc(st_bbox(bbox)), linewidth = 1, fill = NA)
            }
        }
        if (show_sizes) {
            tt <- paste0(names(sfes)[[i]], " μm, I = ",
                         format(rowData(sfes[[i]])[feature1, "moran_sample01"],
                                digits = 3))
        } else {
            tt <- paste0(names(sfes)[[i]], ", I = ",
                         format(rowData(sfes[[i]])[feature1, "moran_sample01"],
                                digits = 3))
        }
        if (!is.null(title) && i == 1L) {
            tt <- paste0(title, "\n", tt)
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

#' Plot select Lee's L curves
#'
#' Also annotate with relevant linear mixed model info
#'
#' @param df_lees Data frame with Lee's L results
#' @param lmm_res Linear mixed model results data frame
#' @param pairs_use Which pairs to plot
#' @param title Title of the plot
#' @return A \code{ggplot2} oboject
#' @export
plotLeeSelect <- function(df_lees, lmm_res, pairs_use, title = NULL) {
    df <- df_lees |>
        filter(pair %in% pairs_use) |>
        left_join(lmm_res, by = "pair")
    if ("type" %in% names(lmm_res))
        df <- df |>
            mutate(label = paste0(pair, " (", type, ")"))
    else df <- df |> dplyr::filter(p_random_adj < 0.05)
    p <- ggplot(df, aes(side, value, group = sample, color = stage)) +
        geom_line() +
        geom_hline(color = "gray", yintercept = 0, linetype = 2) +
        scale_x_continuous(transform = "log2", breaks = scales::breaks_log(n = 10, base = 2)) +
        scale_color_manual(values = ditto_colors[c(27,2,7,1)]) +
        labs(x = "Bin size (μm)", y = "Lee's L", color = "Stage", title = title)
    if ("type" %in% names(lmm_res))
        p <- p + facet_wrap(~ label)
    else p <- p + facet_wrap(~ pair)
    p
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
#' @inheritParams Voyager::plotBivariate
#' @inheritParams patchwork::plot_layout
#' @param sfes A list of SFE objects, whose names must be the bin sizes.
#' @param pair_use Gene pair to plot
#' @param swap_colors Logical, whether to swap the positions of the two genes
#' in `pair_use` in the palette.
#' @param df_lee Data frame from \code{\link{clusterLeeCurves}}, where Lee's L
#' values are extracted to annotate the plot.
#' @return A \code{patchwork} object.
#' @importFrom stringr str_split_fixed
#' @export
plotSFEsBiscale <- function(sfes, pair_use, df_lee, colGeometryName = 1L,
                            bbox = NULL, palette = "BlueGold", exprs_value = "logcounts",
                            swap_rownames = NULL, title = NULL, crop = TRUE, swap_colors = FALSE,
                            show_sizes = TRUE, sep = "(?<!HLA)-", ncol = NULL, nrow = NULL,
                            samples_use = NULL, side_use = NULL, widths = NULL,
                            heights = NULL, design = NULL, legend_index = NULL) {
    # Need to show Lee's L values and side lengths
    genes <- str_split(pair_use, pattern = sep, simplify = TRUE) |> as.vector()
    if (swap_colors) genes <- rev(genes)
    if (!show_sizes) {
        if (length(samples_use) != length(sfes)) stop("samples_use must have the same length as sfes")
        if (is.null(side_use) && "side" %in% names(df_lee)) stop("side_use must be specified")
    }
    plts <- lapply(seq_along(sfes), function(i) {
        if (show_sizes) {
            l <- df_lee |>
                filter(pair == pair_use, side == as.integer(names(sfes)[i])) |>
                pull(lee)
        } else {
            l <- df_lee |>
                filter(pair == pair_use, sample == samples_use[i])
            if (!is.null(side_use)) {
                l <- l |> filter(side == side_use)
            }
            l <- l$lee
        }
        if (show_sizes)
            title_use <- paste0(names(sfes)[i], " μm, L = ", format(l, digits = 3))
        else title_use <- paste0(names(sfes)[i], ", L = ", format(l, digits = 3))
        if (!is.null(title) && i == 1L) {
            title_use <- paste(title, title_use, sep = "\n")
        }
        p <- if (crop || is.null(bbox)) {
            if (is.matrix(bbox)) bbox_use <- bbox[,i] else bbox_use <- bbox
            .gene_biscale(sfes[[i]], genes[1], genes[2], bbox = bbox_use, swap_rownames = swap_rownames,
                          palette = palette, exprs_value = exprs_value, colGeometryName = colGeometryName)
        } else {
            if (is.matrix(bbox)) bbox_use <- bbox[,i] else bbox_use <- bbox
            .gene_biscale(sfes[[i]], genes[1], genes[2], bbox = NULL, swap_rownames = swap_rownames,
                          palette = palette, exprs_value = exprs_value, colGeometryName = colGeometryName) +

                geom_sf(data = st_as_sfc(st_bbox(bbox_use)), linewidth = 0.5, fill = NA)
        }
        p + ggtitle(title_use) +
            theme(plot.title = element_text(vjust = 0, hjust = 0.5))
    })
    if (!is.null(legend_index)) {
        legend <- bi_legend(pal = "BlueGold",
                            dim = 4,
                            xlab = genes[1],
                            ylab = genes[2],
                            size = 14) +
            plot_spacer()
        plts <- c(plts[1:(legend_index-1)], wrap_elements(legend), plts[legend_index:length(plts)])
    }
    wrap_plots(plts, ncol = ncol, nrow = nrow, widths = widths, heights = heights, design = design)
}
