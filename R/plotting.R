ditto_colors <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00",
                  "#CC79A7", "#666666", "#AD7700", "#1C91D4", "#007756", "#D5C711",
                  "#005685", "#A04700", "#B14380", "#4D4D4D", "#FFBE2D", "#80C7EF",
                  "#00F6B3", "#F4EB71", "#06A5FF", "#FF8320", "#D99BBD", "#8C8C8C",
                  "#FFCB57", "#9AD2F2", "#2CFFC6", "#F6EF8E", "#38B7FF", "#FF9B4D",
                  "#E0AFCA", "#A3A3A3", "#8A5F00", "#1674A9", "#005F45", "#AA9F0D",
                  "#00446B", "#803800", "#8D3666", "#3D3D3D")

#' Get the mean and variance of Moran's I
#'
#' This function computes the mean and variance of Moran's I under the null
#' hypothesis of no spatial autocorrelation, in the context of plotting Moran's
#' I curves over bin sizes in \code{\link{plotMoranCurves}}.
#'
#' @param name Whether mean and variance are to be computed for Moran's I (moran)
#' or Geary's C (geary).
#' @param sfes A list of SFE objects of different bin sizes from the same sample
#'   (can be read with \code{\link{readBins}}), required to compute and plot the
#'   mean and variance of Moran's I (or Geary's C). The list must have names
#'   corresponding to bin sizes.
#' @return A data frame with columns side, mean, var, th1, and th2. Column th1
#'   is the 2.5% percentile in the normal distribution given mean and var, and
#'   th2 is the 97.5% percentile, after Bonferroni correction for the number of
#'   genes and number of bin sizes so the interval is rather conservative.
#' @export
#' @seealso [spdep::moran.test()]
getMoranMeanVar <- function(sfes, name = c("moran", "geary")) {
    name <- match.arg(name)
    test_fun <- switch (name,
                        moran = moran.test,
                        geary = geary.test
    )
    mean_vars <- tibble(side = as.numeric(names(sfes)),
                        mvs = lapply(sfes, function(x) {
                            foo <- test_fun(counts(x)[1,],
                                            listw = colGraph(x))
                            foo$estimate[2:3]
                        }),
                        mean = vapply(mvs, function(x) x[1],
                                      FUN.VALUE = numeric(1)),
                        var = vapply(mvs, function(x) x[2],
                                     FUN.VALUE = numeric(1))) |>
        dplyr::select(-mvs)
    mean_vars |>
        mutate(th1 = qnorm(0.025/nrow(sfes[[1]])/length(sfes),
                           mean = mean, sd = sqrt(var),
                           lower.tail = TRUE),
               th2 = qnorm(0.025/nrow(sfes[[1]])/length(sfes),
                           mean = mean, sd = sqrt(var),
                           lower.tail = FALSE))
}

#' Plot Moran's I curves
#'
#' The curves show how Moran's I changes over spatial scale, as indicated by bin
#' size used to aggregate the data. Despite the name, this function can be
#' applied to other univariate spatial results such as Geary's C.
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
#' @param color_name Name to show for colors in the legend.
#' @param name Column name in \code{df} with the spatial results.
#' @param mean_vars Data frame for means and variances of Moran's I under null
#'   hypothesis from \code{\link{getMoranMeanVar}}, required for \code{show_null
#'   = TRUE}.
#' @param title_name Name of the metric to put on the plot title.
#' @importFrom ggplot2 ggplot aes geom_line scale_x_continuous ggtitle
#'   scale_color_viridis_c facet_wrap geom_ribbon labs scale_color_manual
#' @importFrom rlang .data
#' @importFrom spdep moran.test geary.test
#' @importFrom SpatialFeatureExperiment colGraph
#' @importFrom SingleCellExperiment logcounts counts
#' @importFrom scales breaks_log
#' @importFrom stats qnorm
#' @importFrom ggplot2 ggplot geom_line aes scale_color_viridis_c
#'   scale_color_manual scale_x_continuous facet_wrap geom_ribbon labs
#' @return A \code{ggplot} object.
#' @export
plotMoranCurves <- function(df, name = "moran", color_by = NULL, facet_by = NULL,
                            show_null = FALSE, mean_vars = NULL,
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
        if (is.null(mean_vars)) {
            stop("mean_vars must be supplied when show_null = TRUE")
        }
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
    p + labs(title = paste(title_name, "across scales"), x = sprintf("Bin size (\u03BCm)"),
             y = title_name, color = color_name)
}

#' Plot gene expression from multiple SFE objects
#'
#' This is useful to show what this gene looks like across different bin sizes.
#' When there are more than 200,000 cells or bins and no bbox is specified, then
#' scattermore is used to speed up plotting.
#'
#' @inheritParams patchwork::wrap_plots
#' @param sfes A list of SFE objects.
#' @param feature Gene to plot, only one gene
#' @param bbox Bounding box to plot a subarea
#' @param title Title of the entire multi-panel plot
#' @param show_sizes Logical, whether to show bin sizes in the titles of
#'   individual panels. If TRUE, then \code{sfes} must have names that are the
#'   bin sizes
#' @param crop Logical, whether to only plot a cropped area when \code{bbox} is
#'   specified. If FALSE, then the bounding box will be shown as a box on the
#'   plot while the entire tissue section will be shown. This is used to show
#'   where a bounding box is to provide context along side a plot cropping by
#'   the bounding box.
#' @param ... Other arguments to pass to
#'   \code{\link[Voyager]{plotSpatialFeature}}
#' @importFrom Voyager plotSpatialFeature
#' @importFrom ggplot2 ggtitle geom_sf
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
        moran_name <- paste("moran", sampleIDs(sfes[[i]]), sep = "_")
        if (show_sizes) {
            tt <- paste0(names(sfes)[[i]], sprintf(" \u03BCm, I = "),
                         format(rowData(sfes[[i]])[feature1, moran_name],
                                digits = 3))
        } else {
            tt <- paste0(sampleIDs(sfes[[i]]), ", I = ",
                         format(rowData(sfes[[i]])[feature1, moran_name],
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
    p + labs(title = "Lee's L across scales", x = sprintf("Bin size (\u03BCm)"), y = "Lee's L")
}

#' Plot select Lee's L curves
#'
#' Plot Lee's L curves from select gene pairs, coloring by biological conditions
#'
#' @param df Data frame with Lee's L results from all samples, with columns
#'   pair, side, lee, sample, and group (for biological conditions).
#' @param lmm_res Linear mixed model results data frame from
#'   \code{\link{runBinLMM}}
#' @param pairs_use Which pairs to plot
#' @param title Title of the plot
#' @return A \code{ggplot2} oboject
#' @importFrom ggplot2 geom_hline
#' @export
plotLeeSelect <- function(df, lmm_res, pairs_use, title = NULL) {
    df <- df |>
        filter(pair %in% pairs_use) |>
        left_join(lmm_res, by = c("pair" = "feature"))
    if ("type" %in% names(lmm_res))
        df <- df |>
            mutate(label = paste0(feature, " (", type, ")"))
    else df <- df |> dplyr::filter(p_random_adj < 0.05)
    p <- ggplot(df, aes(side, lee, group = sample, color = group)) +
        geom_line() +
        geom_hline(color = "gray", yintercept = 0, linetype = 2) +
        scale_x_continuous(transform = "log2", breaks = scales::breaks_log(n = 10, base = 2)) +
        scale_color_manual(values = ditto_colors) +
        labs(x = sprintf("Bin size (\u03BCm)"), y = "Lee's L", title = title)
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
        labs(x = sprintf("Bin size (\u03BCm)"), y = "Median")
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
#' @param df_lee Data frame from \code{\link{clusterLeeCurves}}, where Lee's L
#'   values are extracted to annotate the plot.
#' @param exprs_value Which assay whose data should be plotted
#' @param side_use Which side length whose Lee's L is to be plotted. It must be
#' specified when plotting the same bin size from multiple samples.
#' @param bbox Named numeric vector specifying a bounding box, either to zoom
#'   into a smaller area or to show the box itself. The names should be "xmin",
#'   "xmax", "ymin", and "ymax" in any order. The same bbox will be used for all
#'   SFE objects in \code{sfes}.
#' @param sep Separator between the two gene symbols in \code{df_lee}
#' @param legend_index Index of position to place the legend. If \code{NULL},
#'   then the legend will not be shown. Unlike in univariate plots, the legend
#'   itself is a \code{ggplot2} object
#' @param crop Logical, whether to crop the samples by \code{bbox}. If
#'   \code{FALSE}, then the \code{bbox} will be shown as a box on the plot to
#'   indicate its location.
#' @return A \code{patchwork} object.
#' @importFrom stringr str_split_fixed
#' @importFrom Voyager plotBivariate
#' @importFrom patchwork plot_spacer wrap_plots wrap_elements
#' @importFrom biscale bi_legend
#' @export
plotSFEsBiscale <- function(sfes, feature1, feature2, df_lee, colGeometryName = 1L,
                            bbox = NULL, palette = "BlueGold", exprs_value = "logcounts",
                            swap_rownames = NULL, title = NULL, crop = TRUE,
                            show_sizes = TRUE, sep = "_", ncol = NULL, nrow = NULL,
                            side_use = NULL, widths = NULL,
                            heights = NULL, design = NULL, legend_index = NULL) {
    # Need to show Lee's L values and side lengths
    pair_use <- paste(sort(c(feature1, feature2)), collapse = sep)
    if (!show_sizes) {
        if (is.null(side_use) && "side" %in% names(df_lee)) stop("side_use must be specified")
    }
    plts <- lapply(seq_along(sfes), function(i) {
        if (show_sizes) {
            l <- df_lee |>
                filter(pair == pair_use, side == as.integer(names(sfes)[i])) |>
                pull(lee)
        } else {
            samples_use <- lapply(sfes, sampleIDs)
            l <- df_lee |>
                filter(pair == pair_use, sample == samples_use[i])
            if (!is.null(side_use)) {
                l <- l |> filter(side == side_use)
            }
            l <- l$lee
        }
        if (show_sizes)
            title_use <- paste0(names(sfes)[i], sprintf(" \u03BCm, L = "), format(l, digits = 3))
        else title_use <- paste0(sampleIDs(sfes[[i]]), ", L = ", format(l, digits = 3))
        if (!is.null(title) && i == 1L) {
            title_use <- paste(title, title_use, sep = "\n")
        }
        p <- if (crop || is.null(bbox)) {
            if (is.matrix(bbox)) bbox_use <- bbox[,i] else bbox_use <- bbox
            plotBivariate(sfes[[i]], feature1, feature2, bbox = bbox_use, swap_rownames = swap_rownames,
                          palette = palette, exprs_value = exprs_value, colGeometryName = colGeometryName)
        } else {
            if (is.matrix(bbox)) bbox_use <- bbox[,i] else bbox_use <- bbox
            plotBivariate(sfes[[i]], feature1, feature2, bbox = NULL, swap_rownames = swap_rownames,
                          palette = palette, exprs_value = exprs_value, colGeometryName = colGeometryName) +

                geom_sf(data = st_as_sfc(st_bbox(bbox_use)), linewidth = 0.5, fill = NA)
        }
        p + ggtitle(title_use) +
            theme(plot.title = element_text(vjust = 0, hjust = 0.5))
    })
    if (!is.null(legend_index)) {
        legend <- bi_legend(pal = palette,
                            dim = 4,
                            xlab = feature1,
                            ylab = feature2,
                            size = 14) +
            plot_spacer()
        # In case legend_index is 1 or length(plts)+1, i.e. legend is 1st or last panel
        if (isTRUE(all.equal(legend_index, 1L)))
            plts <- c(wrap_elements(legend), plts)
        else if (legend_index > length(plts))
            plts <- c(plts, wrap_elements(legend))
        else
            plts <- c(plts[1:(legend_index-1)], wrap_elements(legend), plts[legend_index:length(plts)])
    }
    wrap_plots(plts, ncol = ncol, nrow = nrow, widths = widths, heights = heights, design = design)
}
