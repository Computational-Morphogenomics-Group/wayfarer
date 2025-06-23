#' Plot variance explained across bin sizes
#'
#' This function plots the iconic PCA elbow plot for variance explained or
#' eigenvalue for multiple SCE objects, distinguished by color.
#'
#' @param sfes A list of SFE (SCE and SPE are fine) objects with the dimension.
#' Names of the list are bin sizes.
#' reduction of interest computed and stored in \code{reducedDims}.
#' @param reduction Name of the dimension reduction of interest
#' @param field Field in the attribute of the dimension reduction where variance
#' explained or eigenvalues are stored.
#' @param npcs Number of components whose variance explained or eigenvalues are
#' to be plotted.
#' @param xlab X-axis label
#' @param ylab Y-axis label
#' @return A \code{ggplot} object
#' @importFrom SingleCellExperiment reducedDim
#' @export
plotEigenvalues <- function(sfes, npcs = 20, reduction = "PCA", field = "varExplained",
                            xlab = "Principal component", ylab = "Variance explained") {
    df_ves <- tibble(side = as.integer(names(sfes)),
                     x = log2(side),
                     data = lapply(sfes, function(x) {
                         y <- attr(reducedDim(x, reduction), field)
                         tibble(eig = y,
                                PC = seq_along(y))
                     })) |>
        unnest(cols = data)
    df_ves |>
        filter(PC <= npcs) |>
        ggplot() +
        geom_line(aes(PC, eig, group = side, color = side)) +
        scale_color_continuous(transform = "log2") +
        scale_x_continuous(breaks = scales::breaks_width(2)) +
        labs(x = "Principal component", y = "Variance explained",
             color = "Bin size (μm)")
}

.calculate_angle <- function(v1, v2) {
    # Credit: Joe Rich
    dot_product <- sum(v1 * v2)

    magnitude_vector1 <- sqrt(sum(v1^2))
    magnitude_vector2 <- sqrt(sum(v2^2))

    cos_theta <- dot_product / (magnitude_vector1 * magnitude_vector2)

    if (isTRUE(all.equal(cos_theta, 1))) {
        angle_radians <- 0
    } else if (isTRUE(all.equal(cos_theta, -1))) {
        angle_radians <- pi
    } else {
        angle_radians <- acos(cos_theta)
    }

    return (angle_radians)
}

.flip_pc <- function(v1, v2) {
    theta1_2 <- .calculate_angle(v1, v2)
    theta1_neg2 <- .calculate_angle(v1, -v2)

    if (theta1_neg2 < theta1_2) {
        v2 <- -v2
    }
    v2
}

.flip_pc_sfe <- function(sfes, pc, reduction, field) {
    genes_use <- rownames(sfes[[1]])
    for (i in seq(2, length(sfes), by = 1)) {
        v1 <- attr(reducedDim(sfes[[i-1L]], reduction), field)
        if (is.null(v1)) break
        v1 <- v1[genes_use, pc]
        v2 <- attr(reducedDim(sfes[[i]], reduction), field)[genes_use, pc]

        theta1_2 <- .calculate_angle(v1, v2)
        theta1_neg2 <- .calculate_angle(v1, -v2)

        if (theta1_neg2 < theta1_2) {
            reducedDim(sfes[[i]], reduction)[,pc] <- -reducedDim(sfes[[i]], reduction)[,pc]
            attr(reducedDim(sfes[[i]], reduction), field)[genes_use, pc] <- -v2
        }
    }
    sfes
}

.get_loadings <- function(sfes, pc = 1, reduction = "PCA", field = "rotation",
                          flip = FALSE) {
    # They should have the same genes
    genes_use <- rownames(sfes[[1]])
    pcs <- lapply(sfes, function(x) {
        m <- attr(reducedDim(x, reduction), field)
        tibble(gene = genes_use,
               loadings = as.vector(m[genes_use,pc]))
    })
    if (flip) pcs[[1]]$loadings <- -pcs[[1]]$loadings
    for (i in seq(2, length(pcs), by = 1)) {
        pcs[[i]] <- pcs[[i]] |>
            mutate(loadings = .flip_pc(pcs[[i-1L]]$loadings, loadings))
    }
    tibble(side = as.integer(names(sfes)),
           data = pcs,
           PC = pc) |>
        unnest(data)
}

#' Plot gene loadings of a given principal component across scales
#'
#' This can visualize how the meanings of principal components change with
#' spatial scales.
#'
#' @param sfes A list of SFE objects with names as bin sizes
#' @param pc Index of principal component or eigenvector or loading vector to
#'   plot. Facetting is used when length of this argument is greater than 1.
#' @param reduction Name of dimension reduction to plot
#' @param field Name of the attribute of the dimension reduction that holds the
#'   gene loadings
#' @param color_by Whether to color the curves by Moran's I or loading values.
#' Can only color by Moran's I when multiple PCs are plotted because each PC has
#' different gene loadings.
#' @param color_ind Since one curve can only have one color, choose the index of
#'   the SFE object in \code{sfes} whose values are used for coloring.
#' @param flip Logical, whether to flip the signs of the loadings.
#' @param ... Arguments passed to \code{\link[ggplot2]{facet_wrap}}.
#' @importFrom scico scale_color_scico
#' @importClassesFrom dplyr bind_rows
#' @return A \code{ggplot} object
#' @export
plotLoadingCurves <- function(sfes, pc = 1, reduction = "PCA", field = "rotation",
                              color_by = c("loading", "moran"),
                              color_ind = 1, flip = FALSE, ...) {
    color_by <- match.arg(color_by)
    df_pc <- lapply(pc, .get_loadings, sfes = sfes, reduction = reduction,
                    field = field, flip = flip) |>
        bind_rows()
    sides <- as.integer(names(sfes))
    if (length(pc) > 1L) color_by <- "moran"
    if (color_by == "moran") {
        df_morans <- rowData(sfes[[color_ind]]) |> as.data.frame() |>
            rownames_to_column("gene") |>
            select(gene, moran = starts_with("moran"))
        df_pc <- df_pc |>
            left_join(df_morans, by = "gene")
        p <- ggplot(df_pc, aes(side, loadings, group = gene, color = moran)) +
            scale_color_viridis_c(name = "Moran's I", option = "E", end = 0.95)
    } else {
        df_col <- df_pc |>
            filter(side == sides[color_ind]) |>
            dplyr::select(gene, color = loadings)
        df_pc <- df_pc |>
            left_join(df_col, by = "gene")
        p <- ggplot(df_pc, aes(side, loadings, group = gene, color = color)) +
            scale_color_scico(palette = "roma", direction = -1, midpoint = 0,
                              name = "Loading")
    }
    if (length(pc) > 1L) {
        p <- p + facet_wrap(~ PC, ...)
        title <- NULL
    } else title <- paste0("Gene loadings for PC ", pc)
    p +
        geom_line(alpha = 0.3) +
        scale_x_continuous(transform = "log2", breaks = scales::breaks_log(n = 10, base = 2)) +
        geom_hline(yintercept = 0, color = "gray", linetype = 2) +
        labs(x = "Bin size (μm)", y = "Gene loading",
             title = title)
}

#' Plot Moran's I of bin projection in PCA space
#'
#' Basically see how spatially structured each principal component is and how it
#' may relate to bin size.
#'
#' @inheritParams plotEigenvalues
#' @param sfes A list of SFE objects, with Moran's I computed for bin projections
#' in the dimension reduction of interest.
#' @return A \code{ggplot} object
#' @export
plotMoranPCs <- function(sfes, npcs = 20, reduction = "PCA") {
    df_moran_pcs <- tibble(side = sides,
           x = log2(side),
           data = lapply(sfes, function(x) {
               # Too lazy to change now; it assumes that Moran's I is the first one computed
               y <- reducedDimFeatureData(x, reduction)[seq_len(npcs),1]
               tibble(moran = y,
                      PC = seq_along(y))
           })) |>
        unnest(cols = data) |>
        filter(PC <= npcs)
    df_moran_pcs |>
        ggplot(aes(PC, moran, group = side, color = side)) +
        geom_line() +
        scale_color_continuous(trans = "log2") +
        scale_x_continuous(breaks = scales::breaks_width(2)) +
        labs(x = "Principal component", y = "Moran's I",
             color = "Bin size (μm)")
}

#' Plot dimension reduction in space for multiple SFE objects
#'
#' This function not only plots dimension reductions in space, but also aligns
#' the components for a consistent look because in PCA, the PC's can be flipped
#' in sign.
#'
#' @inheritParams plotLoadingCurves
#' @param sfes A list of SFE objects, for which the dimension reduction of
#'   interest has been computed.
#' @param field Field in the attribute of the dimension reduction where variance
#' explained or eigenvalues are stored. If this field is not found, then the
#' components will not be flipped and flipping is not necessary for dimension
#' reductions such as NMF which only produce non-negative results.
#' @return A \code{patchwork} object
#' @importFrom Voyager spatialReducedDim
#' @importFrom patchwork plot_annotation
#' @importFrom SingleCellExperiment reducedDim<-
#' @export
plotSFEsRedDim <- function(sfes, pc = 1, reduction = "PCA", field = "rotation",
                           divergent = TRUE, diverge_center = 0, bbox = NULL) {
    sfes <- .flip_pc_sfe(sfes, pc = pc, reduction = reduction, field = field)
    plts <- lapply(seq_along(sfes), function(i) {
        p <- if (ncol(sfes[[i]]) > 2e5 && is.null(bbox))
            spatialReducedDim(sfes[[i]], dimred = reduction, components = pc,
                              colGeometryName = "centroids", scattermore = TRUE,
                              divergent = divergent, diverge_center = diverge_center)
        else spatialReducedDim(sfes[[i]], dimred = reduction, components = pc,
                               divergent = divergent, diverge_center = diverge_center,
                               bbox = bbox)
        p + ggtitle(paste0(names(sfes)[[i]], " μm"))
    })
    wrap_plots(plts) + plot_annotation(title = paste0(reduction, " component ", pc))
}
