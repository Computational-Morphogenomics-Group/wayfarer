#' Run basic analyses on aggregated SFE objects of all bin sizes
#'
#' This function runs log normalization, bin-level QC, PCA, adjacency graph,
#' Moran's I, and Lee's L for all binned SFE objects in a directory, given that
#' bins overlapping too little with tissue have been removed.
#'
#' @inheritParams Voyager::calculateUnivariate
#' @inheritParams spdep::poly2nb
#' @inheritParams getBinOverlapProp
#' @param dir Directory where the binned outputs are located. Output of each bin
#'   size must be in a directory named "binx", such as "bin12" for 12 micron
#'   bins.
#' @param ncomponents Number of components to compute for PCA.
#' @param min_props Minimum proportion of each bin overlapping tissue; bins that
#'   don't overlap enough will be removed to mitigate edge effect. This should
#'   be a numeric vector of the same length as the number of bin sizes and can
#'   differ for different bin sizes. If length 1, then the same value will be
#'   applied to all bin sizes. Otherwise, it must have names corresponding to
#'   bin sizes.
#' @param quantiles Numeric vector of quantiles of area of bin overlapping
#'   tissue, can differ for different bin sizes. If not NULL, this will
#'   supercede \code{min_props}.
#' @return Invisibly \code{out_dir}; the log normalization, QC, PCA, adjacency
#'   graph, and Moran's I will be stored in the SFE object; the SFE object with
#'   the results will be written to \code{out_dir} with directory names
#'   "binx_esda". A data frame for Moran's I and  Lee's L across bin sizes will
#'   be written to \code{out_dir} as CSV files.
#' @importFrom scater logNormCounts runPCA addPerCellQC
#' @importFrom Voyager runMoransI calculateBivariate
#' @importFrom alabaster.base readObject
#' @importFrom SpatialFeatureExperiment sampleIDs
#' @importFrom dplyr mutate bind_rows
#' @importFrom tidyr pivot_longer unite
runBinAnalyses <- function(dir, out_dir, tissue_geometry,
                           min_props = 0.9, quantiles = NULL,
                           ncomponents = 30, queen = FALSE,
                           zero.policy = TRUE, p.adjust.method = "BH",
                           BPPARAM = SerialParam(), ...) {
    dir <- normalizePath(dir, mustWork = TRUE)
    bins_dir <- list.files(dir, "^bin\\d+$", include.dirs = TRUE, full.names = TRUE)
    if (!length(bins_dir)) {
        stop("Binned output should be in directories named bin<x> such as bin12")
    }
    out_dir <- normalizePath(out_dir, mustWork = FALSE)
    if (!dir.exists(out_dir)) dir.create(out_dir)
    bin_sizes <- gsub("bin", "", basename(bins_dir))
    if (is.null(quantiles)) {
        if (!is.numeric(min_props) || (length(min_props)!=1L || length(min_props)!=length(bins_dir))) {
            stop("min_props must be a numeric vector of either length 1 or same length as the number of bin sizes")
        }
        if (length(min_props) > 1L) {
            if (!setequal(bin_sizes, names(min_props))) {
                stop("min_props must have the same names as bin sizes")
            }
        } else {
            min_props <- setNames(rep(min_props, length(bin_sizes)), bin_sizes)
        }
    } else {
        if (!is.numeric(quantiles) || (length(quantiles)!=1L || length(quantiles)!=length(bins_dir))) {
            stop("quantiles must be a numeric vector of either length 1 or same length as the number of bin sizes")
        }
        if (length(quantiles) > 1L) {
            if (!setequal(bin_sizes, names(quantiles))) {
                stop("quantiles must have the same names as bin sizes")
            }
        } else {
            quantiles <- setNames(rep(quantiles, length(bin_sizes)), bin_sizes)
        }
    }
    lees_out <- list()
    moran_out <- list()
    for (i in seq_along(bins_dir)) {
        d <- bins_dir[i]
        cat("Reading", basename(d), "\n")
        sfe <- readObject(d)
        cat("Removing edge bins\n")
        bs <- min(1000, round(ncol(sfe)/10))
        overlap_props <- getBinOverlapProp(sfe, tissue_geometry, BPPARAM = BPPARAM,
                                           batch_size = bs)
        sfe <- removeEdgeBins(sfe, overlap_props, min_prop = min_props[bin_sizes[i]],
                              quantile = quantiles[bin_sizes[i]])
        cat("Normalizing data\n")
        sfe <- addPerCellQC(sfe)
        sfe <- logNormCounts(sfe)
        cat("Running PCA\n")
        sfe <- runPCA(sfe, ncomponents = ncomponents, scale = TRUE)
        cat("Running Moran's I\n")
        sfe <- runMoransI(sfe, zero.policy = zero.policy, BPPARAM = BPPARAM)
        nm <- paste("moran", sampleIDs(sfe), sep = "_")
        morans_out[[bin_sizes[i]]] <- data.frame(moran = rowData(sfe)[[nm]],
                                                 side = bin_sizes[i] |> as.integer())
        cat("Saving SFE basis analysis\n")
        saveObject(sfe, file.path(out_dir, paste(basename(d), "esda", sep = "_")))
        cat("Running Lee's L\n")
        lees_out[[bin_sizes[i]]] <- calculateBivariate(sfe, "lee", feature1 = rownames(sfe))
    }
    # Re-format all the Moran's I results into a data frame
    df_moran <- bind_rows(morans_out)
    morans_out$sample <- sampleIDs(sfe)
    df_lee <- .get_df_lee(lees_out)
    df_lee$sample <- sampleIDs(sfe)
    write.csv(df_moran, file.path(out_dir, "df_moran.csv"), quote = FALSE, col.names = TRUE,
              row.names = FALSE)
    write.csv(df_lee, file.path(out_dir, "df_lee.csv"), quote = FALSE, col.names = TRUE,
              row.names = FALSE)
    invisible(out_dir)
}

.get_pairs_df <- function(lees) {
    nr <- nrow(lees[[1]])
    inds_df <- data.frame(j = unlist(lapply(seq_len(nr-1), function(x) rep(x+1, times = x))),
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

.get_df_lee <- function(lees, symbol_df = NULL) {
    lees <- lapply(lees, as.matrix)
    rns <- sort(rownames(lees[[1]]))
    lees <- lapply(lees, function(x) x[rns, rns])
    if (!is.null(symbol_df)) {
        lees <- lapply(lees, function(x) {
            rns2 <- symbol_df[rns, "Symbol"]
            rownames(x) <- colnames(x) <- rns2
            rns2 <- sort(rns2)
            x[rns2, rns2]
        })
    }
    df_lee <- .get_pairs_df(lees)
    df_lee |>
        pivot_longer(-pair, names_to = "name", values_to = "value")
}
