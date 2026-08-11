.make_sims <- function(template, means, side, nsim, BPPARAM = SerialParam(), queen = FALSE) {
    areas <- side^2 * template$overlap_props
    simss <- lapply(means, function(m) {
        m_use <- m * areas
        sims <- lapply(seq_len(nsim), function(x) rpois(length(m_use), m_use))
        sims <- do.call(cbind, sims)
        as(sims, "CsparseMatrix")
    })
    sims <- do.call(cbind, simss)
    sims <- Matrix::t(sims)
    sfe <- SpatialFeatureExperiment(assays = list(counts = sims),
                                    colGeometries = list(bins = colGeometry(template)),
                                    sample_id = sampleIDs(template))
    rowData(sfe)$quantile <- rep(names(means), each = nsim)
    sfe <- sfe[,Matrix::colSums(counts(sfe)) > 0]
    sfe <- normalizeRnaCounts.se(sfe, center = FALSE)
    rownames(sfe) <- seq_len(nrow(sfe))
    cat("Running Moran's I\n")
    colGraph(sfe, "poly2nb") <- findSpatialNeighbors(sfe, type = "bins",
                                                     method = "poly2nb", queen = queen)
    runMoransI(sfe, NAOK = TRUE, BPPARAM = BPPARAM)
}

#' Create Monte Carlo simulations with no spatial structure
#'
#' The purpose of this function is to check if edge effect is causing spurious
#' changes in Moran's I and Lee's L across scales.This function uses geometries
#' from existing aggregations to create Monte Carlo simulations of gene
#' expression with no spatial structure. Poisson samples with mean informed by
#' the total counts of a given gene are taken on small bins, which are then
#' aggregated. Because given complete spatial randomness (CSR), the transcript
#' spots should be a homogeneous Poisson point process within the tissue
#' boundary. The number of counts within any given area follows a Poisson
#' distribution whose mean is the mean of the Poisson point process multiplied
#' by the size of that area. Therefore, to make computation faster, we can take
#' Poisson samples with means scaled by the area of bins in tissue.
#'
#' This function creates uncorrelated non-spatial patterns. They are stored as
#' SFE objects and written to disk. Then spatial metrics can be computed on
#' those patterns. This can be used to check if edge effect is causing artifacts
#' in the multi-scale analysis, and to compute Monte Carlo p-values for the
#' linear mixed models used to compare across biological conditions where the
#' log likelihood ratio test p-values are unreliable due to correlation of
#' spatial metrics between adjacent bin sizes and heteroscedasticity, i.e.
#' variance is higher for larger bins because there are fewer bins.
#'
#' Because edge effect artifacts affect genes with different expression levels
#' differently, \code{nsim} CSR patterns are generated for each quantile of
#' expression levels.
#'
#' @param in_path Path with SFE objects written to subdirectories with
#'   \code{\link{runBinAnalyses}}.
#' @param out_path Output directory
#' @param quantiles Vector of quantiles of total counts of genes to use for the
#'   intensity of the Poisson point process, because edge effect affects genes
#'   with different expression levels differently.
#' @param tissue_boundary \code{sf} or \code{sfc} for tissue boundary. See
#'   \code{\link[SpatialFeatureExperiment]{getTissueBoundaryConcave}} and
#'   \code{\link[SpatialFeatureExperiment]{getTissueBoundaryImg}} on obtaining
#'   the tissue boundary.
#' @return \code{out_path}; the output is written to disk there.
#' @importFrom SummarizedExperiment rowData<-
#' @export
makeMonteCarlo <- function(in_path, out_path, tissue_boundary,
                           quantiles = c(0.25, 0.5, 0.75),
                           nsim = 100, BPPARAM = SerialParam(), queen = FALSE) {
    out_path <- normalizePath(out_path, mustWork = FALSE)
    in_path <- normalizePath(in_path, mustWork = TRUE)
    if (!dir.exists(out_path)) dir.create(out_path)
    dirs <- list.files(in_path, pattern = "^bin\\d+_esda", full.names = TRUE)
    if (!length(dirs)) stop("No output found")
    sides <- str_extract(basename(dirs), "\\d+") |> as.integer()

    template <- readObject(dirs[1])
    tot_counts <- Matrix::rowSums(counts(template))
    qs <- quantile(tot_counts, quantiles)
    ms <- qs / st_area(tissue_boundary)

    lees_out <- list()
    morans_out <- list()
    for (i in seq_along(dirs)) {
        cat("Reading", basename(dirs[i]), "\n")
        if (i > 1) template <- readObject(dirs[i])
        sfe <- .make_sims(template, ms, sides[i], nsim, BPPARAM, queen)
        nm <- paste("moran", sampleIDs(sfe), sep = "_")
        morans_out[[as.character(sides[i])]] <- data.frame(moran = rowData(sfe)[[nm]],
                                                 gene = rownames(sfe),
                                                 side = sides[i] |> as.integer())
        # To do: use a subset of the features in each band, enough for nsims
        cat("Running Lee's L\n")
        lees_out[[as.character(sides[i])]] <- calculateBivariate(sfe, "lee", feature1 = rownames(sfe))
        cat("Saving simulation SFE\n")
        saveObject(sfe, file.path(out_path, paste(basename(dirs[i]), "sim", sep = "_")))
    }

    # Re-format all the Moran's I results into a data frame
    df_moran <- bind_rows(morans_out)
    df_moran$sample <- sampleIDs(sfe)
    df_lee <- .get_df_lee(lees_out)
    df_lee$sample <- sampleIDs(sfe)
    write.csv(df_moran, file.path(out_path, "df_moran.csv"), quote = FALSE,
              row.names = FALSE)
    write.csv(df_lee, file.path(out_path, "df_lee.csv"), quote = FALSE,
              row.names = FALSE)
    invisible(out_path)
}
