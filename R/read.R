#' Read aggregated data from multiple bin sizes for the same sample
#'
#' This function reads from the output of \code{\link{runBinAnalyses}}. SFE
#' objects of the selected bin sizes will be loaded.
#'
#' @param dir Directory where the results are stored. The SFE object for each
#'   bin size must be in a directory whose name begins with "binx_esda" where x
#'   is the bin size, such as "bin12_esda".
#' @param sides Numeric vector of bin sizes to read
#' @return A list of SFE objects whose names are the bin sizes
#' @export
readBins <- function(dir, sides) {
    dir <- normalizePath(dir, mustWork = TRUE)
    sides <- sort(sides)
    fps <- file.path(dir, paste0("bin", sides, "_esda"))
    ex <- dir.exists(fps)
    if (all(!ex)) stop("None of the specified bin sizes found")
    if (any(!ex)) {
        warning("Bin size(s) ", paste(basename(fps[!ex]), collapse = ", "), " not found, skipping")
        fps <- fps[ex]
    }
    sfes <- lapply(fps, readObject)
    names(sfes) <- sides[ex]
    sfes
}

#' Read one bin size from multiple samples
#'
#' This function reads from the output of \code{\link{runBinAnalyses}}. Results
#' from one bin size are read from multiple samples.
#'
#' @param dirs A character vector of directories, one for each sample. The SFE
#'   object for each bin size must be in a directory whose name begins with
#'   "binx_esda" where x is the bin size, such as "bin12_esda".
#' @param side One bin size to read
#' @return A list of SFE objects whose names are base names in \code{dirs}.
#' @export
readBinSamples <- function(dirs, side) {
    dirs <- normalizePath(dirs, mustWork = FALSE)
    fps <- file.path(dirs, paste0("bin", side, "_esda"))
    ex <- dir.exists(fps)
    if (all(!ex)) stop("None of the specified bin sizes found")
    if (any(!ex)) {
        warning(paste0("bin", side, "_esda"), " not found in samples(s) ",
                paste(basename(dirs[!ex]), collapse = ", "), ", skipping")
        fps <- fps[ex]
    }
    sfes <- lapply(fps, readObject)
    names(sfes) <- basename(dirs[ex])
    sfes
}

#' Read multi-scale Moran's I results from multiple samples
#'
#' \code{\link{runBinAnalyses}} should write Moran's I for all genes in all bin
#' sizes to a file \code{df_moran.csv}. This function reads this file from
#' multiple samples and concatenates them.
#'
#' @param dirs A character vector of directories, one for each sample. There
#'   must be a \code{df_moran.csv} file in each directory.
#' @param sample_info Optional data frame with more info about each sample.
#'   There must be a column called "sample". Such info is helpful when plotting.
#' @return A data frame with columns moran, side, sample, and optionally other
#'   columns in \code{sample_info}.
#' @export
readMoranSamples <- function(dirs, sample_info = NULL) {
    dirs <- normalizePath(dirs, mustWork = FALSE)
    fps <- file.path(dirs, "df_moran.csv")
    ex <- file.exists(fps)
    if (all(!ex)) stop("df_moran.csv not found in any of the samples")
    if (any(!ex)) {
        warning("df_moran.csv not found in samples(s) ",
                paste(basename(dirs[!ex]), collapse = ", "), ", skipping")
        fps <- fps[ex]
    }
    if (!is.null(sample_info)) {
        stopifnot(is.data.frame(sample_info))
        stopifnow("sample" %in% names(sample_info))
    }
    dfs <- lapply(fps, read.csv, col.names = TRUE)
    out <- bind_rows(dfs)
    if (!is.null(sample_info)) {
        out <- out |> left_join(sample_info, by = "sample")
    }
    out
}

#' Read multi-sample Lee's L results from multiple samples
#'
#' \code{\link{runBinAnalyses}} should write Lee's L for all genes in all bin
#' sizes to a file \code{df_lee.csv}. This function reads this file from
#' multiple samples and concatenates them.
#'
#' @inheritParams readMoranSamples
#' @param dirs A character vector of directories, one for each sample. There
#'   must be a \code{df_lee.csv} file in each directory.
#' @param cutoff Gene pairs whose absolute values of Lee's L is below this
#'   cutoff in all bin sizes and all samples will not be included in the output.
#'   This will remove genes that are not spatially correlated at any length
#'   scale and any sample.
#' @return A data frame with columns moran, side, sample, and optionally other
#'   columns in \code{sample_info}.
#' @export
readLeeSamples <- function(dirs, sample_info = NULL, cutoff = 0L) {
    dirs <- normalizePath(dirs, mustWork = FALSE)
    fps <- file.path(dirs, "df_lee.csv")
    ex <- file.exists(fps)
    if (all(!ex)) stop("df_lee.csv not found in any of the samples")
    if (any(!ex)) {
        warning("df_lee.csv not found in samples(s) ",
                paste(basename(dirs[!ex]), collapse = ", "), ", skipping")
        fps <- fps[ex]
    }
    if (!is.null(sample_info)) {
        stopifnot(is.data.frame(sample_info))
        stopifnow("sample" %in% names(sample_info))
    }
    dfs <- lapply(fps, read.csv, col.names = TRUE)
    out <- bind_rows(dfs)
    if (cutoff > 0L) {
        pairs_use <- out |>
            group_by(pair, sample) |>
            summarize(max = max(abs(lee))) |>
            filter(max > cutoff) |> pull(pair) |> unique()
        out <- out |>
            filter(pair %in% pairs_use)
    }
    if (!is.null(sample_info)) {
        out <- out |> left_join(sample_info, by = "sample")
    }
    out
}
