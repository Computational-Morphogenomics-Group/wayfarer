.make_dl_function <- function(type) {
    function(sample = c("all", "expanded1", "expanded2",
                        "expanded1_pi7", "expanded2_pi7",
                        "nonexpanded1", "nonexpanded2",
                        "yoda1", "yoda2"),
             bfc = BiocFileCache()) {
        sample <- .validate_sample(sample)
        vapply(sample, .dl_sample1, bfc = bfc, type = type,
               FUN.VALUE = character(1))
    }
}

.validate_sample <- function(sample = c("all", "expanded1", "expanded2",
                                        "expanded1_pi7", "expanded2_pi7",
                                        "nonexpanded1", "nonexpanded2",
                                        "yoda1", "yoda2")) {
    all_samples <- c("expanded1", "expanded2",
                     "expanded1_pi7", "expanded2_pi7",
                     "nonexpanded1", "nonexpanded2",
                     "yoda1", "yoda2")
    sample <- match.arg(sample, several.ok = TRUE)
    if ("all" %in% sample) sample <- all_samples
    sample
}
.dl_sample1 <- function(sample, bfc, type = c("tx_spots", "tissue_boundary",
                                              "binned", "binned_esda")) {
    # Download 1 sample
    type <- match.arg(type)
    url <- switch(type,
                  tx_spots = "https://osf.io/b95ca",
                  tissue_boundary = "https://osf.io/4mxa8",
                  binned = "https://osf.io/kmhdt",
                  binned_esda = "https://osf.io/sk2nq")
    suff <- switch(type,
                   tx_spots = ".csv.gz",
                   tissue_boundary = "_tb.rds",
                   binned = ".tar.gz",
                   binned_esda = "_bin_analyses.tar.gz")
    fnm <- paste0(sample, suff)
    q <- bfcquery(bfc, fnm)
    n <- nrow(q)
    i <- 1
    if (n > 1) {
        message("multiple 'id' hits; using last")
        i <- n
    } else if (n > 0) {
        return(q$rpath[i])
    }
    no <- osf_retrieve_node(url)
    df <- osf_ls_files(no, pattern = fnm, type = "file")
    df <- df[str_detect(df$name, paste0("^", fnm)),]
    dir.create(td <- tempfile())
    osf_download(df, td, recurse=TRUE)
    bfcadd(bfc, fnm, fpath=file.path(td, fnm))
}

#' Download Xenium data from Xue et al.
#'
#' This dataset is from [The mechanotransducer Piezo1 coordinates metabolism and
#' inflammation to promote skin
#' growth](https://www.nature.com/articles/s41467-025-62270-3) by Xue et al. To
#' make the transcript spot files small enough to download as an example
#' dataset, I only kept transcripts assigned to cells and only kept the columns
#' x_location, y_location, cell_id, and feature_name.
#'
#' These are the functions to download different aspects of the  data:
#' \describe{
#' \item{Piezo1TxSpots}{Transcript spot coordinates in csv.gz files}
#' \item{Piezo1TissueBoundary}{Tissue boundary of each sample in RDS files;
#' here the files will be read into R as a `sf` data frame.}
#' \item{Piezo1Binned}{Original binned data without processing or analyses in tar.gz}
#' \item{Piezo1BinAnalyses}{Binned data after removing edge bins, and with PCA,
#' Moran's I, and Lee's L results}
#' }
#' @param sample Which sample(s) to download, can have length greater than 1 to
#'   download multiple samples. If "all", then all samples will be downloaded.
#'   See \code{\link{sample_info}} for information about each sample.
#' @param bfc \code{\link[BiocFileCache]{BiocFileCache}} instance where you can
#'   set where the files will be stored.
#' @importFrom BiocFileCache BiocFileCache bfcquery bfcadd
#' @importFrom osfr osf_retrieve_node osf_ls_files osf_download
#' @return A character vector of file paths, except for
#'   \code{Piezo1TissueBoundary} which returns a `sf` data frame.
#' @name Piezo1-download
#' @importFrom sf st_sf
NULL

#' @rdname Piezo1-download
#' @export
Piezo1TxSpots <- .make_dl_function("tx_spots")

#' @rdname Piezo1-download
#' @export
Piezo1TissueBoundary <- function(sample = c("all", "expanded1", "expanded2",
                                            "expanded1_pi7", "expanded2_pi7",
                                            "nonexpanded1", "nonexpanded2",
                                            "yoda1", "yoda2"),
                                 bfc = BiocFileCache()) {
    sample <- .validate_sample(sample)
    f <- .make_dl_function("tissue_boundary")
    fps <- f(sample, bfc)
    tbs <- lapply(fps, readRDS) |> unlist(recursive = FALSE) |> st_as_sfc()
    st_sf(sample = sample, geometry = tbs)
}

#' @rdname Piezo1-download
#' @export
Piezo1Binned <- .make_dl_function("binned")

#' @rdname Piezo1-download
#' @export
Piezo1BinAnalyses <- .make_dl_function("binned_esda")

#' Sample info for the Piezo1 Xenium data
#'
#' This is a data frame with two columns: sample and condition. Condition is the
#' experimental condition. PI14 stands for 14 days post complete injection, and
#' PI7 is 7 days. Injection means saline was injected into an implant beneath
#' the skin to expand the skin.
#'
#' @format A data frame
#' @source https://www.nature.com/articles/s41467-025-62270-3
"sample_info"
