#' Aggregate transcript spots into bins of varying sizes
#'
#' This function aggregates transcript spots from smFISH-based spatial
#' transcriptomics technologies into bins of various sizes. In the resulting
#' gene count matrix, the count of eah gene in each bin is the number of
#' transcript spots of that gene intersecting with each bin. QC of the original
#' dataset is strongly recommended before this spatial aggregation.
#'
#' If the original file for the transcript spot coordinates is not in parquet
#' format, then it will be converted to a multi-file parquet where each file is
#' for a gene, which will greatly improve the speed of aggregation.
#'
#' @inheritParams SpatialFeatureExperiment::aggregateTx
#' @param sides Side length (for square bins) or hexagon diameter in microns.
#'   This specifies the size of the bins used to aggregate transcripts.
#' @param tx_file File with transcript spot coordinates. Can be a directory for
#'   multi-file parquet.
#' @param out_path Output directory
#' @param tech Technology whose standard output the transcript file is from. If
#'   it's not "other", then arguments \code{spatialCoordsNames} and
#'   \code{gene_col} will be ignored as the column names from the standard
#'   output will be used instead.
#' @param tissue_boundary Optional but recommended. A \code{sf}, \code{sfc}, or
#'   \code{sfg} of the tissue boundary polygon. See
#'   \code{\link[SpatialFeatureExperiment]{getTissueBoundaryImg}} and
#'   \code{\link[SpatialFeatureExperiment]{getTissueBoundaryConcave}} on getting
#'   the tissue boundary with SFE, after removing debris in QC. The tissue
#'   boundary will be used to remove bins that don't overlap with tissue. If QC
#'   is not performed before aggregation, then you'll have to remove debris in
#'   every single aggregation.
#' @param ... Other arguments passed to
#'   \code{\link[SpatialFeatureExperiment]{aggregateTx}}. This excludes
#'   \code{cellsize} and \code{save_memory}.
#' @param tx_parquet_path If the input is not a Parquet file (e.g. a CSV file),
#'   it will be re-formatted into a multi-file Parquet to improve speed of
#'   computation. The reformatted files can be written to a path specified in
#'   this argument. If it's \code{NULL}, then the multi-file Parquet will be
#'   written to a temporary directory.
#' @return Invisibly the output directory; the output is written to disk with
#'   \code{alabaster.sfe} as on disk serializations of
#'   \code{SpatialFeatureExperiment} objects, with one object for each bin size.
#'   The SFE object for each bin size will be written to a subdirectory of
#'   \code{out_path} named "binx" where "x" is the bin size, e.g. "bin12" for 12
#'   micron bins.
#' @export
#' @importFrom SpatialFeatureExperiment .check_tx_file getTechTxFields
#'   aggregateTx crop colGeometry
#' @importFrom alabaster.sfe saveObject
#' @importFrom arrow open_dataset write_dataset
#' @importFrom rlang .data
#' @importFrom dplyr select group_by filter
#' @importFrom tools file_ext
#' @importFrom zeallot %<-%
#' @importFrom BiocParallel SerialParam bplapply
makeAggregates <- function(tx_file, out_path,
                           sides = sort(c(2^(3:8), 12*2^(0:5))),
                           sample_id = "sample01",
                           tech = c("other", "Vizgen", "Xenium", "CosMX"),
                           spatialCoordsNames = c("X", "Y"),
                           gene_col = "gene", phred_col = "qv", min_phred = 20,
                           tissue_boundary = NULL, flip_geometry = TRUE,
                           BPPARAM = SerialParam(), tx_parquet_path = NULL, ...) {
    tx_file <- normalizePath(tx_file, mustWork = TRUE)
    out_path <- normalizePath(out_path, mustWork = FALSE)
    tech <- match.arg(tech)
    if (!dir.exists(out_path)) dir.create(out_path)
    if (!is.null(tissue_boundary) && !any(class(tissue_boundary) %in% c("sf", "sfc", "sfg"))) {
        stop("tissue_boundary must be sf, sfc, or sfg")
    }
    if (tech != "other") {
        c(spatialCoordsNames, gene_col, cell_col, fn) %<-%
            getTechTxFields(tech, NULL)
    }
    is_parquet <- grepl("\\.parquet$", tx_file) |
        (dir.exists(tx_file) &
             all(grepl("\\.parquet$", list.files(tx_file, recursive = TRUE,
                                                 include.dirs = FALSE))))

    # Convert to multi-file parquet
    if (!is_parquet) {
        if (is.null(tx_parquet_path))
            fp <- tempfile(pattern = "tx")
        else {
            fp <- normalizePath(tx_parquet_path, mustWork = FALSE)
            if (!dir.exists(fp)) dir.create(fp)
        }
        toMultifileParquet(tx_file, tx_parquet_path = fp, tech = tech,
                           spatialCoordsNames = spatialCoordsNames,
                           gene_col = gene_col, phred_col = phred_col,
                           min_phred = min_phred)
        spatialCoordsNames <- c("x", "y")
        gene_col <-  "gene"
    } else fp <- tx_file

    for (s in sides) {
        cat("Processing", s, "\n")
        sfe <- aggregateTx(fp, cellsize = s, sample_id = sample_id,
                           spatialCoordsNames = spatialCoordsNames[1:2],
                           flip_geometry = flip_geometry, gene_col = gene_col,
                           save_memory = TRUE, BPPARAM = BPPARAM, ...)
        if (!is.null(tissue_boundary))
            sfe <- crop(sfe, tissue_boundary, keep_whole = "col")
        saveObject(sfe, file.path(out_path, paste0("bin", s)))
        gc()
    }
    invisible(out_path)
}

#' Find portion of each bin overlapping the tissue or cells
#'
#' If tissue boundary is used, then it should allow for holes, because holes can
#' also cause edge effect.
#'
#' @param sfe SFE object with the aggregated bins, generated with
#'   \code{\link{makeAggregates}}.
#' @param tissue_geometry Either tissue boundary (with holes if present in
#'   tissue) or cell segmentation polygons. Area of overlap of each bin with
#'   this geometry will be computed.
#' @param BPPARAM A \code{\link[BiocParallel]{bpparam}}, to parallelize over
#'   batches of bins.
#' @param batch_size Number of bins in each batch.
#' @param prop Logical, whether to return proportions of bin area in tissue
#'   instead of actual area.
#' @return A numeric vector same length as \code{ncol(sfe)}.
#' @importFrom sf st_union st_intersects st_covered_by st_area st_intersection
#'   st_overlaps st_as_sfc st_bbox
#' @export
getBinOverlapProp <- function(sfe, tissue_geometry, BPPARAM = SerialParam(),
                              batch_size = 1000, prop = TRUE) {
    bins <- colGeometry(sfe, "bins")
    bins$index <- seq_len(ncol(sfe))
    bins$batch <- ceiling(bins$index/batch_size)
    AREA <- st_area(st_geometry(bins)[1]) # should be the same for all bins
    areas <- bplapply(unique(bins$batch), function(i) {
        g <- st_geometry(bins)[bins$batch == i]
        out <- numeric(length(g))
        xc <- st_union(st_geometry(tissue_geometry)[st_intersects(st_as_sfc(st_bbox(g)), tissue_geometry, sparse = FALSE)])
        out[st_covered_by(g, xc, sparse = FALSE) |> as.vector()] <- AREA
        inds_comp <- st_overlaps(g, xc, sparse = FALSE) |> as.vector()
        out[inds_comp] <- st_area(st_intersection(g, xc))
        out
    }, BPPARAM = BPPARAM) |> unlist()
    areas[is.na(areas)] <- 0
    if (prop) areas <- areas/AREA
    areas
}

#' Remove bins that overlap too little with the tissue boundary
#'
#' This will help mitigate edge effect when analyzing the binned data,
#' especially for larger bins. Bins with the lowest percentiles in area
#' overlapping tissue will be removed. However, this will not perfectly remove
#' edge effect in downstream analyses; it will only mitigate edge effect.
#'
#' @inheritParams getBinOverlapProp
#' @param overlap_props A numeric vector with proportion of area of each bin in
#'   tissue/cells, or a single string for a column in \code{colData(sfe)} with
#'   such proportions.
#' @param min_prop Minimum proportion of bin area in tissue; bins with lower
#'   proportions will be removed. The default of 0.9 is suitable for small bins
#'   with bin side less than 24 microns; a smaller number should be used for
#'   larger bins, such as 0.8 for 32 microns, 0.5 for 64 microns, and 0.2 for
#'   192 microns and above.
#' @param quantile Bins with overlapping proportion below this quantile will be
#'   removed. If not NULL, it will supercede \code{min_prop}. You can check the
#'   histogram of \code{overlap_props} before deciding \code{min_prop} and
#'   \code{quantile}.
#' @return A SFE object with the bins that overlap too little with tissue
#'   removed.
#' @export
#' @importFrom SummarizedExperiment colData
removeEdgeBins <- function(sfe, overlap_props, min_prop = 0.9, quantile = NULL) {
    if (is.numeric(overlap_props) && length(overlap_props) != ncol(sfe)) {
        stop("Length of numeric overlap_props must be the same as ncol(sfe)")
    }
    if (!is.null(quantile)) {
        stopifnot(is.numeric(quantile))
        stopifnot(length(quantile) == 1L)
        stopifnot(quantile > 0 & quantile < 1)
    }
    if (is.character(overlap_props)) {
        if (length(overlap_props) > 1L) stop("Character overlap_props must have length 1")
        if (!overlap_props %in% names(colData(sfe))) {
            stop(overlap_props, " is not in colData(sfe)")
        }
        overlap_props <- colData(sfe)[,overlap_props]
    }
    if (is.null(quantile)) {
        sfe[,overlap_props > min_prop]
    } else {
        sfe[,overlap_props > quantile(overlap_props, probs = quantile)]
    }
}
