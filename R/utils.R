#' Convert transcript spot files to multi-file Parquet
#'
#' This greatly improves speed of spatial aggregation.
#'
#' @inheritParams makeAggregates
#' @return Nothing into the R session; the reformatted files are written to
#'   disk, to the directory specified in \code{tx_parquet_path}.
#' @export
toMultifileParquet <- function(tx_file, tx_parquet_path,
                               tech = c("other", "Vizgen", "Xenium", "CosMX"),
                               spatialCoordsNames = c("X", "Y", "Z"),
                               gene_col = "gene", phred_col = "qv", min_phred = 20) {
    tx_file <- normalizePath(tx_file, mustWork = TRUE)
    tx_parquet_path <- normalizePath(tx_parquet_path, mustWork = FALSE)
    tech <- match.arg(tech)
    if (tech != "other") {
        c(spatialCoordsNames, gene_col, cell_col, fn) %<-%
            getTechTxFields(tech, NULL)
    }
    orig_format <- file_ext(tx_file)
    ds <- open_dataset(tx_file, format = "csv")
    ds <- ds |>
        select(gene = .data[[gene_col]], x = .data[[spatialCoordsNames[1]]],
               y = .data[[spatialCoordsNames[2]]]) |>
        group_by(gene)
    if (phred_col %in% names(ds)) {
        ds <- df |>
            filter(.data[[phred_col]] > min_phred)
    }
    write_dataset(ds, path = tx_parquet_path)
}
