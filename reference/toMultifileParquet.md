# Convert transcript spot files to multi-file Parquet

This greatly improves speed of spatial aggregation.

## Usage

``` r
toMultifileParquet(
  tx_file,
  tx_parquet_path,
  tech = c("other", "Vizgen", "Xenium", "CosMX"),
  spatialCoordsNames = c("X", "Y", "Z"),
  gene_col = "gene",
  phred_col = "qv",
  min_phred = 20
)
```

## Arguments

- tx_file:

  File with transcript spot coordinates. Can be a directory for
  multi-file parquet.

- tx_parquet_path:

  If the input is not a Parquet file (e.g. a CSV file), it will be
  re-formatted into a multi-file Parquet to improve speed of
  computation. The reformatted files can be written to a path specified
  in this argument. If it's `NULL`, then the multi-file Parquet will be
  written to a temporary directory.

- tech:

  Technology whose standard output the transcript file is from. If it's
  not "other", then arguments `spatialCoordsNames` and `gene_col` will
  be ignored as the column names from the standard output will be used
  instead.

- spatialCoordsNames:

  Column names for the x, y, and optionally z coordinates of the spots.
  The defaults are for Vizgen.

- gene_col:

  Column name for genes.

- phred_col:

  Column name for Phred scores of the spots.

- min_phred:

  Minimum Phred score to keep spot. By default 20, the conventional
  threshold indicating "acceptable", meaning that there's 1 chance that
  the spot was decoded in error.

## Value

Nothing into the R session; the reformatted files are written to disk,
to the directory specified in `tx_parquet_path`.
