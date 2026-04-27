# Aggregate transcript spots into bins of varying sizes

This function aggregates transcript spots from smFISH-based spatial
transcriptomics technologies into bins of various sizes. In the
resulting gene count matrix, the count of eah gene in each bin is the
number of transcript spots of that gene intersecting with each bin. QC
of the original dataset is strongly recommended before this spatial
aggregation.

## Usage

``` r
makeAggregates(
  tx_file,
  out_path,
  sides = sort(c(2^(3:8), 12 * 2^(0:5))),
  sample_id = "sample01",
  tech = c("other", "Vizgen", "Xenium", "CosMX"),
  spatialCoordsNames = c("X", "Y"),
  gene_col = "gene",
  phred_col = "qv",
  min_phred = 20,
  tissue_boundary = NULL,
  flip_geometry = TRUE,
  BPPARAM = SerialParam(),
  tx_parquet_path = NULL,
  ...
)
```

## Arguments

- tx_file:

  File with transcript spot coordinates. Can be a directory for
  multi-file parquet.

- out_path:

  Output directory

- sides:

  Side length (for square bins) or hexagon diameter in microns. This
  specifies the size of the bins used to aggregate transcripts.

- sample_id:

  Which sample in the SFE object the transcript spots should be added
  to.

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

- tissue_boundary:

  Optional but recommended. A `sf`, `sfc`, or `sfg` of the tissue
  boundary polygon. See
  [`getTissueBoundaryImg`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/getTissueBoundaryImg.html)
  and
  [`getTissueBoundaryConcave`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/getTissueBoundaryConcave.html)
  on getting the tissue boundary with SFE, after removing debris in QC.
  The tissue boundary will be used to remove bins that don't overlap
  with tissue. If QC is not performed before aggregation, then you'll
  have to remove debris in every single aggregation.

- flip_geometry:

  Logical, whether to flip the transcript spot geometries to match the
  images if added later.

- BPPARAM:

  bpparam object to specify parallel computing over genes. If a lot of
  memory is used, then stick to \`SerialParam()\`. Using multicore in
  the R level when \`save_memory = TRUE\` works and results into faster
  computation.

- tx_parquet_path:

  If the input is not a Parquet file (e.g. a CSV file), it will be
  re-formatted into a multi-file Parquet to improve speed of
  computation. The reformatted files can be written to a path specified
  in this argument. If it's `NULL`, then the multi-file Parquet will be
  written to a temporary directory.

- ...:

  Other arguments passed to
  [`aggregateTx`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/aggregateTx.html).
  This excludes `cellsize` and `save_memory`.

## Value

Invisibly the output directory; the output is written to disk with
`alabaster.sfe` as on disk serializations of `SpatialFeatureExperiment`
objects, with one object for each bin size. The SFE object for each bin
size will be written to a subdirectory of `out_path` named "binx" where
"x" is the bin size, e.g. "bin12" for 12 micron bins.

## Details

If the original file for the transcript spot coordinates is not in
parquet format, then it will be converted to a multi-file parquet where
each file is for a gene, which will greatly improve the speed of
aggregation.
