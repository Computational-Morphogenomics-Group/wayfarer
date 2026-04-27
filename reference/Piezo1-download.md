# Download Xenium data from Xue et al.

This dataset is from \[The mechanotransducer Piezo1 coordinates
metabolism and inflammation to promote skin
growth\](https://www.nature.com/articles/s41467-025-62270-3) by Xue et
al. To make the transcript spot files small enough to download as an
example dataset, I only kept transcripts assigned to cells and only kept
the columns x_location, y_location, cell_id, and feature_name.

## Usage

``` r
Piezo1TxSpots(
  sample = c("all", "expanded1", "expanded2", "expanded1_pi7", "expanded2_pi7",
    "nonexpanded1", "nonexpanded2", "yoda1", "yoda2"),
  bfc = BiocFileCache()
)

Piezo1TissueBoundary(
  sample = c("all", "expanded1", "expanded2", "expanded1_pi7", "expanded2_pi7",
    "nonexpanded1", "nonexpanded2", "yoda1", "yoda2"),
  bfc = BiocFileCache()
)

Piezo1Binned(
  sample = c("all", "expanded1", "expanded2", "expanded1_pi7", "expanded2_pi7",
    "nonexpanded1", "nonexpanded2", "yoda1", "yoda2"),
  bfc = BiocFileCache()
)

Piezo1BinAnalyses(
  sample = c("all", "expanded1", "expanded2", "expanded1_pi7", "expanded2_pi7",
    "nonexpanded1", "nonexpanded2", "yoda1", "yoda2"),
  bfc = BiocFileCache()
)
```

## Arguments

- sample:

  Which sample(s) to download, can have length greater than 1 to
  download multiple samples. If "all", then all samples will be
  downloaded. See [`sample_info`](sample_info.md) for information about
  each sample.

- bfc:

  [`BiocFileCache`](https://rdrr.io/pkg/BiocFileCache/man/BiocFileCache-class.html)
  instance where you can set where the files will be stored.

## Value

A character vector of file paths, except for `Piezo1TissueBoundary`
which returns a \`sf\` data frame.

## Details

These are the functions to download different aspects of the data:

- Piezo1TxSpots:

  Transcript spot coordinates in csv.gz files

- Piezo1TissueBoundary:

  Tissue boundary of each sample in RDS files; here the files will be
  read into R as a \`sf\` data frame.

- Piezo1Binned:

  Original binned data without processing or analyses in tar.gz

- Piezo1BinAnalyses:

  Binned data after removing edge bins, and with PCA, Moran's I, and
  Lee's L results
