# Download CellChat database

Taken from \[CellChat's source
code\](https://github.com/jinworks/CellChat/tree/main/data), to use
without installing CellChat, which can't be a imported or suggested by
this package because CellChat is not on CRAN or Bioconductor.

## Usage

``` r
getCellChatDB(species = c("human", "mouse"), bfc = BiocFileCache())
```

## Arguments

- species:

  Either human or mouse

- bfc:

  [`BiocFileCache`](https://rdrr.io/pkg/BiocFileCache/man/BiocFileCache-class.html)
  instance where you can set where the files will be stored.

## Value

A list holding the database
