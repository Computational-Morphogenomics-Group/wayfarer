# Filter Lee's L results for gene pairs in CellChatDB

This function takes CellChat's database without using the CellChat
package itself. It filter the Lee's L LMM results (see
[`runBinLMM`](runBinLMM.md)) to only include gene pairs in the database
for known interactions, and adds info about the interactions from the
database such as the KEGG pathway ID. Note that gene pairs that are not
in the database can still be interesting because they tell something
about cell neighborhoods.

## Usage

``` r
getCellChatInfo(
  df_lees_lmm,
  genes,
  species = c("human", "mouse"),
  bfc = BiocFileCache()
)
```

## Arguments

- df_lees_lmm:

  Output for Lee's L from [`runBinLMM`](runBinLMM.md)

- genes:

  A vector of gene symbols for the genes to consider; for smFISH-based
  technology, the list of gene of interest is often much smaller than
  the number of genes in the database.

- species:

  Either "human" or "mouse".

- bfc:

  A
  [`BiocFileCache`](https://rdrr.io/pkg/BiocFileCache/man/BiocFileCache-class.html)
  instance for where the CellChat database will be downloaded.

## Value

A data frame with the LMM results and additional columns from CellChat
database for more info about the interaction.
