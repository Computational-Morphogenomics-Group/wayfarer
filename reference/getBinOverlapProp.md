# Find portion of each bin overlapping the tissue or cells

If tissue boundary is used, then it should allow for holes, because
holes can also cause edge effect.

## Usage

``` r
getBinOverlapProp(
  sfe,
  tissue_geometry,
  BPPARAM = SerialParam(),
  batch_size = 1000,
  prop = TRUE
)
```

## Arguments

- sfe:

  SFE object with the aggregated bins, generated with
  [`makeAggregates`](makeAggregates.md).

- tissue_geometry:

  Either tissue boundary (with holes if present in tissue) or cell
  segmentation polygons. Area of overlap of each bin with this geometry
  will be computed.

- BPPARAM:

  A [`bpparam`](https://rdrr.io/pkg/BiocParallel/man/register.html), to
  parallelize over batches of bins.

- batch_size:

  Number of bins in each batch.

- prop:

  Logical, whether to return proportions of bin area in tissue instead
  of actual area.

## Value

A numeric vector same length as `ncol(sfe)`.
