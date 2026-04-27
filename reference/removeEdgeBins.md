# Remove bins that overlap too little with the tissue boundary

This will help mitigate edge effect when analyzing the binned data,
especially for larger bins. Bins with the lowest percentiles in area
overlapping tissue will be removed. However, this will not perfectly
remove edge effect in downstream analyses; it will only mitigate edge
effect.

## Usage

``` r
removeEdgeBins(sfe, overlap_props, min_prop = 0.9, quantile = NULL)
```

## Arguments

- sfe:

  SFE object with the aggregated bins, generated with
  [`makeAggregates`](makeAggregates.md).

- overlap_props:

  A numeric vector with proportion of area of each bin in tissue/cells,
  or a single string for a column in `colData(sfe)` with such
  proportions.

- min_prop:

  Minimum proportion of bin area in tissue; bins with lower proportions
  will be removed. The default of 0.9 is suitable for small bins with
  bin side less than 24 microns; a smaller number should be used for
  larger bins, such as 0.8 for 32 microns, 0.5 for 64 microns, and 0.2
  for 192 microns and above.

- quantile:

  Bins with overlapping proportion below this quantile will be removed.
  If not NULL, it will supercede `min_prop`. You can check the histogram
  of `overlap_props` before deciding `min_prop` and `quantile`.

## Value

A SFE object with the bins that overlap too little with tissue removed.
