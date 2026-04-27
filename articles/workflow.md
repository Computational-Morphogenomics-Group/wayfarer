# workflow

## Introduction

In this vignette, we apply the multi-scale spatial analysis workflow as
shown in the [Wayfarer
paper](https://www.biorxiv.org/content/10.64898/2026.02.16.706245v1.abstract),
but on a different dataset. This Xenium dataset is from the paper [The
mechanotransducer Piezo1 coordinates metabolism and inflammation to
promote skin
growth](https://www.nature.com/articles/s41467-025-62270-3), which
investigates glycolysis and inflammation in mouse skin expansion,
mediated by Piezo1. Skin expansion was performed by injecting saline
into a subdermal implant. This is relevant because skin expansion occurs
in pregnancy and obesity, and is performed to create more skin for skin
grafting. This dataset has two samples (I believe biological replica)
per condition:

- Non-expanded – the implants were inserted but not expanded
- Expanded, 14 days post inflation (PI14)
- Expanded control (PI7)
- Expanded skin with topical Yoda1 treatment (PI7) which activates
  Piezo1

In this vignette, we perform a multi-scale spatial analysis with
Wayfarer to identify genes and gene co-expression whose multi-scale
spatial patterns differ across conditions.

``` r
library(Wayfarer)
library(Voyager)
library(SpatialFeatureExperiment)
library(sf)
library(BiocParallel)
library(alabaster.sfe)
library(dplyr)
library(readr)
library(ggplot2)
library(stringr)
library(purrr)
library(tidyr)
library(scater)
library(R.utils)
theme_set(theme_bw())
```

## Single sample

### Creating spatial bins

First, we apply spatial aggregation of various bin sizes to one sample
and perform some basic analyses. You can do the same analyses and
visualizations on all the other samples. The tissue boundary was
computed from a concave hull of cell centroids as the images were not
provided by the authors. See
[`SpatialFeatureExperiment::getTissueBoundaryConcave()`](https://pachterlab.github.io/SpatialFeatureExperiment/reference/getTissueBoundaryConcave.html).
For the purpose of demonstration, we can start with sample `expanded1`,
one of the expanded skin samples 14 days post injection.

The example dataset has been uploaded to OSF and can be downloaded with
functions in the Wayfarer package.

``` r
tb <- Piezo1TissueBoundary(sample = "expanded1")
plot(st_geometry(tb))
```

![](workflow_files/figure-html/unnamed-chunk-2-1.png)

Also download the transcript spots

``` r
(tx_path <- Piezo1TxSpots("expanded1"))
#>                                                          expanded1 
#> "/home/runner/.cache/R/BiocFileCache/e0ccda281f8_expanded1.csv.gz"
```

With the transcript spots, we can create spatial aggregates with a range
of bin sizes (microns)

``` r
(sides <- sort(c(2^(3:8), 12 * 2^(0:4))))
#>  [1]   8  12  16  24  32  48  64  96 128 192 256
```

The tissue boundary is used to remove bins that are entirely outside
tissue. Often there are transcript spots detected outside tissue. Also
unlike in the Wayfarer paper, hexagonal bins rather than square ones are
used here to show that hexagonal bins also work.

``` r
binned_path <- makeAggregates(tx_path,
                              out_path = "expanded1", sample_id = "expanded1",
                              tech = "Xenium", tissue_boundary = tb, sides = sides,
                              flip_geometry = FALSE,
                              square = FALSE, BPPARAM = MulticoreParam(2, progressbar = TRUE))
```

Since this takes a while to run, for the purpose of rendering this
vignette on GitHub Actions, the results have been uploaded to OSF and
can be downloaded

``` r
(binned_path <- Piezo1Binned("expanded1"))
#>                                                           expanded1 
#> "/home/runner/.cache/R/BiocFileCache/e0cc2108640f_expanded1.tar.gz"
```

Decompress it into the current directory

``` r
untar(binned_path)
```

``` r
min_props <- c(0.9, 0.9, 0.85, 0.85, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.3)
names(min_props) <- sides
```

Bins that overlap too little with tissue will be outliers in analyses
and cause artifacts, so here bins that overlap the tissue below a
proportion will be removed. These proportions were found by manual
tuning in the LUAD data used in the Wayfarer paper. Percentiles can also
be used with the `quantiles` argument. A value should be supplied to
each bin size.

The `runBinAnalyses` performs the basic analyses on each bin size after
removing bins that overlap too little with the tissue:

1.  Log normalize the data, using proportion of each bin overlapping
    tissue as size factor instead of total counts. Using the proportion
    vs. total counts makes a big difference in downstream analyses.
2.  Non-spatial PCA on log normalized data
3.  Moran’s I on each gene, with log normalized data
4.  Lee’s L on all gene pairs, with log normalized data

``` r
runBinAnalyses("expanded1", "expanded1/bin_analyses", tissue_geometry = tb,
               min_props = min_props,
               BPPARAM = MulticoreParam(2, progressbar = TRUE))
```

Since the analyses can take a while to run, the results can also be
downloaded from OSF

``` r
(bin_analyses_path <- Piezo1BinAnalyses("expanded1"))
#>                                                                        expanded1 
#> "/home/runner/.cache/R/BiocFileCache/e0cc1598bc79_expanded1_bin_analyses.tar.gz"
untar(bin_analyses_path)
```

### PCA

Load one of the bin sizes to take a look, here 16 micron bin diameter

``` r
(sfe <- readObject("expanded1/bin_analyses/bin16_esda/"))
#> >>> Reading SpatialExperiment
#> >>> Reading colgeometries
#> >>> Reading spatial graphs
#> Warning in sn2listw(df, style = style, zero.policy = zero.policy,
#> from_mat2listw = TRUE): neighbour object has 5 sub-graphs
#> Warning in mat2listw(as(altReadObject(dddd), "CsparseMatrix"), style =
#> method$args$style, : neighbour object has 5 sub-graphs
#> class: SpatialFeatureExperiment 
#> dim: 100 9748 
#> metadata(0):
#> assays(2): counts logcounts
#> rownames(100): Acaca Adgre1 ... Gsdmc Ngf
#> rowData names(2): moran_expanded1 K_expanded1
#> colnames(9748): 956 1113 ... 78406 78407
#> colData names(6): sample_id overlap_props ... total sizeFactor
#> reducedDimNames(1): PCA
#> mainExpName: NULL
#> altExpNames(0):
#> spatialCoords names(2) : X Y
#> imgData names(4): sample_id image_id data scaleFactor
#> 
#> unit: micron
#> Geometries:
#> colGeometries: bins (POLYGON) 
#> 
#> Graphs:
#> expanded1: col: poly2nb
```

``` r
ElbowPlot(sfe)
```

![](workflow_files/figure-html/unnamed-chunk-12-1.png)

``` r
plotReducedDim(sfe, "PCA", ncomponents = 4) + geom_density2d()
```

![](workflow_files/figure-html/unnamed-chunk-13-1.png)

That circle in PC1 and PC2 looks interesting

``` r
plotDimLoadings(sfe)
```

![](workflow_files/figure-html/unnamed-chunk-14-1.png)

Rotate the tissue to save space when plotting

``` r
sfe <- rotateMinRect(sfe)
```

``` r
plotSpatialFeature(sfe, "Col17a1")
```

![](workflow_files/figure-html/unnamed-chunk-16-1.png)

``` r
spatialReducedDim(sfe, "PCA", ncomponents = 4, ncol = 1, divergent = TRUE, diverge_center = 0)
```

![](workflow_files/figure-html/unnamed-chunk-17-1.png)

### Moran’s I curves

Next see how Moran’s I changes with bin sizes in this sample

``` r
df_moran <- read_csv("expanded1/bin_analyses/df_moran.csv")
#> Rows: 1100 Columns: 4
#> ── Column specification ────────────────────────────────────────────────────────
#> Delimiter: ","
#> chr (2): gene, sample
#> dbl (2): moran, side
#> 
#> ℹ Use `spec()` to retrieve the full column specification for this data.
#> ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```

``` r
plotMoranCurves(df_moran)
```

![](workflow_files/figure-html/unnamed-chunk-19-1.png)

Here each gene has a curve. We see different kinds of curves: some genes
have higher Moran’s I at small scales that decrease at larger scales.
Some form a peak between 16 and 32 microns in bin diameter. Some form a
peak between 64 and 128 microns. We can cluster the curves to better
visualize these different patterns:

``` r
df_moran2 <- clusterMoranCurves(df_moran)
names(df_moran2)
#> [1] "moran"        "gene"         "side"         "sample"       "hclust"      
#> [6] "leiden"       "hclust_diffs" "leiden_diffs"
```

Hierarchical clustering and Leiden are used to cluster the curves, both
the original values and differences between adjacent bin sizes
(hclust_diffs and leiden_diffs).

Plot the clusters

``` r
plotMoranCurves(df_moran2, facet_by = "leiden")
```

![](workflow_files/figure-html/unnamed-chunk-21-1.png)

Plot the cluster medians

``` r
plotClusterMedians(df_moran2, "leiden")
```

![](workflow_files/figure-html/unnamed-chunk-22-1.png)

Plot the Moran’s I curves with the mean and variance of Moran’s I under
the null hypothesis of no spatial autocorrelation, and the 2.5% and
97.5% quantiles after Bonferroni correction for the number of genes and
bin sizes

``` r
# Need to read all bin sizes for the spatial neighborhood graphs
sfes <- readBins("expanded1/bin_analyses", sides = sides)
#> >>> Reading SpatialExperiment
#> >>> Reading colgeometries
#> >>> Reading spatial graphs
#> >>> Reading SpatialExperiment
#> >>> Reading colgeometries
#> >>> Reading spatial graphs
#> >>> Reading SpatialExperiment
#> >>> Reading colgeometries
#> >>> Reading spatial graphs
#> Warning in sn2listw(df, style = style, zero.policy = zero.policy,
#> from_mat2listw = TRUE): neighbour object has 5 sub-graphs
#> Warning in mat2listw(as(altReadObject(dddd), "CsparseMatrix"), style =
#> method$args$style, : neighbour object has 5 sub-graphs
#> >>> Reading SpatialExperiment
#> >>> Reading colgeometries
#> >>> Reading spatial graphs
#> Warning in sn2listw(df, style = style, zero.policy = zero.policy,
#> from_mat2listw = TRUE): neighbour object has 2 sub-graphs
#> Warning in mat2listw(as(altReadObject(dddd), "CsparseMatrix"), style =
#> method$args$style, : neighbour object has 2 sub-graphs
#> >>> Reading SpatialExperiment
#> >>> Reading colgeometries
#> >>> Reading spatial graphs
#> >>> Reading SpatialExperiment
#> >>> Reading colgeometries
#> >>> Reading spatial graphs
#> >>> Reading SpatialExperiment
#> >>> Reading colgeometries
#> >>> Reading spatial graphs
#> >>> Reading SpatialExperiment
#> >>> Reading colgeometries
#> >>> Reading spatial graphs
#> >>> Reading SpatialExperiment
#> >>> Reading colgeometries
#> >>> Reading spatial graphs
#> >>> Reading SpatialExperiment
#> >>> Reading colgeometries
#> >>> Reading spatial graphs
#> >>> Reading SpatialExperiment
#> >>> Reading colgeometries
#> >>> Reading spatial graphs
mean_vars <- getMoranMeanVar(sfes, "moran")
```

``` r
plotMoranCurves(df_moran2, facet_by = "leiden", mean_vars = mean_vars,
                show_null = TRUE)
```

![](workflow_files/figure-html/unnamed-chunk-24-1.png)

### Lee’s L curves

Also read Lee’s L; to make clustering easier, a cutoff of 0.2 is set so
that gene pairs with Lee’s L lower than 0.2 in all bin sizes will be
removed.

``` r
df_lee <-readLeeSamples("expanded1/bin_analyses", cutoff = 0.2)
```

Due to the large number of gene pairs, a smaller subset can be plotted
to prevent the plot from becoming a solid block

``` r
plotLeeCurves(df_lee, sample_n = 1000)
```

![](workflow_files/figure-html/unnamed-chunk-26-1.png)

Here each gene pair has a curve. Some gene pairs have higher Lee’s L at
smaller scales which decrease at larger scales. Some peak between 16 and
32 microns. Some peak between 64 and 128 microns. Some have low Lee’s L
at small scales which increase at larger scales. Some have negative
Lee’s L at small scales which become positive at larger scales. Some
have weakly negative Lee’s L at smaller scales which becomes lower at
intermediate scales.

Cluster Lee’s L curves to better visualize these patterns

``` r
df_lee2 <- clusterLeeCurves(df_lee)
```

``` r
plotLeeCurves(df_lee2, facet_by = "hclust_diffs", show_median = TRUE) +
    geom_hline(yintercept = 0, linetype = 2, color = "gray")
```

![](workflow_files/figure-html/unnamed-chunk-28-1.png)

``` r
plotClusterMedians(df_lee2, "leiden")
```

![](workflow_files/figure-html/unnamed-chunk-29-1.png)

Unfotunately we can’t plot that interval under null hypothesis like in
Moran’s I because the mean and variance of Lee’s L under the null
hypothesis of no spatial autocorrelation differ for each gene pair,
depending on their Pearson correlation. When there’s no spatial pattern,
gene pairs with high Pearson correlation can still get a moderate value
away from 0. In addition, computing such mean and variance requires a
dense square matrix, so they can only be computed for larger bin sizes
where the number of bins is smaller.

### Plotting multiple bin sizes

Here we plot gene expression in space in different bin sizes

``` r
# Reduce empty space
for (i in seq_along(sfes)) {
    sfes[[i]] <- rotateMinRect(sfes[[i]])
}
```

Randomly select a gene from Leiden cluster 4

``` r
set.seed(29)
gene_use <- df_moran2 |> 
    filter(leiden == 4) |> 
    pull(gene) |> 
    sample(1)
```

``` r
plotSFEs(sfes[c("16", "48", "96", "256")], gene_use, show_sizes = TRUE, ncol = 1)
```

![](workflow_files/figure-html/unnamed-chunk-32-1.png)

We can also make this kind of plot for Lee’s L, with a bivariate
palette. Randomly select a gene pair

``` r
set.seed(29)
pair_use <- df_lee2 |> 
    filter(hclust_diffs == 4) |> 
    pull(pair) |> 
    sample(1) |> 
    str_split(pattern = "_", simplify = TRUE) |> 
    as.vector()
```

``` r
plotSFEsBiscale(sfes[c("16", "48", "96", "256")], pair_use[1], pair_use[2], 
                df_lee = df_lee, ncol = 1, legend_index = 1,
                heights = rep(1, 5))
```

![](workflow_files/figure-html/unnamed-chunk-34-1.png)

## Multiple samples

After seeing the Moran’s I and Lee’s L curves for one sample, one may
say, “Interesting, so what?” So we’ll compare these curves across
multiple biological conditions and see what genes and gene pairs differ
in this multi-scale manner that may not be apparent when only one
spatial scale is analyzed. Here we get which sample is from which
biological condition

``` r
data("sample_info")
sample_info
#> # A tibble: 8 × 2
#>   sample        condition   
#>   <chr>         <chr>       
#> 1 expanded1     PI14        
#> 2 expanded2     PI14        
#> 3 expanded1_pi7 PI7         
#> 4 expanded2_pi7 PI7         
#> 5 nonexpanded1  non-expanded
#> 6 nonexpanded2  non-expanded
#> 7 yoda1         Yoda1       
#> 8 yoda2         Yoda1
```

Run the whole pipeline for all the other samples; because this takes a
while to run (1 and a half hours on my laptop with 4 cores), we can
download the results from OSF to render the vignette fast.

``` r
tbs <- Piezo1TissueBoundary("all")
tx_paths <- Piezo1TxSpots("all")
for (s in sample_info$sample[-1]) {
    tb <- tbs$geometry[tbs$sample == s]
    makeAggregates(tx_paths[s],
               out_path = s, sample_id = s,
               tech = "Xenium", tissue_boundary = tb, sides = sides,
               flip_geometry = FALSE,
               square = FALSE, BPPARAM = MulticoreParam(2, progressbar = TRUE))
    runBinAnalyses(s, file.path(s, "bin_analyses"), tissue_geometry = tb,
               min_props = min_props,
               BPPARAM = MulticoreParam(2, progressbar = TRUE))
}
```

``` r
# Download from OSF
bin_analyses_dirs <- Piezo1BinAnalyses(sample = "all")
for (s in bin_analyses_dirs[-1]) untar(s)
```

### LMM analysis for Moran’s I

Here we use linear mixed models (LMM) with a spline term to model the
Moran’s I curve, with a random effect in intercept and slope so that
each biological condition has its own intercept and slope.

``` r
dirs_use <- file.path(sample_info$sample, "bin_analyses")
df_morans <- readMoranSamples(dirs_use, sample_info = sample_info)
```

Because the non-expanded samples have many small pieces and too few
cells, their Moran’s I and Lee’s L curves are rather erratic. So they’re
excluded from the following analyses.

``` r
df_morans_lmm <- df_morans |> 
    filter(condition != "non-expanded") |> 
    dplyr::rename(value = moran, group = condition, feature = gene) |> 
    runBinLMM(BPPARAM = SerialParam())
#> Warning: There were 25 warnings in `mutate()`.
#> The first warning was:
#> ℹ In argument: `pvals = bplapply(...)`.
#> Caused by warning:
#> ! Model failed to converge with 1 negative eigenvalue: -2.4e-03
#> ℹ Run `dplyr::last_dplyr_warnings()` to see the 24 remaining warnings.
```

What proportion of genes have significant random effects?

``` r
mean(df_morans_lmm$p_random_adj < 0.05)
#> [1] 0.6
```

What proportion have significant random slopes?

``` r
mean(df_morans_lmm$p_slope_adj < 0.05)
#> [1] 0.28
```

Here we plot the top 9 genes with the most significant random effects,
meaning either the intercept or the slope differ among different
conditions.

``` r
genes_use <- df_morans_lmm |> 
    filter(p_random_adj < 0.05) |> 
    slice_max(log_p_random_adj, n = 9) |> 
    pull(feature)
```

``` r
data("ditto_colors")
```

``` r
df_morans |> 
    filter(gene %in% genes_use, condition != "non-expanded") |> 
    ggplot(aes(side, moran, group = sample)) +
    geom_line(aes(color = condition)) +
    scale_color_manual(values = ditto_colors) +
    scale_x_continuous(transform = "log2") +
    facet_wrap(~ gene) +
    labs(x = sprintf("Bin size (\u03BCm)"), y = "Moran's I")
```

![](workflow_files/figure-html/unnamed-chunk-44-1.png)

Also plot genes with the most significant random slopes, i.e. slopes of
the spline terms differ among conditions

``` r
genes_use <- df_morans_lmm |> 
    filter(p_slope_adj < 0.05) |> 
    slice_max(log_p_slope_adj, n = 9) |> 
    pull(feature)
```

``` r
df_morans |> 
    filter(gene %in% genes_use, condition != "non-expanded") |> 
    ggplot(aes(side, moran, group = sample)) +
    geom_line(aes(color = condition)) +
    scale_color_manual(values = ditto_colors) +
    scale_x_continuous(transform = "log2") +
    facet_wrap(~ gene) +
    labs(x = sprintf("Bin size (\u03BCm)"), y = "Moran's I")
```

![](workflow_files/figure-html/unnamed-chunk-46-1.png)

The difference seems small for some of the genes; the LMM here does not
take into account the different variances of Moran’s I at different bin
sizes. Plot one of these genes in space in different samples

``` r
sfes <- readBinSamples(file.path(sample_info$sample, "bin_analyses"),
                       side = 32)
```

``` r
# Reduce empty space
for (i in seq_along(sfes)) {
    sfes[[i]] <- rotateMinRect(sfes[[i]])
}
```

``` r
plotSFEs(sfes, "Lgr6", ncol = 2, widths = c(1,1), show_sizes = FALSE) &
    theme(legend.position = "none")
```

![](workflow_files/figure-html/unnamed-chunk-49-1.png)

### LMM analysis for Lee’s L

Also run LMM for Lee’s L curves, to see if the multi-scale behaviors of
Lee’s L differ among conditions. With the cutoff of 0.2 here, gene pairs
that have absolute values of Lee’s L below 0.2 for all bin sizes and all
samples will be removed.

``` r
df_lees <- readLeeSamples(dirs_use, sample_info = sample_info, cutoff = 0.2)
```

It takes a few minutes to run

``` r
df_lees_lmm <- df_lees |> 
    filter(condition != "non-expanded") |> 
    dplyr::rename(value = lee, group = condition, feature = pair) |> 
    runBinLMM(BPPARAM = SerialParam())
#> Warning: There were 1064 warnings in `mutate()`.
#> The first warning was:
#> ℹ In argument: `pvals = bplapply(...)`.
#> Caused by warning:
#> ! Model failed to converge with 1 negative eigenvalue: -7.0e-02
#> ℹ Run `dplyr::last_dplyr_warnings()` to see the 1063 remaining warnings.
```

What proportion of gene pairs have significant random effects?

``` r
mean(df_lees_lmm$p_random_adj < 0.05)
#> [1] 0.6265244
```

What proportion have significant random slopes?

``` r
mean(df_lees_lmm$p_slope_adj < 0.05)
#> [1] 0.245122
```

Plot the gene pairs with the most significant random effects

``` r
pairs_use <- df_lees_lmm |> 
    slice_max(log_p_random_adj, n = 9) |> 
    pull(feature)
```

``` r
df_lees |> 
    filter(condition != "non-expanded") |> 
    dplyr::rename(group = condition) |> 
    plotLeeSelect(df_lees_lmm, pairs_use)
```

![](workflow_files/figure-html/unnamed-chunk-55-1.png)

The differences are more pronounced here than in Moran’s I. It can be
due to the differences between PI7 (Yoda1 is also at PI7) and PI14, but
sometimes it seems that Yoda1 treatment makes the gene co-expression at
PI7 more like that at PI14.

Also plot gene pairs with the most significant random slopes

``` r
pairs_use <- df_lees_lmm |> 
    slice_max(log_p_slope_adj, n = 9) |> 
    pull(feature)
```

``` r
df_lees |> 
    filter(condition != "non-expanded") |> 
    dplyr::rename(group = condition) |> 
    plotLeeSelect(df_lees_lmm, pairs_use)
```

![](workflow_files/figure-html/unnamed-chunk-57-1.png)

Plot one of the gene pairs in space

``` r
sfes <- readBinSamples(file.path(sample_info$sample, "bin_analyses"),
                       side = 128)
```

``` r
# get rid of that one annoying bin from a smaller piece of tissue
sfes[[6]] <- findDebrisCells(sfes[[6]], distance_cutoff = 300)
sfes[[6]] <- sfes[[6]][,!sfes[[6]]$is_debris]
```

``` r
# Reduce empty space
for (i in seq_along(sfes)) {
    sfes[[i]] <- rotateMinRect(sfes[[i]])
}
```

``` r
plotSFEsBiscale(sfes, "Ldhb", "Pfkm", df_lees, ncol = 2, widths = c(1,1), 
                show_sizes = FALSE, side_use = 64, legend_index = 6)
```

![](workflow_files/figure-html/unnamed-chunk-61-1.png)

### Cell Chat

Here we find which gene pairs with significant Lee’s L random effects or
random slopes that are known interactions in CellChat’s database,
without actually using the CellChat package (which is not on CRAN or
Bioconductor).

``` r
df_lee_pathways <- getCellChatInfo(df_lees_lmm, genes = rownames(sfes[[1]]),
                                   species = "mouse")
df_lee_pathways
#> # A tibble: 9 × 39
#>   pair            p_slope p_random p_main p_slope_adj  p_random_adj p_main_adj
#>   <chr>             <dbl>    <dbl>  <dbl>       <dbl>         <dbl>      <dbl>
#> 1 Itgb1_Tnc    0.00000210 3.87e-10 0.0465   0.0000888 0.00000000509     0.0963
#> 2 Itgb1_Tnc    0.00000210 3.87e-10 0.0465   0.0000888 0.00000000509     0.0963
#> 3 Cdh1_Itgb1   0.00409    8.45e- 3 0.0637   0.0233    0.0156            0.113 
#> 4 Col4a1_Itgb1 0.350      5.17e- 5 0.0345   0.588     0.000167          0.0842
#> 5 Col4a1_Itgb1 0.350      5.17e- 5 0.0345   0.588     0.000167          0.0842
#> 6 Col4a1_Itgb1 0.350      5.17e- 5 0.0345   0.588     0.000167          0.0842
#> 7 Col4a1_Itgb1 0.350      5.17e- 5 0.0345   0.588     0.000167          0.0842
#> 8 Col4a1_Itgb1 0.350      5.17e- 5 0.0345   0.588     0.000167          0.0842
#> 9 Col4a1_Itgb1 0.350      5.17e- 5 0.0345   0.588     0.000167          0.0842
#> # ℹ 32 more variables: log_p_slope_adj <dbl>, log_p_random_adj <dbl>,
#> #   log_p_main_adj <dbl>, interaction <int>, interaction_name <chr>,
#> #   pathway_name <chr>, ligand <chr>, receptor <chr>, agonist <chr>,
#> #   antagonist <chr>, co_A_receptor <chr>, co_I_receptor <chr>,
#> #   annotation <chr>, interaction_name_2 <chr>, evidence <chr>,
#> #   is_neurotransmitter <chr>, ligand.symbol <chr>, ligand.family <chr>,
#> #   ligand.location <chr>, ligand.keyword <chr>, ligand.secreted_type <chr>, …
```

Plot these gene pairs that are in CellChat’s database

``` r
pairs_use <- df_lee_pathways |> 
    pull(pair) |> unique()
```

``` r
df_lees |> 
    filter(condition != "non-expanded") |> 
    dplyr::rename(group = condition) |> 
    plotLeeSelect(df_lees_lmm, pairs_use)
```

![](workflow_files/figure-html/unnamed-chunk-64-1.png)

Plot one of them in space

``` r
plotSFEsBiscale(sfes, "Itgb1", "Tnc", df_lees, ncol = 2, widths = c(1,1), 
                show_sizes = FALSE, side_use = 64, legend_index = 6)
```

![](workflow_files/figure-html/unnamed-chunk-65-1.png)

``` r
# Clean up
unlink(sample_info$sample, recursive = TRUE)
```

## Session info

``` r
sessionInfo()
#> R Under development (unstable) (2026-04-26 r89963)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.4 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8          LC_NUMERIC=C             
#>  [3] LC_TIME=C.UTF-8           LC_COLLATE=C.UTF-8       
#>  [5] LC_MONETARY=C.UTF-8       LC_MESSAGES=C.UTF-8      
#>  [7] LC_PAPER=C.UTF-8          LC_NAME=C.UTF-8          
#>  [9] LC_ADDRESS=C.UTF-8        LC_TELEPHONE=C.UTF-8     
#> [11] LC_MEASUREMENT=C.UTF-8    LC_IDENTIFICATION=C.UTF-8
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats4    stats     graphics  grDevices utils     datasets  methods  
#> [8] base     
#> 
#> other attached packages:
#>  [1] R.utils_2.13.0                  R.oo_1.27.1                    
#>  [3] R.methodsS3_1.8.2               scater_1.39.4                  
#>  [5] scuttle_1.21.6                  SingleCellExperiment_1.33.2    
#>  [7] SummarizedExperiment_1.41.1     Biobase_2.71.0                 
#>  [9] GenomicRanges_1.63.2            Seqinfo_1.1.0                  
#> [11] IRanges_2.45.0                  S4Vectors_0.49.3               
#> [13] BiocGenerics_0.57.1             generics_0.1.4                 
#> [15] MatrixGenerics_1.23.0           matrixStats_1.5.0              
#> [17] tidyr_1.3.2                     purrr_1.2.2                    
#> [19] stringr_1.6.0                   ggplot2_4.0.3                  
#> [21] readr_2.2.0                     dplyr_1.2.1                    
#> [23] alabaster.sfe_1.3.0             alabaster.base_1.11.4          
#> [25] BiocParallel_1.45.0             sf_1.1-0                       
#> [27] Voyager_1.13.1                  SpatialFeatureExperiment_1.13.2
#> [29] Wayfarer_0.99.0                 BiocStyle_2.39.0               
#> 
#> loaded via a namespace (and not attached):
#>   [1] fs_2.1.0                  spatialreg_1.4-3         
#>   [3] bitops_1.0-9              EBImage_4.53.0           
#>   [5] httr_1.4.8                RColorBrewer_1.1-3       
#>   [7] numDeriv_2016.8-1.1       tools_4.7.0              
#>   [9] backports_1.5.1           utf8_1.2.6               
#>  [11] R6_2.6.1                  HDF5Array_1.39.1         
#>  [13] rhdf5filters_1.23.3       withr_3.0.2              
#>  [15] sp_2.2-1                  gridExtra_2.3            
#>  [17] cli_3.6.6                 textshaping_1.0.5        
#>  [19] RBioFormats_1.11.0        isoband_0.3.0            
#>  [21] sandwich_3.1-1            labeling_0.4.3           
#>  [23] marginaleffects_0.32.0    alabaster.se_1.11.0      
#>  [25] sass_0.4.10               mvtnorm_1.3-7            
#>  [27] arrow_23.0.1.2            S7_0.2.2                 
#>  [29] proxy_0.4-29              pkgdown_2.2.0            
#>  [31] systemfonts_1.3.2         scico_1.5.0              
#>  [33] limma_3.67.3              RSQLite_2.4.6            
#>  [35] httpcode_0.3.0            vroom_1.7.1              
#>  [37] spdep_1.4-2               Matrix_1.7-5             
#>  [39] ggbeeswarm_0.7.3          abind_1.4-8              
#>  [41] terra_1.9-11              lifecycle_1.0.5          
#>  [43] multcomp_1.4-30           yaml_2.3.12              
#>  [45] edgeR_4.9.9               rhdf5_2.55.16            
#>  [47] SparseArray_1.11.13       BiocFileCache_3.1.0      
#>  [49] grid_4.7.0                blob_1.3.0               
#>  [51] dqrng_0.4.1               crayon_1.5.3             
#>  [53] alabaster.spatial_1.11.1  lattice_0.22-9           
#>  [55] beachmat_2.27.5           magick_2.9.1             
#>  [57] zeallot_0.2.0             pillar_1.11.1            
#>  [59] knitr_1.51                rjson_0.2.23             
#>  [61] osfr_0.2.9                boot_1.3-32              
#>  [63] sfarrow_0.4.1             codetools_0.2-20         
#>  [65] wk_0.9.5                  glue_1.8.1               
#>  [67] data.table_1.18.2.1       memuse_4.2-3             
#>  [69] urltools_1.7.3.1          vctrs_0.7.3              
#>  [71] png_0.1-9                 Rdpack_2.6.6             
#>  [73] gtable_0.3.6              assertthat_0.2.1         
#>  [75] cachem_1.1.0              xfun_0.57                
#>  [77] rbibutils_2.4.1           S4Arrays_1.11.1          
#>  [79] DropletUtils_1.31.1       coda_0.19-4.1            
#>  [81] reformulas_0.4.4          survival_3.8-6           
#>  [83] sfheaders_0.4.5           rJava_1.0-18             
#>  [85] units_1.0-1               statmod_1.5.1            
#>  [87] bluster_1.21.1            TH.data_1.1-5            
#>  [89] nlme_3.1-169              bit64_4.8.0              
#>  [91] alabaster.ranges_1.11.0   filelock_1.0.3           
#>  [93] bslib_0.10.0              irlba_2.3.7              
#>  [95] vipor_0.4.7               KernSmooth_2.23-26       
#>  [97] spData_2.3.4              DBI_1.3.0                
#>  [99] tidyselect_1.2.1          bit_4.6.0                
#> [101] compiler_4.7.0            curl_7.1.0               
#> [103] httr2_1.2.2               BiocNeighbors_2.5.4      
#> [105] h5mread_1.3.3             xml2_1.5.2               
#> [107] desc_1.4.3                DelayedArray_0.37.1      
#> [109] triebeard_0.4.1           bookdown_0.46            
#> [111] scales_1.4.0              classInt_0.4-11          
#> [113] rappdirs_0.3.4            tiff_0.1-12              
#> [115] SpatialExperiment_1.21.0  digest_0.6.39            
#> [117] fftwtools_0.9-11          minqa_1.2.8              
#> [119] alabaster.matrix_1.11.0   rmarkdown_2.31           
#> [121] XVector_0.51.0            htmltools_0.5.9          
#> [123] pkgconfig_2.0.3           jpeg_0.1-11              
#> [125] lme4_2.0-1                sparseMatrixStats_1.23.0 
#> [127] dbplyr_2.5.2              fastmap_1.2.0            
#> [129] rlang_1.2.0               htmlwidgets_1.6.4        
#> [131] DelayedMatrixStats_1.33.0 farver_2.1.2             
#> [133] jquerylib_0.1.4           zoo_1.8-15               
#> [135] biscale_1.1.0             jsonlite_2.0.0           
#> [137] BiocSingular_1.27.1       RCurl_1.98-1.18          
#> [139] magrittr_2.0.5            s2_1.1.9                 
#> [141] patchwork_1.3.2           Rhdf5lib_1.99.6          
#> [143] Rcpp_1.1.1-1.1            ggnewscale_0.5.2         
#> [145] viridis_0.6.5             stringi_1.8.7            
#> [147] alabaster.schemas_1.11.0  MASS_7.3-65              
#> [149] parallel_4.7.0            ggrepel_0.9.8            
#> [151] deldir_2.0-4              splines_4.7.0            
#> [153] hms_1.1.4                 locfit_1.5-9.12          
#> [155] igraph_2.3.0              ScaledMatrix_1.19.0      
#> [157] LearnBayes_2.15.2         crul_1.6.0               
#> [159] evaluate_1.0.5            BiocManager_1.30.27      
#> [161] tzdb_0.5.0                nloptr_2.2.1             
#> [163] alabaster.sce_1.11.0      rsvd_1.0.5               
#> [165] e1071_1.7-17              RSpectra_0.16-2          
#> [167] viridisLite_0.4.3         class_7.3-23             
#> [169] ragg_1.5.2                tibble_3.3.1             
#> [171] lmerTest_3.2-1            memoise_2.0.1            
#> [173] beeswarm_0.4.0            cluster_2.1.8.2
```
