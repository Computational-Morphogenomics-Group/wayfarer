# Here I make an example dataset from transcripts from mouse skin
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE289263

# 1. See the size reduction when I write it to Parquet
# 2. Shall I rotate the samples to minimize their bbox? What about different pieces in one replicate? I won't separate them
# 3. Find convex hulls and centroids of cells based on label in tx files
# 4. Create SFE objects, for cells and bins
# 5. Upload Parquet tx files and SFE objects to OSF, and add function to download the data

library(tidyverse)
library(arrow)
library(sf)
library(sfheaders)
library(R.utils)
library(SpatialFeatureExperiment)

ds <- open_dataset("expanded2_transcripts.csv.gz", format = "csv")
write_dataset(ds, "expanded2", format = "feather") # file is not smaller

df <- read_csv("expanded1_transcripts.csv.gz")
df <- df |> filter(cell_id != "UNASSIGNED", qv > 20) |>
    dplyr::select(cell_id, x_location, y_location)
write_parquet(df, "expanded2.parquet")
# A lot smaller after getting rid of the other stuff, but csv.gz is smaller for the purpose of downloading

tx_slim <- function(fn) {
    name <- str_remove(fn, "_transcripts.csv.gz$")
    out_name <- paste0(name, ".csv.gz")
    if (file.exists(out_name)) return(out_name)
    df <- read_csv(fn)
    df <- df |> filter(cell_id != "UNASSIGNED", qv > 20) |>
        dplyr::select(cell_id, feature_name, x_location, y_location)
    out_name2 <- paste0(name, ".csv")
    write_csv(df, out_name2)
    gzip(out_name2)
}

fns <- list.files(pattern = "*_transcripts.csv.gz")
for (f in fns) tx_slim(f)
# Great, the file sizes are greatly reduced. I think they're small enough for
# a vignette that has to be run on GitHub Actions and Bioc's website.

# Next, make the cell polygons and centroids
make_cells_tx <- function(fn) {
    df <- read_csv(fn)
    name <- str_remove(fn, "\\.csv\\.gz$")
    out_name <- paste0(name, "_cells.rds")
    if (file.exists(out_name)) return(out_name)
    df <- df |>
        group_by(cell_id) |>
        group_nest() |>
        mutate(data = map(data,
                          function(df) {
                              inds <- chull(as.matrix(df))
                              df[c(inds,inds[1]),]
                          }),
               nrow = map_int(data, nrow)) |>
        filter(nrow >= 4) |>
        select(-nrow) |>
        unnest(cols = data) |>
        sf_polygon(x = "x_location", y = "y_location", polygon_id = "cell_id")

    df$centroids <- st_centroid(df$geometry)
    saveRDS(df, out_name)
}

fns2 <- str_remove(fns, "_transcripts")
for (f in fns2) make_cells_tx(f)

# Next, create the tissue outlines

ggplot(nonexpanded2_cells, aes(geometry = centroids)) +
    geom_sf(size = 0.1)

tb <- getTissueBoundaryConcave(nonexpanded2_cells$centroids, multiple_pieces = TRUE,
                               ratio = 0.03, allow_holes = FALSE)
plot(tb)
# OK, 0.03 seems like a good ratio. I'll create and write the tissue boundaries for all samples

df <- read_csv("nonexpanded1_transcripts.csv.gz")
ggplot(df, aes(x_location, y_location)) +
    geom_point(size = 0.1)

name <- str_remove(fns, "_transcripts.csv.gz")
for (n in name) {
    f <- paste0(n, "_cells.rds")
    df <- readRDS(f)
    tb <- getTissueBoundaryConcave(df$centroids, multiple_pieces = TRUE, ratio = 0.03)
    saveRDS(tb, paste0(n, "_tb.rds"))
}

# Great. Next I'll make the aggregations (done in the vignette)

setwd("vignettes")
dirs <- list.dirs(recursive = FALSE)
cmds <- paste0("tar --exclude='bin_analyses' -czvf ", dirs, ".tar.gz ", dirs)
for (cmd in cmds) system(cmd)

cmds2 <- paste0("tar -czvf ", dirs, "_bin_analyses.tar.gz ", dirs, "/bin_analyses")
for (cmd in cmds2) system(cmd)
