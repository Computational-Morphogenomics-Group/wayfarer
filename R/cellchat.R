#' Filter Lee's L results for gene pairs in CellChatDB
#'
#' This function takes CellChat's database without using the CellChat package
#' itself. It filter the Lee's L LMM results (see \code{\link{runBinLMM}}) to
#' only include gene pairs in the database for known interactions, and adds info
#' about the interactions from the database such as the KEGG pathway ID. Note
#' that gene pairs that are not in the database can still be interesting because
#' they tell something about cell neighborhoods.
#'
#' @param df_lees_lmm Output for Lee's L from \code{\link{runBinLMM}}
#' @param genes A vector of gene symbols for the genes to consider; for
#'   smFISH-based technology, the list of gene of interest is often much smaller
#'   than the number of genes in the database.
#' @param species Either "human" or "mouse".
#' @param bfc A \code{\link{BiocFileCache}} instance for where the CellChat
#'   database will be downloaded.
#' @return A data frame with the LMM results and additional columns from
#'   CellChat database for more info about the interaction.
#' @importFrom stringr str_split str_to_sentence str_to_upper
#' @export
getCellChatInfo <- function(df_lees_lmm, genes, species = c("human", "mouse"),
                            bfc = BiocFileCache()) {
    species <- match.arg(species)
    ccdb <- getCellChatDB(species, bfc)

    cc_int <- str_split(ccdb$interaction$interaction_name, "_")
    if (species == "mouse")
        cc_int <- lapply(cc_int, str_to_sentence)
    cc_int_unlist <- unlist(cc_int) |> unique()
    cc_int_sub <- cc_int[map_lgl(cc_int, ~ any(genes %in% .))]
    df_lee_pathways <- df_lees_lmm |>
        separate_wider_delim(feature, "_", names = paste0("gene", 1:2)) |>
        filter(gene1 %in% cc_int_unlist, gene2 %in% cc_int_unlist)
    df_lee_pathways <- df_lee_pathways |>
        mutate(interaction = map2(gene1, gene2, function(x,y)
            which(map_lgl(cc_int_sub, function(p) x %in% p & y %in% p)))) |>
        mutate(n_pathways = lengths(interaction)) |>
        filter(n_pathways > 0) |>
        select(-n_pathways) |>
        unnest(interaction)
    cc_int_sub2 <- map_chr(cc_int_sub, paste, collapse = "_")
    if (species == "mouse")
        cc_int_sub2 <- str_to_upper(cc_int_sub2)
    df_lee_pathways |>
        mutate(interaction_name = cc_int_sub2[interaction]) |>
        left_join(CellChatDB.mouse$interaction, by = "interaction_name") |>
        unite("pair", gene1, gene2) |>
        filter(p_slope_adj < 0.05 | p_random_adj < 0.05 | p_main_adj < 0.05) |>
        arrange(p_slope_adj)
}
