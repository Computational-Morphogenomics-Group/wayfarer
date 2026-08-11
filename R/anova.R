#' Perform ANOVA to find whether spatial metric is associated with scale and
#' group
#'
#' Because the spatial metrics for the same gene or gene pairs at different bin
#' sizes are not independeent, repeated measures ANOVA is used as if the
#' different bin sizes are repeated measures. The group in biological condition
#' is the between subject factor. This function will find whether the spatial
#' metric is generally associated with scale or group, or the interaction
#' between scale and group.
#'
#' False positives are often caused by increasing variance of Moran's I and
#' Lee's L with smaller number of bins for larger bins. To mitigate this type of
#' false positives, genes/gene pairs whose ANOVA model shows significantly
#' non-normal residuals can be removed.
#'
#' @importFrom afex aov_ez
#' @importFrom emmeans emmeans
#' @importFrom performance check_normality
#' @inheritParams runBinLMM
#' @export
#' @return A data frame with p-values for bin size (side), group, and
#'   interaction of the two. The model and the estimated marginal means are also
#'   saved. The column \code{p_normality} is the p-value in the normality test.
#'   If it's small, it indicates non-normal residuals. You can inspect by taking
#'   the relevant item in the \code{model} column with the function
#'   \code{\link[performance]{check_model}}.
runAnova <- function(df_res, degree = 2, BPPARAM = SerialParam()) {
    df_res |>
        group_by(feature) |>
        group_nest() |>
        mutate(res = bplapply(data, function(d) {
            s <- aov_ez("sample", "value", d, between = "group", within = "side")
            m <- emmeans(s, "group", by = "side")
            ps <- s$anova_table$`Pr(>F)`
            tibble(p_group = ps[1],
                   p_side = ps[2],
                   p_group_side = ps[3],
                   model = list(s),
                   mmeans = list(m))
        }, BPPARAM = BPPARAM)) |>
        dplyr::select(feature, res) |>
        unnest(col = res) |>
        mutate(across(starts_with("p_"), ~ p.adjust(.x, method = "BH"), .names = "{.col}_adj"),
               across(ends_with("_adj"), ~ -log10(.x), .names = "log_{.col}"),
               p_normality = vapply(model, \(m) {
                   nn <- check_normality(m)
                   attributes(nn) <- NULL
                   nn
               }, FUN.VALUE = numeric(1)),
               p_normality_adj = p.adjust(p_normality, "BH"))
}
