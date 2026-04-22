#' Use linear mixed model to compare among groups
#'
#' This function uses a linear mixed model with a spline term to model the
#' Moran's I or Lee's L curves, so each group has its own intercept and slope,
#' distinguishing between variability within and between groups. Then the entire
#' random effects term or the random slope is dropped, and a likelihood ratio
#' test is used to compared the full model and the reduced one to see if the
#' random effects or random slope is significant. If it is significant, then it
#' indicates that the curves differ among groups, though it does not indicate
#' which groups are different.
#'
#' @inheritParams splines::bs
#' @param df_res A data frame with Moran's I or Lee's L across bin sizes and
#'   samples, such as from \code{\link{readMoranSamples}} or
#'   \code{\link{readLeeSamples}}. The columns should be renamed so that the
#'   gene or gene pair column is named "feature", the Moran's I or Lee's L
#'   values column is named "value", and the group column from additional sample
#'   info is named "group".
#' @param BPPARAM A \code{bpparam} object to parallelize computation over
#'   features.
#' @return A data frame with p-values and adjusted p-values for each feature
#' @importFrom splines bs
#' @importFrom lmerTest ranova lmer
#' @importFrom dplyr group_nest across
#' @importFrom tidyr unnest pivot_wider
#' @export
runBinLMM <- function(df_res, degree = 2, BPPARAM = SerialParam()) {
    df_res2 <- df_res |>
        mutate(x = log2(side)) |>
        group_by(feature) |>
        group_nest()
    df_res2 <- df_res2 |>
        mutate(
            pvals = bplapply(data, function(d) {
                s1 <- lmer(value ~ bs(x, degree = 2) + (1+bs(x, degree = 2) | group), data = d)
                res_slope <- ranova(s1, reduce.terms = TRUE)
                res_int <- ranova(s1, reduce.terms = FALSE)
                res_main <- anova(s1)
                tibble(p_slope = res_slope$`Pr(>Chisq)`[2],
                       p_random = res_int$`Pr(>Chisq)`[2],
                       p_main = res_main$`Pr(>F)`)
            }, BPPARAM = BPPARAM))
    df_res2 <- df_res2 |>
        select(feature, pvals) |>
        unnest(pvals)
    df_res2 <- df_res2 |>
        mutate(across(starts_with("p_"), ~ p.adjust(.x, method = "BH"), .names = "{.col}_adj"),
               across(ends_with("_adj"), ~ -log10(.x), .names = "log_{.col}"))
    df_res2
}
