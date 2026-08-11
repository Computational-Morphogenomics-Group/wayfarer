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
#'   info is named "group". There should also be a column "weights" due to
#'   heteroscadiscity.
#' @param BPPARAM A \code{bpparam} object to parallelize computation over
#'   features.
#' @param save_models Logical, whether to save the fitted models.
#' @return A data frame with p-values, adjusted p-values, and log likelihood
#'   ratios for each feature. If \code{save_models = TRUE}, then there will also
#'   be a list column for the fitted models and a list column \code{data} for
#'   the data used to fit the models, which can be used to get the predicted
#'   curves.
#' @importFrom splines bs
#' @importFrom lmerTest ranova lmer
#' @importFrom dplyr group_nest across
#' @importFrom tidyr unnest pivot_wider
#' @export
runBinLMM <- function(df_res, degree = 2, BPPARAM = SerialParam(), save_models = FALSE) {
    df_res2 <- df_res |>
        mutate(x = log2(side)) |>
        group_by(feature) |>
        group_nest()
    df_res2 <- df_res2 |>
        mutate(
            pvals = bplapply(data, function(d) {
                s1 <- lmer(value ~ bs(x, degree = degree) + (1+bs(x, degree = degree) | group), data = d,
                           weights = weights)
                res_slope <- ranova(s1, reduce.terms = TRUE)
                res_int <- ranova(s1, reduce.terms = FALSE)
                res_main <- anova(s1)
                o <- tibble(p_slope = res_slope$`Pr(>Chisq)`[2],
                            p_random = res_int$`Pr(>Chisq)`[2],
                            p_main = res_main$`Pr(>F)`,
                            LRT_slope = res_slope$LRT[2],
                            LRT_random = res_int$LRT[2])
                if (save_models) {
                    data_pred <- d |>
                        select(side, sample, group)
                    o <- o |>
                        mutate(model = list(s1),
                               data = list(data_pred))
                }
                o
            }, BPPARAM = BPPARAM))
    df_res2 <- df_res2 |>
        select(feature, pvals) |>
        unnest(pvals)
    df_res2 <- df_res2 |>
        mutate(across(starts_with("p_"), ~ p.adjust(.x, method = "BH"), .names = "{.col}_adj"),
               across(ends_with("_adj"), ~ -log10(.x), .names = "log_{.col}"))
    df_res2
}

#' Get fitted curves per group
#'
#' After running \code{\link{runBinLMM}} with \code{save_models = TRUE}, use this function
#' to get the predicted curves for each group.
#'
#' @param df_res Data frame output from \code{\link{runBinLMM}}, with list columns "model" and "data".
#' @return A data frame with column "pred" for the predicted curves.
#' @export
getPredCurves <- function(df_res) {
    df_res |>
        select(feature, model, data) |>
        # Sometimes values with too low a weight are not used in fitting
        mutate(pred = lapply(model, predict)) |>
        mutate(re = map2_lgl(pred, data, \(p, d) length(p) != nrow(d)),
               pred = case_when(re ~ map2(model, data, \(m, d) predict(m, newdata=d |> mutate(x=log2(side)))),
                                TRUE ~ pred)) |>
        select(-re, -model) |>
        unnest(c(pred, data)) |>
        mutate(pred = unname(pred)) |>
        select(-sample) |> distinct()
}
