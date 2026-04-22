#' Example preranked vector for FGSEA-based COMPASS analysis
#'
#' A small FGSEA-ready ranked numeric vector derived from
#' \code{drugs_U251MG_withLINCsignatures.cleaned.rds} in the COMPASS paper
#' reproducibility workflow.
#'
#' The dataset corresponds to the U251MG drug perturbation signature
#' \code{AZD-8055} at \code{0.12 uM}. The resulting object is intended as a
#' lightweight package example for \code{compass(..., algorithm = "fgsea")}.
#'
#' Duplicated gene symbols in the selected source signature were excluded
#' before constructing the final ranked vector.
#'
#' @format A named numeric vector of ranked gene-level statistics.
#'
#' @details
#' Source context: COMPASS paper reproducibility workflow, Figure 4e.
#'
#' @examples
#' data(example_u251_azd8055_rank, package = "protivity")
#' length(example_u251_azd8055_rank)
#'
#' @source
#' Derived from \code{drugs_U251MG_withLINCsignatures.cleaned.rds}, based on
#' the COMPASS paper reproducibility workflow.
"example_u251_azd8055_rank"