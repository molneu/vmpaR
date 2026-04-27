#' Kebir glioblastoma example dataset
#'
#' A small glioblastoma example dataset stored as a Biobase ExpressionSet.
#' It contains gene-level expression values and sample-level metadata from
#' Kebir et al. and can be used for package examples and tutorials.
#'
#' Expression values can be accessed with `Biobase::exprs(kebir_gb)`.
#' Sample metadata can be accessed with `Biobase::pData(kebir_gb)`.
#'
#' @format A Biobase ExpressionSet with gene-level expression values in
#'   `Biobase::exprs(kebir_gb)` and sample annotations in
#'   `Biobase::pData(kebir_gb)`.
#'
#' @source Kebir et al.; GSE145128.
#'
#' @examples
#' if (requireNamespace("Biobase", quietly = TRUE)) {
#'   data(kebir_gb)
#'
#'   expr_mat <- Biobase::exprs(kebir_gb)
#'   sample_info <- Biobase::pData(kebir_gb)
#'
#'   dim(expr_mat)
#'   head(sample_info)
#' }
"kebir_gb"