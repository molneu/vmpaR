#' Kebir glioblastoma example dataset
#'
#' A small glioblastoma example dataset stored as a `Biobase::ExpressionSet`.
#' It contains gene-level expression values and sample-level metadata from
#' Kebir et al. and can be used for package examples and tutorials.
#'
#' Expression values can be accessed with `Biobase::exprs(kebir_gb)`.
#' Sample metadata can be accessed with `Biobase::pData(kebir_gb)`.
#'
#' @format A `Biobase::ExpressionSet` with gene-level expression values in
#'   `Biobase::exprs(kebir_gb)` and sample annotations in
#'   `Biobase::pData(kebir_gb)`.
#'
#' @source Kebir S, Ullrich V, Berger P, Dobersalske C, Langer S,
#' Rauschenbach L, et al. A Sequential Targeting Strategy Interrupts
#' AKT-Driven Subclone-Mediated Progression in Glioblastoma.
#' Clinical Cancer Research. 2023;29(2):488-500.
#' doi:10.1158/1078-0432.CCR-22-0611. PMID:36239995.
#' GEO accession: GSE145128.
#'
#' @examples
#' if (requireNamespace("Biobase", quietly = TRUE)) {
#'   data(kebir_gb, package = "vmpaR")
#'
#'   expr_mat <- Biobase::exprs(kebir_gb)
#'   sample_info <- Biobase::pData(kebir_gb)
#'
#'   dim(expr_mat)
#'   head(sample_info)
#' }
"kebir_gb"