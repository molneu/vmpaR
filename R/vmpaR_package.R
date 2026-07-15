#' vmpaR: context-matched protein activity inference
#'
#' vmpaR implements VMPA, a context-matched protein activity scoring
#' workflow for transcriptomic data.
#'
#' Use [vmpa()] to run VMPA for a selected cancer context with either
#' GSVA or FGSEA. For `algorithm = "gsva"`, the input is a gene-by-sample
#' expression matrix or data frame. For `algorithm = "fgsea"`, the input is
#' a named numeric vector of gene-level statistics or rankings.
#'
#' Use [vmpa_gsc()] to retrieve VMPA gene-set collections for a selected
#' context without running downstream scoring.
#'
#' @section Reference data:
#' The context-specific VMPA reference subsets are available from Figshare at
#' <https://doi.org/10.6084/m9.figshare.32060643>.
#'
#' @section Example data:
#' The package includes the example dataset [kebir_gb], a small glioblastoma
#' `Biobase::ExpressionSet`. Expression values can be accessed with
#' `Biobase::exprs(kebir_gb)` and used as input for GSVA-based VMPA scoring.
#'
#' @section Citation:
#' Please cite the VMPA paper when using vmpaR or VMPA-derived scores
#' in scientific work.
#'
#' @import cmapR
#' @keywords internal
"_PACKAGE"
