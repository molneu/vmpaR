#' protivity: context-matched protein activity inference
#'
#' protivity implements COMPASS, a context-matched protein activity scoring
#' workflow for transcriptomic data.
#'
#' Use [compass()] to run COMPASS for a selected cancer context with either
#' GSVA or FGSEA. For `algorithm = "gsva"`, the input is a gene-by-sample
#' expression matrix or data frame. For `algorithm = "fgsea"`, the input is
#' a named numeric vector of gene-level statistics or rankings.
#'
#' Use [compass_gsc()] to retrieve COMPASS gene-set collections for a selected
#' context without running downstream scoring.
#'
#' @section Example data:
#' The package includes the example dataset [kebir_gb], a small glioblastoma
#' `Biobase::ExpressionSet`. Expression values can be accessed with
#' `Biobase::exprs(kebir_gb)` and used as input for GSVA-based COMPASS scoring.
#'
#' @section Citation:
#' Please cite the COMPASS paper and the protivity package when using COMPASS
#' scores in scientific work. Citation information can be retrieved with
#' `citation("protivity")`.
#'
#' @keywords internal
"_PACKAGE"