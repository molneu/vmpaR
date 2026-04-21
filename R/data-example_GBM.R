#' Example glioblastoma expression matrix from GSE145128
#'
#' A small GSVA-ready example expression matrix derived from
#' \code{ExpressionSet_GSE145128.rds} in the COMPASS paper reproducibility
#' workflow.
#'
#' The dataset represents a filtered Figure-6-like relapse comparison subset:
#' samples with \code{relapse_TYPE == "eR"} were excluded, and patients
#' \code{B046} and \code{B078} were removed. The resulting example contains
#' 10 samples from 5 paired patients, with one \code{n} and one \code{cR}
#' sample per patient.
#'
#' The matrix is intended as a lightweight package example for
#' \code{compass(..., algorithm = "gsva")}.
#'
#' @format A numeric matrix with genes in rows and samples in columns.
#'
#' @details
#' Source GEO accession: \code{GSE145128}.
#'
#' The corresponding sample metadata are available in
#' \code{\link{example_gbm_sample_info}}.
#'
#' @examples
#' data(example_gbm_expr, package = "protivity")
#' dim(example_gbm_expr)
#'
#' @source
#' Derived from GEO accession \code{GSE145128}, based on the COMPASS paper
#' reproducibility workflow.
"example_gbm_expr"


#' Example sample metadata for the glioblastoma expression example
#'
#' Sample-level metadata corresponding to \code{\link{example_gbm_expr}}.
#'
#' The metadata describe the 10 filtered samples retained from GEO accession
#' \code{GSE145128} for the package example dataset. The subset contains
#' 5 paired patients, each represented by one \code{n} sample and one
#' \code{cR} sample.
#'
#' @format A data frame with 10 rows and 6 variables:
#' \describe{
#'   \item{sample_id}{Sample identifier matching the columns of
#'   \code{example_gbm_expr}.}
#'   \item{patient_ID}{Patient identifier.}
#'   \item{relapse_TYPE}{Relapse group label (\code{n} or \code{cR}).}
#'   \item{MGMT_STATUS}{MGMT status as provided in the source object.}
#'   \item{geo_accession}{GEO accession/sample accession metadata from the
#'   source object.}
#'   \item{sample_type}{Sample type field derived from
#'   \code{sample type:ch1} in the source phenotype data.}
#' }
#'
#' @details
#' Source GEO accession: \code{GSE145128}.
#'
#' The matching expression matrix is available in
#' \code{\link{example_gbm_expr}}.
#'
#' @examples
#' data(example_gbm_sample_info, package = "protivity")
#' example_gbm_sample_info
#'
#' @source
#' Derived from GEO accession \code{GSE145128}, based on the COMPASS paper
#' reproducibility workflow.
"example_gbm_sample_info"