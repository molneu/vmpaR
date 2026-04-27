#' Run COMPASS for one context and one analysis algorithm
#'
#' Main user-facing function of the protivity package.
#'
#' `compass()` resolves the subset database for the selected cancer context,
#' builds COMPASS gene sets internally, and then applies either GSVA or FGSEA
#' to the user-provided query input.
#'
#' @param input Main query input.
#'   - For `algorithm = "gsva"`: numeric matrix or data.frame with genes in rows
#'     and samples in columns.
#'   - For `algorithm = "fgsea"`: named numeric vector of gene-level statistics
#'     or rankings.
#' @param context Character; preferred context, one of:
#'   `"glioma"`, `"melanoma"`, `"nsclc"`, `"gastric"`, `"ovarian"`,
#'   `"crc"`, `"breast"`, `"prostate"`, `"pdac"`, `"headneck"`.
#' @param algorithm Character scalar. Either `"gsva"` or `"fgsea"`.
#'   Default: `"gsva"`.
#' @param gsva_score_scaling Character scalar. Post-processing applied to GSVA
#'   score matrices. One of `"none"`, `"sample_wise_zscore"`,
#'   `"sample_wise_population_sd"`, or `"signature_wise_zscore"`.
#'   Default: `"none"`. Only used when `algorithm = "gsva"`.
#' @param subset_dir Optional character scalar. Local directory containing
#'   subset files such as `glioma_subset.rds`. Checked before cache/download.
#' @param cache_dir Optional character scalar. Override for the package cache
#'   directory. If `NULL`, a package-specific user cache directory is used.
#' @param download_if_missing Logical. If `TRUE`, missing subset files are
#'   downloaded and cached automatically. Default: `TRUE`.
#' @param n Integer. Number of bottom-ranked genes to include per reference
#'   signature when building COMPASS gene sets. Default: `250L`.
#' @param min_conf Integer. Minimum `cps_conf_total` required for a
#'   reference signature to be included. Must be one of `1`, `2`, or `3`.
#'   Default: `1L`.
#' @param targets Optional character vector. If provided, only signatures with
#'   matching targets are retained. Default: `NULL`.
#' @param driver_filter Logical. If `TRUE`, only signatures with
#'   `cancer_driver_summary != "None"` are retained. This is a broad filter:
#'   any non-`"None"` driver-related annotation is kept, not only canonical
#'   drivers. Default: `FALSE`.
#' @param return_gene_sets Logical. If `TRUE`, return the internally built
#'   gene sets together with the COMPASS result. Default: `FALSE`.
#' @param verbose Logical. If `TRUE`, print progress messages for subset
#'   resolution, download steps, and optional GSVA score scaling. Default: `TRUE`.
#'
#' @return
#' If `return_gene_sets = FALSE`:
#' - for `algorithm = "gsva"`: a numeric matrix of COMPASS scores
#' - for `algorithm = "fgsea"`: a data.frame of FGSEA results
#'
#' If `return_gene_sets = TRUE`, a list with:
#' - `compass_result`
#' - `gene_sets`
#' - `context`
#' - `algorithm`
#'
#' @examples
#' if (requireNamespace("Biobase", quietly = TRUE)) {
#'   data(kebir_gb, package = "protivity")
#'
#'   df <- Biobase::exprs(kebir_gb)
#'
#'   gsva_result_example <- compass(
#'     input = df,
#'     context = "glioma",
#'     algorithm = "gsva"
#'   )
#' }
#'
#' @export
compass <- function(input,
                    context,
                    algorithm = c("gsva", "fgsea"),
                    gsva_score_scaling = c(
                      "none",
                      "sample_wise_zscore",
                      "sample_wise_population_sd",
                      "signature_wise_zscore"
                    ),
                    subset_dir = NULL,
                    cache_dir = NULL,
                    download_if_missing = TRUE,
                    n = 250L,
                    min_conf = 1L,
                    targets = NULL,
                    driver_filter = FALSE,
                    return_gene_sets = FALSE,
                    verbose = TRUE) {
  algorithm <- match.arg(algorithm)
  gsva_score_scaling <- match.arg(gsva_score_scaling)
  context <- .compass_validate_context(context)

  if (!is.logical(download_if_missing) || length(download_if_missing) != 1L || is.na(download_if_missing)) {
    stop("`download_if_missing` must be TRUE or FALSE.", call. = FALSE)
  }

  if (!is.logical(return_gene_sets) || length(return_gene_sets) != 1L || is.na(return_gene_sets)) {
    stop("`return_gene_sets` must be TRUE or FALSE.", call. = FALSE)
  }

  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("`verbose` must be TRUE or FALSE.", call. = FALSE)
  }

  if (!is.numeric(n) || length(n) != 1L || is.na(n) || n <= 0 || n != as.integer(n)) {
    stop("`n` must be a single positive integer.", call. = FALSE)
  }
  n <- as.integer(n)

  if (!is.numeric(min_conf) || length(min_conf) != 1L || is.na(min_conf) ||
      min_conf != as.integer(min_conf) || !min_conf %in% c(1L, 2L, 3L)) {
    stop("`min_conf` must be one of: 1, 2, 3.", call. = FALSE)
  }
  min_conf <- as.integer(min_conf)

  if (!is.null(targets) && !is.character(targets)) {
    stop("`targets` must be NULL or a character vector.", call. = FALSE)
  }

  if (!is.logical(driver_filter) || length(driver_filter) != 1L || is.na(driver_filter)) {
    stop("`driver_filter` must be TRUE or FALSE.", call. = FALSE)
  }

  if (algorithm != "gsva" && gsva_score_scaling != "none") {
    stop(
      "`gsva_score_scaling` is only available when `algorithm = \"gsva\"`.",
      call. = FALSE
    )
  }

  # 1) Resolve subset file (local -> cache -> optional download)
  subset_file <- .compass_resolve_subset_file(
    context = context,
    subset_dir = subset_dir,
    cache_dir = cache_dir,
    download_if_missing = download_if_missing,
    verbose = verbose
  )

  # 2) Read subset GCT and build COMPASS gene sets internally
  gct <- .compass_read_subset_gct(subset_file)

  gene_set_list <- .compass_build_gene_sets(
    gct = gct,
    n = n,
    min_conf = min_conf,
    targets = targets,
    driver_filter = driver_filter
  )

  if (length(gene_set_list) == 0L) {
    stop("No COMPASS gene sets available after filtering.", call. = FALSE)
  }

  # 3) Run analysis with the selected algorithm
  if (algorithm == "gsva") {
    compass_result <- .compass_run_gsva(
      expr_mat = input,
      gene_set_list = gene_set_list,
      score_scaling = gsva_score_scaling,
      verbose = verbose
    )
  } else if (algorithm == "fgsea") {
    compass_result <- .compass_run_fgsea(
      stats_vec = input,
      gene_set_list = gene_set_list,
      context = context,
      gs_size = n,
      conf = min_conf,
      min_size = 10L,
      max_size = 500L,
      n_perm_simple = 5000L
    )
  }

  # 4) Return result
  if (isTRUE(return_gene_sets)) {
    return(list(
      compass_result = compass_result,
      gene_sets = gene_set_list,
      context = context,
      algorithm = algorithm
    ))
  }

  compass_result
}
