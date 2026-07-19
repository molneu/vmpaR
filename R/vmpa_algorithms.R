# Internal analysis-algorithm utilities for vmpaR
# These functions are not user-facing and should not be exported.

.vmpa_validate_seed <- function(seed) {
  if (is.null(seed)) {
    return(NULL)
  }

  if (!is.numeric(seed) || length(seed) != 1L || is.na(seed) ||
      !is.finite(seed) || seed != trunc(seed) ||
      seed < -.Machine$integer.max || seed > .Machine$integer.max) {
    stop("`seed` must be NULL or a single non-missing integer.", call. = FALSE)
  }

  as.integer(seed)
}

.vmpa_run_gsva <- function(expr_mat,
                              gene_set_list,
                              min_size = 1L,
                              score_scaling = c(
                                "none",
                                "sample_z",
                                "sample_pop_sd",
                                "signature_z"
                              ),
                              verbose = TRUE) {
  score_scaling <- match.arg(score_scaling)
  min_size <- .vmpa_validate_min_size(min_size, "gsva_min_size")

  if (!requireNamespace("GSVA", quietly = TRUE)) {
    stop("Package `GSVA` must be installed for `algorithm = \"gsva\"`.", call. = FALSE)
  }

  expr_mat <- .vmpa_prepare_gsva_input(expr_mat)

  if (!is.list(gene_set_list)) {
    stop("`gene_set_list` must be a list.", call. = FALSE)
  }

  if (length(gene_set_list) == 0L) {
    stop("`gene_set_list` is empty.", call. = FALSE)
  }

  if (is.null(names(gene_set_list)) || any(names(gene_set_list) == "")) {
    stop("`gene_set_list` must be a named list.", call. = FALSE)
  }

  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("`verbose` must be TRUE or FALSE.", call. = FALSE)
  }

  gsva_par <- GSVA::gsvaParam(
    exprData = expr_mat,
    geneSets = gene_set_list,
    minSize = min_size,
    maxDiff = TRUE
  )

  vmpa_result <- GSVA::gsva(gsva_par, verbose = verbose)

  vmpa_result <- as.matrix(vmpa_result)
  attr(vmpa_result, "geneSets") <- NULL

  if (score_scaling == "none") {
    .vmpa_msg(verbose, "No GSVA score scaling applied.")

  } else if (score_scaling == "sample_z") {
    .vmpa_msg(verbose, "Applying sample-wise z-score scaling to GSVA scores.")
    .vmpa_msg(verbose, "Interpretation: emphasizes relative activity patterns within each sample.")
    vmpa_result <- base::scale(vmpa_result)
    vmpa_result <- as.matrix(vmpa_result)

  } else if (score_scaling == "sample_pop_sd") {
    .vmpa_msg(verbose, "Applying sample-wise population-SD scaling to GSVA scores.")
    .vmpa_msg(verbose, "Interpretation: emphasizes relative activity patterns within each sample using population-SD scaling.")
    N <- nrow(vmpa_result)
    vmpa_result <- base::scale(
      vmpa_result,
      scale = apply(vmpa_result, 2, stats::sd) * sqrt((N - 1) / N)
    )
    vmpa_result <- as.matrix(vmpa_result)

  } else if (score_scaling == "signature_z") {
    .vmpa_msg(verbose, "Applying signature-wise z-score scaling to GSVA scores.")
    .vmpa_msg(verbose, "Interpretation: emphasizes relative differences for each signature across samples.")
    vmpa_result <- t(base::scale(t(vmpa_result)))
    vmpa_result <- as.matrix(vmpa_result)
  }

  vmpa_result
}

.vmpa_run_fgsea <- function(stats_vec,
                               gene_set_list,
                               context,
                               gs_size,
                               conf,
                               min_size = 10L,
                               max_size = 500L,
                               n_perm_simple = 5000L,
                               seed = 123L) {
  if (!requireNamespace("fgsea", quietly = TRUE)) {
    stop("Package `fgsea` must be installed for `algorithm = \"fgsea\"`.", call. = FALSE)
  }

  stats_vec <- .vmpa_prepare_fgsea_input(stats_vec)

  if (!is.list(gene_set_list)) {
    stop("`gene_set_list` must be a list.", call. = FALSE)
  }

  if (length(gene_set_list) == 0L) {
    stop("`gene_set_list` is empty.", call. = FALSE)
  }

  if (is.null(names(gene_set_list)) || any(names(gene_set_list) == "")) {
    stop("`gene_set_list` must be a named list.", call. = FALSE)
  }

  if (!is.character(context) || length(context) != 1L || is.na(context) || context == "") {
    stop("`context` must be a single non-empty character string.", call. = FALSE)
  }

  if (!is.numeric(gs_size) || length(gs_size) != 1L || is.na(gs_size) ||
      gs_size <= 0 || gs_size != as.integer(gs_size)) {
    stop("`gs_size` must be a single positive integer.", call. = FALSE)
  }
  gs_size <- as.integer(gs_size)

  if (!is.numeric(conf) || length(conf) != 1L || is.na(conf) ||
      conf != as.integer(conf) || !conf %in% c(1L, 2L, 3L)) {
    stop("`conf` must be one of: 1, 2, 3.", call. = FALSE)
  }
  conf <- as.integer(conf)

  min_size <- .vmpa_validate_min_size(min_size, "fgsea_min_size")

  if (!is.numeric(max_size) || length(max_size) != 1L || is.na(max_size) || max_size <= 0) {
    stop("`max_size` must be a single positive number.", call. = FALSE)
  }

  if (!is.numeric(n_perm_simple) || length(n_perm_simple) != 1L || is.na(n_perm_simple) || n_perm_simple <= 0) {
    stop("`n_perm_simple` must be a single positive number.", call. = FALSE)
  }

  max_size <- as.integer(max_size)
  n_perm_simple <- as.integer(n_perm_simple)
  seed <- .vmpa_validate_seed(seed)

  if (min_size > max_size) {
    stop("`min_size` must be <= `max_size`.", call. = FALSE)
  }

  run_fgsea <- function() {
    fgsea::fgseaMultilevel(
      pathways = gene_set_list,
      stats = stats_vec,
      minSize = min_size,
      maxSize = max_size,
      nPermSimple = n_perm_simple
    )
  }

  fgsea_result <- if (is.null(seed)) {
    run_fgsea()
  } else {
    withr::with_seed(seed, run_fgsea())
  }

  vmpa_result <- as.data.frame(fgsea_result, stringsAsFactors = FALSE)
  vmpa_result$context <- context
  vmpa_result$gs_size <- as.integer(gs_size)
  vmpa_result$conf <- as.integer(conf)
  vmpa_result$ref_id <- sub("_c[0-9]+$", "", vmpa_result$pathway)
  vmpa_result$conf_total <- suppressWarnings(
    as.integer(sub(".*_c([0-9]+)$", "\\1", vmpa_result$pathway))
  )

  vmpa_result <- vmpa_result[
    order(vmpa_result$padj, -abs(vmpa_result$NES)),
    ,
    drop = FALSE
  ]

  rownames(vmpa_result) <- NULL

  vmpa_result
}
