# Internal analysis-mode utilities for protivity
# These functions are not user-facing and should not be exported.

.compass_run_gsva <- function(expr_mat,
                              gene_set_list,
                              scale = FALSE) {
  if (!requireNamespace("GSVA", quietly = TRUE)) {
    stop("Package `GSVA` must be installed for `mode = \"gsva\"`.", call. = FALSE)
  }
  
  if (!(is.matrix(expr_mat) || is.data.frame(expr_mat))) {
    stop(
      "`mode = \"gsva\"` requires `input` to be a matrix or data.frame ",
      "with genes in rows and samples in columns.",
      call. = FALSE
    )
  }
  
  expr_mat <- as.matrix(expr_mat)
  
  if (!is.numeric(expr_mat)) {
    stop("`input` must contain numeric expression values.", call. = FALSE)
  }
  
  if (is.null(rownames(expr_mat))) {
    stop("`input` must have gene names as rownames.", call. = FALSE)
  }
  
  if (is.null(colnames(expr_mat))) {
    stop("`input` must have sample names as colnames.", call. = FALSE)
  }
  
  if (anyNA(expr_mat) || any(!is.finite(expr_mat))) {
    stop("`input` contains NA, NaN, or infinite values.", call. = FALSE)
  }
  
  if (!is.list(gene_set_list)) {
    stop("`gene_set_list` must be a list.", call. = FALSE)
  }
  
  if (length(gene_set_list) == 0L) {
    stop("`gene_set_list` is empty.", call. = FALSE)
  }
  
  if (is.null(names(gene_set_list)) || any(names(gene_set_list) == "")) {
    stop("`gene_set_list` must be a named list.", call. = FALSE)
  }
  
  if (!is.logical(scale) || length(scale) != 1L || is.na(scale)) {
    stop("`scale` must be TRUE or FALSE.", call. = FALSE)
  }
  
  gsva_par <- GSVA::gsvaParam(
    exprData = expr_mat,
    geneSets = gene_set_list,
    maxDiff = TRUE
  )
  
  compass_result <- GSVA::gsva(gsva_par)
  
  compass_result <- as.matrix(compass_result)
  attr(compass_result, "geneSets") <- NULL
  
  if (isTRUE(scale)) {
    compass_result <- t(base::scale(t(compass_result)))
  }
  
  compass_result
}


.compass_run_fgsea <- function(stats_vec,
                               gene_set_list,
                               context,
                               gs_size,
                               conf,
                               min_size = 10L,
                               max_size = 500L,
                               n_perm_simple = 5000L) {
  if (!requireNamespace("fgsea", quietly = TRUE)) {
    stop("Package `fgsea` must be installed for `mode = \"fgsea\"`.", call. = FALSE)
  }
  
  if (!is.numeric(stats_vec)) {
    stop(
      "`mode = \"fgsea\"` requires `input` to be a named numeric vector.",
      call. = FALSE
    )
  }
  
  if (is.null(names(stats_vec))) {
    stop(
      "`mode = \"fgsea\"` requires a named numeric vector ",
      "(gene names as names).",
      call. = FALSE
    )
  }
  
  if (anyNA(names(stats_vec)) || any(names(stats_vec) == "")) {
    stop("`input` contains missing or empty gene names.", call. = FALSE)
  }
  
  if (anyDuplicated(names(stats_vec)) > 0L) {
    dup_names <- unique(names(stats_vec)[duplicated(names(stats_vec))])
    stop(
      "`input` contains duplicated gene names. Examples: ",
      paste(utils::head(dup_names, 10), collapse = ", "),
      call. = FALSE
    )
  }
  
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
  
  if (!is.numeric(min_size) || length(min_size) != 1L || is.na(min_size) || min_size <= 0) {
    stop("`min_size` must be a single positive number.", call. = FALSE)
  }
  
  if (!is.numeric(max_size) || length(max_size) != 1L || is.na(max_size) || max_size <= 0) {
    stop("`max_size` must be a single positive number.", call. = FALSE)
  }
  
  if (!is.numeric(n_perm_simple) || length(n_perm_simple) != 1L || is.na(n_perm_simple) || n_perm_simple <= 0) {
    stop("`n_perm_simple` must be a single positive number.", call. = FALSE)
  }
  
  min_size <- as.integer(min_size)
  max_size <- as.integer(max_size)
  n_perm_simple <- as.integer(n_perm_simple)
  
  if (min_size > max_size) {
    stop("`min_size` must be <= `max_size`.", call. = FALSE)
  }
  
  stats_vec <- stats_vec[is.finite(stats_vec) & !is.na(stats_vec)]
  
  if (length(stats_vec) == 0L) {
    stop("`input` contains no finite numeric values after filtering.", call. = FALSE)
  }
  
  stats_vec <- sort(stats_vec, decreasing = TRUE)
  
  fgsea_result <- fgsea::fgseaMultilevel(
    pathways = gene_set_list,
    stats = stats_vec,
    minSize = min_size,
    maxSize = max_size,
    nPermSimple = n_perm_simple
  )
  
  compass_result <- as.data.frame(fgsea_result, stringsAsFactors = FALSE)
  compass_result$context <- context
  compass_result$gs_size <- as.integer(gs_size)
  compass_result$conf <- as.integer(conf)
  compass_result$ref_id <- sub("_c[0-9]+$", "", compass_result$pathway)
  compass_result$conf_total <- suppressWarnings(
    as.integer(sub(".*_c([0-9]+)$", "\\1", compass_result$pathway))
  )
  
  compass_result <- compass_result[
    order(compass_result$padj, -abs(compass_result$NES)),
    ,
    drop = FALSE
  ]
  
  rownames(compass_result) <- NULL
  
  compass_result
}