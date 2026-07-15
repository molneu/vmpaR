# Internal utilities for reducing multiple VMPA signatures per target.
# These functions are not user-facing and should not be exported.

.vmpa_select_unique_candidates <- function(gene_set_list) {

  if (!is.list(gene_set_list)) {
    stop("`gene_set_list` must be a list.", call. = FALSE)
  }

  if (length(gene_set_list) == 0L) {
    return(gene_set_list)
  }

  gene_set_names <- names(gene_set_list)

  if (is.null(gene_set_names)) {
    stop("`gene_set_list` must be a named list.", call. = FALSE)
  }

  if (anyNA(gene_set_names) || any(gene_set_names == "")) {
    stop("`gene_set_list` must not contain missing or empty names.", call. = FALSE)
  }

  if (anyDuplicated(gene_set_names) > 0L) {
    dup_names <- unique(gene_set_names[duplicated(gene_set_names)])
    stop(
      "`gene_set_list` contains duplicated names. Examples: ",
      paste(utils::head(dup_names, 10), collapse = ", "),
      call. = FALSE
    )
  }

  signature_metadata <- attr(gene_set_list, "signature_metadata")

  if (is.null(signature_metadata) || !is.data.frame(signature_metadata)) {
    stop(
      "`signature_metadata` is required for `unique = TRUE`.",
      call. = FALSE
    )
  }

  required_cols <- c("gene_set_name", "cmap_name", "cps_conf_total")
  missing_cols <- setdiff(required_cols, colnames(signature_metadata))

  if (length(missing_cols) > 0L) {
    stop(
      "`signature_metadata` is missing required columns for `unique = TRUE`: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  if (anyNA(signature_metadata$gene_set_name) ||
      any(signature_metadata$gene_set_name == "")) {
    stop(
      "`signature_metadata$gene_set_name` contains missing or empty values.",
      call. = FALSE
    )
  }

  if (anyDuplicated(signature_metadata$gene_set_name) > 0L) {
    dup_names <- unique(signature_metadata$gene_set_name[
      duplicated(signature_metadata$gene_set_name)
    ])
    stop(
      "`signature_metadata$gene_set_name` contains duplicated values. Examples: ",
      paste(utils::head(dup_names, 10), collapse = ", "),
      call. = FALSE
    )
  }

  if (anyNA(signature_metadata$cmap_name) ||
      any(signature_metadata$cmap_name == "")) {
    stop("`cmap_name` contains missing or empty values.", call. = FALSE)
  }

  cps_conf_total <- suppressWarnings(as.numeric(signature_metadata$cps_conf_total))

  if (anyNA(cps_conf_total) ||
      any(!is.finite(cps_conf_total)) ||
      any(cps_conf_total != floor(cps_conf_total))) {
    stop(
      "`cps_conf_total` contains missing, non-finite, or non-integer values.",
      call. = FALSE
    )
  }

  signature_metadata$cps_conf_total <- as.integer(cps_conf_total)

  if (!setequal(signature_metadata$gene_set_name, gene_set_names)) {
    stop(
      "`signature_metadata$gene_set_name` must match names(gene_set_list).",
      call. = FALSE
    )
  }

  signature_metadata <- signature_metadata[
    match(gene_set_names, signature_metadata$gene_set_name),
    ,
    drop = FALSE
  ]

  validation_col <- .vmpa_find_validation_column(signature_metadata)

  if (!is.null(validation_col)) {
    signature_metadata$.validated <- .vmpa_as_validation_flag(
      signature_metadata[[validation_col]]
    )
  } else {
    signature_metadata$.validated <- rep(FALSE, nrow(signature_metadata))
  }

  split_meta <- split(
    signature_metadata,
    factor(
      signature_metadata$cmap_name,
      levels = unique(signature_metadata$cmap_name)
    ),
    drop = TRUE
  )

  selected_meta <- do.call(
    rbind,
    lapply(split_meta, function(x) {
      x <- x[order(x$gene_set_name), , drop = FALSE]
      
      # 1) Prefer validated signatures if validation metadata are available.
      if (any(x$.validated, na.rm = TRUE)) {
        x <- x[x$.validated %in% TRUE, , drop = FALSE]
      }
      
      # 2) Keep the highest-confidence candidates.
      max_conf <- max(x$cps_conf_total, na.rm = TRUE)
      
      if (is.finite(max_conf)) {
        x <- x[x$cps_conf_total == max_conf, , drop = FALSE]
      }
      
      # 3) Keep all equally prioritized candidates.
      x
    })
  )

  rownames(selected_meta) <- selected_meta$gene_set_name

  selected_gene_sets <- gene_set_list[selected_meta$gene_set_name]

  selected_meta$.validated <- NULL

  attr(selected_gene_sets, "signature_metadata") <- selected_meta
  attr(selected_gene_sets, "unique_selection") <- selected_meta

  selected_gene_sets
}

.vmpa_reduce_unique_gsva_scores <- function(score_mat,
                                               signature_metadata) {
  if (!(is.matrix(score_mat) || is.data.frame(score_mat))) {
    stop("`score_mat` must be a matrix or data.frame.", call. = FALSE)
  }

  score_mat <- as.matrix(score_mat)

  if (!is.numeric(score_mat)) {
    stop("`score_mat` must contain numeric values.", call. = FALSE)
  }

  if (nrow(score_mat) == 0L || ncol(score_mat) == 0L) {
    stop("`score_mat` must not be empty.", call. = FALSE)
  }

  if (anyNA(score_mat) || any(!is.finite(score_mat))) {
    stop("`score_mat` contains NA, NaN, or infinite values.", call. = FALSE)
  }

  if (is.null(rownames(score_mat)) ||
      anyNA(rownames(score_mat)) ||
      any(rownames(score_mat) == "")) {
    stop("`score_mat` must have non-missing row names.", call. = FALSE)
  }

  if (anyDuplicated(rownames(score_mat)) > 0L) {
    dup_names <- unique(rownames(score_mat)[duplicated(rownames(score_mat))])
    stop(
      "`score_mat` contains duplicated row names. Examples: ",
      paste(utils::head(dup_names, 10), collapse = ", "),
      call. = FALSE
    )
  }

  if (is.null(colnames(score_mat)) ||
      anyNA(colnames(score_mat)) ||
      any(colnames(score_mat) == "")) {
    stop("`score_mat` must have non-missing column names.", call. = FALSE)
  }

  if (anyDuplicated(colnames(score_mat)) > 0L) {
    dup_names <- unique(colnames(score_mat)[duplicated(colnames(score_mat))])
    stop(
      "`score_mat` contains duplicated column names. Examples: ",
      paste(utils::head(dup_names, 10), collapse = ", "),
      call. = FALSE
    )
  }

  if (is.null(signature_metadata) || !is.data.frame(signature_metadata)) {
    stop(
      "`signature_metadata` is required to reduce GSVA scores with `unique = TRUE`.",
      call. = FALSE
    )
  }

  required_cols <- c("gene_set_name", "cmap_name")
  missing_cols <- setdiff(required_cols, colnames(signature_metadata))

  if (length(missing_cols) > 0L) {
    stop(
      "`signature_metadata` is missing required columns to reduce GSVA scores: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  if (anyNA(signature_metadata$gene_set_name) ||
      any(signature_metadata$gene_set_name == "")) {
    stop(
      "`signature_metadata$gene_set_name` contains missing or empty values.",
      call. = FALSE
    )
  }

  if (anyDuplicated(signature_metadata$gene_set_name) > 0L) {
    dup_names <- unique(signature_metadata$gene_set_name[
      duplicated(signature_metadata$gene_set_name)
    ])
    stop(
      "`signature_metadata$gene_set_name` contains duplicated values. Examples: ",
      paste(utils::head(dup_names, 10), collapse = ", "),
      call. = FALSE
    )
  }

  if (anyNA(signature_metadata$cmap_name) ||
      any(signature_metadata$cmap_name == "")) {
    stop(
      "`signature_metadata$cmap_name` contains missing or empty values.",
      call. = FALSE
    )
  }

  if (!setequal(rownames(score_mat), signature_metadata$gene_set_name)) {
    stop(
      "`signature_metadata$gene_set_name` must match rownames(score_mat).",
      call. = FALSE
    )
  }

  signature_metadata <- signature_metadata[
    match(rownames(score_mat), signature_metadata$gene_set_name),
    ,
    drop = FALSE
  ]

  target_order <- unique(signature_metadata$cmap_name)

  reduced <- lapply(target_order, function(target) {
    gene_set_names <- signature_metadata$gene_set_name[
      signature_metadata$cmap_name == target
    ]

    if (length(gene_set_names) == 1L) {
      return(as.numeric(score_mat[gene_set_names, , drop = TRUE]))
    }

    colMeans(score_mat[gene_set_names, , drop = FALSE])
  })

  reduced_mat <- do.call(rbind, reduced)
  rownames(reduced_mat) <- target_order
  colnames(reduced_mat) <- colnames(score_mat)

  reduced_mat
}

.vmpa_find_validation_column <- function(signature_metadata) {
  candidate_cols <- c(
    "validated",
    "validation",
    "validation_mark",
    "is_validated",
    "validated_signature"
  )

  found_cols <- candidate_cols[candidate_cols %in% colnames(signature_metadata)]

  if (length(found_cols) == 0L) {
    return(NULL)
  }

  found_cols[[1L]]
}

.vmpa_as_validation_flag <- function(x) {
  if (is.logical(x)) {
    return(!is.na(x) & x)
  }

  if (is.numeric(x) || is.integer(x)) {
    return(!is.na(x) & x != 0)
  }

  if (is.character(x) || is.factor(x)) {
    x_clean <- tolower(trimws(as.character(x)))
    return(x_clean %in% c("true", "t", "yes", "y", "1", "validated", "valid"))
  }

  rep(FALSE, length(x))
}
