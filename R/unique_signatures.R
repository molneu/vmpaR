# Internal utilities for reducing multiple COMPASS signatures per target.
# These functions are not user-facing and should not be exported.

.compass_select_unique_candidates <- function(gene_set_list,
                                              mode = c("aggregate", "representative")) {
  mode <- match.arg(mode)
  
  if (!is.list(gene_set_list)) {
    stop("`gene_set_list` must be a list.", call. = FALSE)
  }
  
  if (length(gene_set_list) == 0L) {
    return(gene_set_list)
  }
  
  if (is.null(names(gene_set_list)) || any(names(gene_set_list) == "")) {
    stop("`gene_set_list` must be a named list.", call. = FALSE)
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
  
  if (!all(signature_metadata$gene_set_name %in% names(gene_set_list))) {
    stop(
      "Some `signature_metadata$gene_set_name` values are not present in `gene_set_list`.",
      call. = FALSE
    )
  }
  
  signature_metadata$cps_conf_total <- as.integer(signature_metadata$cps_conf_total)
  
  validation_col <- .compass_find_validation_column(signature_metadata)
  
  if (!is.null(validation_col)) {
    signature_metadata$.validated <- .compass_as_validation_flag(signature_metadata[[validation_col]])
  } else {
    signature_metadata$.validated <- rep(FALSE, nrow(signature_metadata))
  }
  
  split_meta <- split(signature_metadata, signature_metadata$cmap_name)
  
  selected_meta <- do.call(
    rbind,
    lapply(split_meta, function(x) {
      x <- x[order(x$gene_set_name), , drop = FALSE]
      
      # 1) Prefer validated signatures if validation metadata are available.
      if (any(x$.validated, na.rm = TRUE)) {
        x <- x[x$.validated %in% TRUE, , drop = FALSE]
      }
      
      # 2) Otherwise / afterward, keep the highest-confidence candidates.
      max_conf <- suppressWarnings(max(x$cps_conf_total, na.rm = TRUE))
      
      if (is.finite(max_conf)) {
        x <- x[x$cps_conf_total == max_conf, , drop = FALSE]
      }
      
      # 3) For aggregation, keep all equally prioritized candidates.
      #    For representative output, choose one deterministically.
      if (mode == "representative" && nrow(x) > 1L) {
        x <- x[order(x$gene_set_name), , drop = FALSE]
        x <- x[1L, , drop = FALSE]
      }
      
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

.compass_reduce_unique_gsva_scores <- function(score_mat,
                                               signature_metadata) {
  if (!(is.matrix(score_mat) || is.data.frame(score_mat))) {
    stop("`score_mat` must be a matrix or data.frame.", call. = FALSE)
  }
  
  score_mat <- as.matrix(score_mat)
  
  if (is.null(rownames(score_mat)) || any(rownames(score_mat) == "")) {
    stop("`score_mat` must have row names.", call. = FALSE)
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
  
  signature_metadata <- signature_metadata[
    signature_metadata$gene_set_name %in% rownames(score_mat),
    ,
    drop = FALSE
  ]
  
  if (nrow(signature_metadata) == 0L) {
    stop(
      "No signature metadata rows match the GSVA score matrix row names.",
      call. = FALSE
    )
  }
  
  target_order <- unique(signature_metadata$cmap_name)
  
  reduced <- lapply(target_order, function(target) {
    gene_set_names <- signature_metadata$gene_set_name[
      signature_metadata$cmap_name == target
    ]
    
    gene_set_names <- gene_set_names[gene_set_names %in% rownames(score_mat)]
    
    if (length(gene_set_names) == 1L) {
      return(as.numeric(score_mat[gene_set_names, , drop = TRUE]))
    }
    
    colMeans(score_mat[gene_set_names, , drop = FALSE], na.rm = TRUE)
  })
  
  reduced_mat <- do.call(rbind, reduced)
  rownames(reduced_mat) <- target_order
  colnames(reduced_mat) <- colnames(score_mat)
  
  attr(reduced_mat, "unique_selection") <- signature_metadata
  
  reduced_mat
}

.compass_find_validation_column <- function(signature_metadata) {
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

.compass_as_validation_flag <- function(x) {
  if (is.logical(x)) {
    return(!is.na(x) & x)
  }
  
  if (is.numeric(x) || is.integer(x)) {
    return(!is.na(x) & x != 0)
  }
  
  if (is.character(x)) {
    x_clean <- tolower(trimws(x))
    return(x_clean %in% c("true", "t", "yes", "y", "1", "validated", "valid"))
  }
  
  rep(FALSE, length(x))
}