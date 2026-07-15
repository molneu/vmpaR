# Internal result-formatting utilities for vmpaR.
# These functions are not user-facing unless explicitly exported.

.vmpa_format_conf <- function(x) {
  x <- suppressWarnings(as.integer(x))
  out <- ifelse(is.na(x), NA_character_, paste0("c", x))
  out
}

.vmpa_collapse_conf <- function(x) {
  x <- suppressWarnings(as.integer(x))
  x <- x[!is.na(x)]

  if (length(x) == 0L) {
    return(NA_character_)
  }

  paste0("c", sort(unique(x)), collapse = ";")
}

.vmpa_get_metadata_col <- function(x, col) {
  if (col %in% colnames(x)) {
    return(x[[col]])
  }

  rep(NA, nrow(x))
}

.vmpa_new_result <- function(x,
                                context = NULL,
                                algorithm = NULL,
                                unique = NULL,
                                signature_metadata = NULL,
                                raw_result = NULL) {
  if (!is.data.frame(x)) {
    stop("`x` must be a data.frame.", call. = FALSE)
  }

  attr(x, "context") <- context
  attr(x, "algorithm") <- algorithm
  attr(x, "unique") <- unique
  attr(x, "signature_metadata") <- signature_metadata
  attr(x, "raw_result") <- raw_result

  class(x) <- c("vmpa_result", setdiff(class(x), "vmpa_result"))

  x
}

.vmpa_format_gsva_result <- function(score_mat,
                                        signature_metadata,
                                        context,
                                        algorithm,
                                        unique) {
  if (!(is.matrix(score_mat) || is.data.frame(score_mat))) {
    stop("`score_mat` must be a matrix or data.frame.", call. = FALSE)
  }

  score_mat <- as.matrix(score_mat)

  if (!is.numeric(score_mat)) {
    stop("`score_mat` must contain numeric values.", call. = FALSE)
  }

  if (is.null(rownames(score_mat)) ||
      anyNA(rownames(score_mat)) ||
      any(rownames(score_mat) == "")) {
    stop("`score_mat` must have non-missing row names.", call. = FALSE)
  }

  if (is.null(colnames(score_mat)) ||
      anyNA(colnames(score_mat)) ||
      any(colnames(score_mat) == "")) {
    stop("`score_mat` must have non-missing column names.", call. = FALSE)
  }

  score_df <- as.data.frame(
    score_mat,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (is.null(signature_metadata) || !is.data.frame(signature_metadata)) {
    meta_df <- data.frame(
      target = rownames(score_mat),
      conf = NA_character_,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

    out <- cbind(meta_df, score_df)
    rownames(out) <- NULL

    return(.vmpa_new_result(
      out,
      context = context,
      algorithm = algorithm,
      unique = unique,
      signature_metadata = signature_metadata,
      raw_result = score_mat
    ))
  }

  required_cols <- c("gene_set_name", "cmap_name", "cps_conf_total")
  missing_cols <- setdiff(required_cols, colnames(signature_metadata))

  if (length(missing_cols) > 0L) {
    stop(
      "`signature_metadata` is missing required columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  if (isTRUE(unique)) {
    meta_list <- lapply(rownames(score_mat), function(target) {
      x <- signature_metadata[
        signature_metadata$cmap_name == target,
        ,
        drop = FALSE
      ]

      if (nrow(x) == 0L) {
        stop(
          "Missing signature metadata for target `", target, "`.",
          call. = FALSE
        )
      }

      data.frame(
        target = target,
        conf = .vmpa_collapse_conf(x$cps_conf_total),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    })

    meta_df <- do.call(rbind, meta_list)

  } else {
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

    signature_id <- .vmpa_get_metadata_col(signature_metadata, "id")

    meta_df <- data.frame(
      target = signature_metadata$cmap_name,
      conf = .vmpa_format_conf(signature_metadata$cps_conf_total),
      signature = signature_id,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }

  out <- cbind(meta_df, score_df)
  rownames(out) <- NULL

  .vmpa_new_result(
    out,
    context = context,
    algorithm = algorithm,
    unique = unique,
    signature_metadata = signature_metadata,
    raw_result = score_mat
  )
}

.vmpa_format_fgsea_result <- function(vmpa_result,
                                         context,
                                         algorithm,
                                         unique) {
  if (!is.data.frame(vmpa_result)) {
    stop("`vmpa_result` must be a data.frame.", call. = FALSE)
  }

  if (nrow(vmpa_result) == 0L) {
    return(.vmpa_new_result(
      vmpa_result,
      context = context,
      algorithm = algorithm,
      unique = unique,
      raw_result = vmpa_result
    ))
  }

  target <- if ("target" %in% colnames(vmpa_result)) {
    vmpa_result$target
  } else {
    rep(NA_character_, nrow(vmpa_result))
  }

  conf_source <- if ("cps_conf_total" %in% colnames(vmpa_result)) {
    vmpa_result$cps_conf_total
  } else if ("conf_total" %in% colnames(vmpa_result)) {
    vmpa_result$conf_total
  } else {
    rep(NA_integer_, nrow(vmpa_result))
  }

  pathway <- if ("pathway" %in% colnames(vmpa_result)) {
    vmpa_result$pathway
  } else {
    rep(NA_character_, nrow(vmpa_result))
  }

  signature <- if ("signature_id" %in% colnames(vmpa_result)) {
    vmpa_result$signature_id
  } else if ("ref_id" %in% colnames(vmpa_result)) {
    vmpa_result$ref_id
  } else {
    pathway
  }

  meta_df <- data.frame(
    target = target,
    conf = .vmpa_format_conf(conf_source),
    pathway = pathway,
    signature = signature,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  stat_cols <- c(
    "NES",
    "pval",
    "padj",
    "ES",
    "size",
    "log2err"
  )

  stat_cols <- stat_cols[stat_cols %in% colnames(vmpa_result)]

  out <- cbind(
    meta_df,
    vmpa_result[, stat_cols, drop = FALSE]
  )

  rownames(out) <- NULL

  .vmpa_new_result(
    out,
    context = context,
    algorithm = algorithm,
    unique = unique,
    raw_result = vmpa_result
  )
}

#' @export
print.vmpa_result <- function(x, ...) {
  context <- attr(x, "context")
  algorithm <- attr(x, "algorithm")
  unique <- attr(x, "unique")

  cat("VMPA result\n")

  if (!is.null(context)) {
    cat("Context: ", context, "\n", sep = "")
  }

  if (!is.null(algorithm)) {
    cat("Algorithm: ", algorithm, "\n", sep = "")
  }

  if (!is.null(unique)) {
    cat("Unique: ", unique, "\n", sep = "")
  }

  cat("Rows: ", nrow(x), "\n", sep = "")
  cat("Columns: ", ncol(x), "\n\n", sep = "")

  print.data.frame(utils::head(as.data.frame(x), 10), ...)
  invisible(x)
}
