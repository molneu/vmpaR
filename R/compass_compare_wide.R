#' Pivot COMPASS comparison results to a wide signature table
#'
#' `compass_compare_wide()` converts signature-level output from
#' `compass_compare()` into a wide table so that multiple comparisons can be read
#' side by side for the same signature.
#'
#' Typical use case:
#' - compare `Temsi_50 vs DMSO` and `Temsi_500 vs DMSO` side by side
#' - inspect whether the same signature shifts similarly or differently
#' - create a practical work table for downstream filtering/export
#'
#' @param compare_obj Either:
#'   - the full list returned by `compass_compare()`, or
#'   - a signature-level data.frame like `cmp$results_signature`
#' @param comparison_ids Optional character vector. If provided, restrict to
#'   selected comparisons in this order.
#' @param id_cols Character vector defining the identifier columns that should
#'   remain on the left side of the wide table. Default:
#'   `c("pathway", "protein", "conf_total")`.
#' @param value_cols Character vector of columns to pivot wide. Default:
#'   `c("delta")`.
#'   Useful examples:
#'   - `c("delta")`
#'   - `c("delta", "direction")`
#'   - `c("delta", "mean_case", "mean_control")`
#' @param comparison_col Character scalar. Column holding comparison labels.
#'   Default: `"comparison_id"`.
#' @param add_delta_summary Logical. If `TRUE` and `"delta"` is among
#'   `value_cols`, add summary columns such as `max_abs_delta`,
#'   `mean_abs_delta`, `n_delta_pos`, `n_delta_neg`, and `delta_pattern`.
#'   Default: `TRUE`.
#' @param sort_by Character scalar. One of:
#'   - `"max_abs_delta"`
#'   - `"mean_abs_delta"`
#'   - `"none"`
#'   Default: `"max_abs_delta"`.
#' @param decreasing Logical. Sort descending if `TRUE`. Default: `TRUE`.
#'
#' @return A wide data.frame.
#'
#' @export
compass_compare_wide <- function(compare_obj,
                                 comparison_ids = NULL,
                                 id_cols = c("pathway", "protein", "conf_total"),
                                 value_cols = c("delta"),
                                 comparison_col = "comparison_id",
                                 add_delta_summary = TRUE,
                                 sort_by = c("max_abs_delta", "mean_abs_delta", "none"),
                                 decreasing = TRUE) {
  sort_by <- match.arg(sort_by)
  
  df <- .compass_compare_wide_extract_signature_df(compare_obj)
  
  if (!is.character(id_cols) || length(id_cols) == 0L) {
    stop("`id_cols` must be a non-empty character vector.", call. = FALSE)
  }
  
  if (!is.character(value_cols) || length(value_cols) == 0L) {
    stop("`value_cols` must be a non-empty character vector.", call. = FALSE)
  }
  
  if (!is.character(comparison_col) || length(comparison_col) != 1L || is.na(comparison_col)) {
    stop("`comparison_col` must be a single character string.", call. = FALSE)
  }
  
  if (!comparison_col %in% colnames(df)) {
    stop("`comparison_col` not found in input data.", call. = FALSE)
  }
  
  missing_id <- setdiff(id_cols, colnames(df))
  if (length(missing_id) > 0L) {
    stop(
      "The following `id_cols` are missing from input data: ",
      paste(missing_id, collapse = ", "),
      call. = FALSE
    )
  }
  
  missing_val <- setdiff(value_cols, colnames(df))
  if (length(missing_val) > 0L) {
    stop(
      "The following `value_cols` are missing from input data: ",
      paste(missing_val, collapse = ", "),
      call. = FALSE
    )
  }
  
  if (!is.logical(add_delta_summary) || length(add_delta_summary) != 1L || is.na(add_delta_summary)) {
    stop("`add_delta_summary` must be TRUE or FALSE.", call. = FALSE)
  }
  
  if (!is.logical(decreasing) || length(decreasing) != 1L || is.na(decreasing)) {
    stop("`decreasing` must be TRUE or FALSE.", call. = FALSE)
  }
  
  df[[comparison_col]] <- as.character(df[[comparison_col]])
  
  if (!is.null(comparison_ids)) {
    if (!is.character(comparison_ids)) {
      stop("`comparison_ids` must be NULL or a character vector.", call. = FALSE)
    }
    
    bad_ids <- setdiff(comparison_ids, unique(df[[comparison_col]]))
    if (length(bad_ids) > 0L) {
      stop(
        "The following `comparison_ids` are not present in the input data: ",
        paste(bad_ids, collapse = ", "),
        call. = FALSE
      )
    }
    
    df <- df[df[[comparison_col]] %in% comparison_ids, , drop = FALSE]
    df[[comparison_col]] <- factor(df[[comparison_col]], levels = comparison_ids)
  } else {
    comparison_ids <- unique(df[[comparison_col]])
    df[[comparison_col]] <- factor(df[[comparison_col]], levels = comparison_ids)
  }
  
  .compass_compare_wide_check_duplicates(
    df = df,
    id_cols = id_cols,
    comparison_col = comparison_col
  )
  
  wide_list <- lapply(value_cols, function(val_col) {
    .compass_compare_wide_single_value(
      df = df,
      id_cols = id_cols,
      comparison_col = comparison_col,
      value_col = val_col,
      comparison_ids = comparison_ids
    )
  })
  
  out <- wide_list[[1]]
  if (length(wide_list) > 1L) {
    for (i in 2:length(wide_list)) {
      out <- merge(out, wide_list[[i]], by = id_cols, all = TRUE, sort = FALSE)
    }
  }
  
  if (isTRUE(add_delta_summary) && "delta" %in% value_cols) {
    out <- .compass_compare_add_delta_summary(
      wide_df = out,
      comparison_ids = comparison_ids
    )
  }
  
  out <- .compass_compare_wide_reorder_columns(
    wide_df = out,
    id_cols = id_cols,
    comparison_ids = comparison_ids,
    value_cols = value_cols,
    add_delta_summary = add_delta_summary && "delta" %in% value_cols
  )
  
  if (sort_by != "none" && sort_by %in% colnames(out)) {
    ord <- order(out[[sort_by]], decreasing = decreasing, na.last = TRUE)
    out <- out[ord, , drop = FALSE]
  }
  
  rownames(out) <- NULL
  out
}

# Internal helpers -----------------------------------------------------------

.compass_compare_wide_extract_signature_df <- function(compare_obj) {
  if (is.list(compare_obj) && "results_signature" %in% names(compare_obj)) {
    df <- compare_obj$results_signature
  } else if (is.data.frame(compare_obj)) {
    df <- compare_obj
  } else {
    stop(
      "`compare_obj` must be either the full output of `compass_compare()` ",
      "or a signature-level data.frame.",
      call. = FALSE
    )
  }
  
  if (!is.data.frame(df)) {
    stop("Could not extract a valid signature-level data.frame.", call. = FALSE)
  }
  
  required_cols <- c("comparison_id", "pathway", "protein", "conf_total")
  missing_cols <- setdiff(required_cols, colnames(df))
  
  if (length(missing_cols) > 0L) {
    stop(
      "Input data is missing required columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }
  
  df
}

.compass_compare_wide_check_duplicates <- function(df,
                                                   id_cols,
                                                   comparison_col) {
  key_df <- df[, c(id_cols, comparison_col), drop = FALSE]
  key_chr <- apply(key_df, 1, function(x) paste(x, collapse = "|||"))
  
  if (anyDuplicated(key_chr) > 0L) {
    dup_rows <- unique(key_chr[duplicated(key_chr)])
    stop(
      "Input contains duplicated rows for the same signature/comparison combination. ",
      "Examples: ",
      paste(utils::head(dup_rows, 10), collapse = " ; "),
      call. = FALSE
    )
  }
  
  invisible(NULL)
}

.compass_compare_wide_single_value <- function(df,
                                               id_cols,
                                               comparison_col,
                                               value_col,
                                               comparison_ids) {
  sub_df <- df[, c(id_cols, comparison_col, value_col), drop = FALSE]
  
  wide_df <- stats::reshape(
    sub_df,
    idvar = id_cols,
    timevar = comparison_col,
    direction = "wide"
  )
  
  old_names <- colnames(wide_df)
  new_names <- old_names
  
  for (cmp in comparison_ids) {
    old_nm <- paste0(value_col, ".", cmp)
    new_nm <- paste0(value_col, "__", cmp)
    new_names[old_names == old_nm] <- new_nm
  }
  
  colnames(wide_df) <- new_names
  rownames(wide_df) <- NULL
  wide_df
}

.compass_compare_add_delta_summary <- function(wide_df, comparison_ids) {
  delta_cols <- paste0("delta__", comparison_ids)
  delta_cols <- delta_cols[delta_cols %in% colnames(wide_df)]
  
  if (length(delta_cols) == 0L) {
    return(wide_df)
  }
  
  delta_mat <- as.matrix(wide_df[, delta_cols, drop = FALSE])
  mode(delta_mat) <- "numeric"
  
  safe_apply <- function(x, fun) {
    apply(x, 1, function(row) {
      row <- row[is.finite(row)]
      if (length(row) == 0L) return(NA_real_)
      fun(row)
    })
  }
  
  wide_df$max_abs_delta <- safe_apply(delta_mat, function(row) max(abs(row)))
  wide_df$mean_abs_delta <- safe_apply(delta_mat, function(row) mean(abs(row)))
  wide_df$n_delta_pos <- apply(delta_mat, 1, function(row) sum(row > 0, na.rm = TRUE))
  wide_df$n_delta_neg <- apply(delta_mat, 1, function(row) sum(row < 0, na.rm = TRUE))
  wide_df$n_delta_nonzero <- apply(delta_mat, 1, function(row) sum(row != 0, na.rm = TRUE))
  
  wide_df$delta_pattern <- apply(delta_mat, 1, function(row) {
    row <- row[is.finite(row)]
    
    if (length(row) == 0L) {
      return(NA_character_)
    }
    
    n_pos <- sum(row > 0)
    n_neg <- sum(row < 0)
    
    if (n_pos > 0L && n_neg == 0L) {
      return("all_positive")
    }
    
    if (n_neg > 0L && n_pos == 0L) {
      return("all_negative")
    }
    
    if (n_pos == 0L && n_neg == 0L) {
      return("all_zero")
    }
    
    "mixed"
  })
  
  wide_df
}

.compass_compare_wide_reorder_columns <- function(wide_df,
                                                  id_cols,
                                                  comparison_ids,
                                                  value_cols,
                                                  add_delta_summary = TRUE) {
  ordered_cols <- id_cols
  
  if (isTRUE(add_delta_summary)) {
    summary_cols <- c(
      "max_abs_delta",
      "mean_abs_delta",
      "n_delta_pos",
      "n_delta_neg",
      "n_delta_nonzero",
      "delta_pattern"
    )
    summary_cols <- summary_cols[summary_cols %in% colnames(wide_df)]
    ordered_cols <- c(ordered_cols, summary_cols)
  }
  
  value_wide_cols <- unlist(lapply(value_cols, function(val_col) {
    paste0(val_col, "__", comparison_ids)
  }), use.names = FALSE)
  value_wide_cols <- value_wide_cols[value_wide_cols %in% colnames(wide_df)]
  
  ordered_cols <- c(ordered_cols, value_wide_cols)
  
  extra_cols <- setdiff(colnames(wide_df), ordered_cols)
  cbind_cols <- c(ordered_cols, extra_cols)
  
  wide_df[, cbind_cols, drop = FALSE]
}