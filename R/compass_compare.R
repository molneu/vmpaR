#' Compare COMPASS score patterns across experimental groups
#'
#' `compass_compare()` performs explicit group-wise comparisons on COMPASS score
#' matrices, typically produced by `compass(..., algorithm = "gsva")`.
#'
#' The main output is a transparent signature-level result table. Optional
#' protein-level aggregation can be added explicitly, but is never silently used
#' as the main result.
#'
#' Comparison specification can be done in two ways:
#'
#' 1. **Directly via `contrast_table`** (recommended for complex designs)
#' 2. **Convenience modes**:
#'    - `"pairwise"`: one case group versus one control group
#'    - `"vs_control"`: one control group versus several case groups
#'
#' @param compass_res Numeric matrix-like COMPASS result object with
#'   pathways/signatures in rows and samples in columns.
#' @param sample_info Data frame describing the samples in `compass_res`.
#'   Must contain at least:
#'   - one sample identifier column
#'   - one grouping column used for comparison
#' @param sample_col Character scalar. Column in `sample_info` containing sample
#'   names that must match `colnames(compass_res)`.
#' @param group_col Character scalar. Column in `sample_info` defining the
#'   comparison groups.
#' @param contrast_table Optional data frame defining comparisons explicitly.
#'   Must contain:
#'   - `case_group`
#'   - `control_group`
#'   Optional:
#'   - `comparison_id`
#'
#'   If `contrast_table` is supplied, it takes precedence over `mode`,
#'   `contrast`, `control_group`, and `case_groups`.
#' @param mode Character scalar. One of `"pairwise"` or `"vs_control"`.
#'   Used only if `contrast_table` is `NULL`.
#' @param contrast Character vector of length 2 for `mode = "pairwise"`.
#'   Order matters:
#'   - `contrast[1]` = case
#'   - `contrast[2]` = control
#' @param control_group Character scalar for `mode = "vs_control"`.
#' @param case_groups Optional character vector of case groups for
#'   `mode = "vs_control"`. If `NULL`, all groups except `control_group` are used.
#' @param keep_conf Optional numeric/integer vector of confidence levels to keep,
#'   e.g. `c(2, 3)`. Default: `NULL` (keep all).
#' @param include_effect_size Logical. If `TRUE`, compute a simple standardized
#'   effect size based on pooled SD where possible. Default: `TRUE`.
#' @param protein_aggregation Character scalar. One of:
#'   - `"none"`: no protein aggregation
#'   - `"by_conf"`: aggregate by protein *and* confidence level
#'   - `"all_conf"`: aggregate by protein across all confidence levels
#'   Default: `"none"`.
#' @param include_summary Logical. If `TRUE`, print a compact console summary.
#' @param print_n Integer. Number of top positive/negative rows shown per
#'   comparison in the console summary. Default: 10.
#'
#' @return A list with:
#' - `compare_info`: data frame describing the realized comparisons
#' - `results_signature`: full signature-level comparison table
#' - `results_protein`: optional protein-level aggregation table or `NULL`
#'
#' @export
compass_compare <- function(compass_res,
                            sample_info,
                            sample_col = "sample",
                            group_col = "group",
                            contrast_table = NULL,
                            mode = c("pairwise", "vs_control"),
                            contrast = NULL,
                            control_group = NULL,
                            case_groups = NULL,
                            keep_conf = NULL,
                            include_effect_size = TRUE,
                            protein_aggregation = c("none", "by_conf", "all_conf"),
                            include_summary = TRUE,
                            print_n = 10L) {
  mode <- match.arg(mode)
  protein_aggregation <- match.arg(protein_aggregation)
  
  if (!(is.matrix(compass_res) || is.data.frame(compass_res))) {
    stop("`compass_res` must be a matrix or data.frame.", call. = FALSE)
  }
  
  compass_mat <- as.matrix(compass_res)
  
  if (!is.numeric(compass_mat)) {
    stop("`compass_res` must contain numeric values.", call. = FALSE)
  }
  
  if (is.null(rownames(compass_mat))) {
    stop("`compass_res` must have pathway/signature names as rownames.", call. = FALSE)
  }
  
  if (is.null(colnames(compass_mat))) {
    stop("`compass_res` must have sample names as colnames.", call. = FALSE)
  }
  
  if (!is.data.frame(sample_info)) {
    stop("`sample_info` must be a data.frame.", call. = FALSE)
  }
  
  if (!is.character(sample_col) || length(sample_col) != 1L || is.na(sample_col)) {
    stop("`sample_col` must be a single character string.", call. = FALSE)
  }
  
  if (!is.character(group_col) || length(group_col) != 1L || is.na(group_col)) {
    stop("`group_col` must be a single character string.", call. = FALSE)
  }
  
  if (!sample_col %in% colnames(sample_info)) {
    stop("`sample_col` not found in `sample_info`.", call. = FALSE)
  }
  
  if (!group_col %in% colnames(sample_info)) {
    stop("`group_col` not found in `sample_info`.", call. = FALSE)
  }
  
  if (!is.logical(include_effect_size) || length(include_effect_size) != 1L || is.na(include_effect_size)) {
    stop("`include_effect_size` must be TRUE or FALSE.", call. = FALSE)
  }
  
  if (!is.logical(include_summary) || length(include_summary) != 1L || is.na(include_summary)) {
    stop("`include_summary` must be TRUE or FALSE.", call. = FALSE)
  }
  
  if (!is.numeric(print_n) || length(print_n) != 1L || is.na(print_n) || print_n <= 0) {
    stop("`print_n` must be a single positive number.", call. = FALSE)
  }
  print_n <- as.integer(print_n)
  
  if (!is.null(keep_conf)) {
    if (!is.numeric(keep_conf)) {
      stop("`keep_conf` must be NULL or a numeric/integer vector.", call. = FALSE)
    }
    if (any(!keep_conf %in% c(1, 2, 3))) {
      stop("`keep_conf` values must be among 1, 2, 3.", call. = FALSE)
    }
    keep_conf <- as.integer(unique(keep_conf))
  }
  
  sample_info <- .compass_compare_prepare_sample_info(
    sample_info = sample_info,
    sample_col = sample_col,
    group_col = group_col,
    sample_names = colnames(compass_mat)
  )
  
  compare_info <- .compass_compare_build_compare_info(
    sample_info = sample_info,
    sample_col = sample_col,
    group_col = group_col,
    contrast_table = contrast_table,
    mode = mode,
    contrast = contrast,
    control_group = control_group,
    case_groups = case_groups
  )
  
  pathway_vec <- rownames(compass_mat)
  protein_vec <- .compass_compare_derive_protein(pathway_vec)
  conf_vec <- .compass_compare_parse_conf(pathway_vec)
  
  if (!is.null(keep_conf)) {
    keep_idx <- conf_vec %in% keep_conf
    compass_mat <- compass_mat[keep_idx, , drop = FALSE]
    pathway_vec <- pathway_vec[keep_idx]
    protein_vec <- protein_vec[keep_idx]
    conf_vec <- conf_vec[keep_idx]
    
    if (nrow(compass_mat) == 0L) {
      stop("No signatures remain after filtering by `keep_conf`.", call. = FALSE)
    }
  }
  
  signature_results <- .compass_compare_signature_level(
    compass_mat = compass_mat,
    pathway_vec = pathway_vec,
    protein_vec = protein_vec,
    conf_vec = conf_vec,
    compare_info = compare_info,
    include_effect_size = include_effect_size
  )
  
  protein_results <- NULL
  
  if (protein_aggregation != "none") {
    protein_results <- .compass_compare_protein_level(
      signature_results = signature_results,
      protein_aggregation = protein_aggregation
    )
  }
  
  if (isTRUE(include_summary)) {
    .compass_compare_print_summary(
      compare_info = compare_info,
      signature_results = signature_results,
      protein_results = protein_results,
      print_n = print_n
    )
  }
  
  list(
    compare_info = compare_info,
    results_signature = signature_results,
    results_protein = protein_results
  )
}

# Internal helpers -----------------------------------------------------------

.compass_compare_prepare_sample_info <- function(sample_info,
                                                 sample_col,
                                                 group_col,
                                                 sample_names) {
  df <- sample_info
  
  df[[sample_col]] <- as.character(df[[sample_col]])
  df[[group_col]] <- as.character(df[[group_col]])
  
  df <- df[!is.na(df[[sample_col]]) & df[[sample_col]] != "", , drop = FALSE]
  df <- df[!is.na(df[[group_col]]) & df[[group_col]] != "", , drop = FALSE]
  
  if (anyDuplicated(df[[sample_col]]) > 0L) {
    dup <- unique(df[[sample_col]][duplicated(df[[sample_col]])])
    stop(
      "`sample_info` contains duplicated sample identifiers. Examples: ",
      paste(utils::head(dup, 10), collapse = ", "),
      call. = FALSE
    )
  }
  
  missing_in_info <- setdiff(sample_names, df[[sample_col]])
  if (length(missing_in_info) > 0L) {
    stop(
      "The following samples are present in `compass_res` but missing in `sample_info`: ",
      paste(utils::head(missing_in_info, 20), collapse = ", "),
      call. = FALSE
    )
  }
  
  df <- df[match(sample_names, df[[sample_col]]), , drop = FALSE]
  rownames(df) <- NULL
  df
}

.compass_compare_build_compare_info <- function(sample_info,
                                                sample_col = "sample",
                                                group_col = "group",
                                                contrast_table = NULL,
                                                mode = c("pairwise", "vs_control"),
                                                contrast = NULL,
                                                control_group = NULL,
                                                case_groups = NULL) {
  mode <- match.arg(mode)
  groups_present <- unique(sample_info[[group_col]])
  
  if (!is.null(contrast_table)) {
    return(.compass_compare_build_compare_info_from_table(
      sample_info = sample_info,
      sample_col = sample_col,
      group_col = group_col,
      contrast_table = contrast_table
    ))
  }
  
  if (mode == "pairwise") {
    return(.compass_compare_build_compare_info_pairwise(
      sample_info = sample_info,
      sample_col = sample_col,
      group_col = group_col,
      contrast = contrast,
      groups_present = groups_present
    ))
  }
  
  .compass_compare_build_compare_info_vs_control(
    sample_info = sample_info,
    sample_col = sample_col,
    group_col = group_col,
    control_group = control_group,
    case_groups = case_groups,
    groups_present = groups_present
  )
}

.compass_compare_build_compare_info_from_table <- function(sample_info,
                                                           sample_col,
                                                           group_col,
                                                           contrast_table) {
  if (!is.data.frame(contrast_table)) {
    stop("`contrast_table` must be a data.frame.", call. = FALSE)
  }
  
  required_cols <- c("case_group", "control_group")
  missing_cols <- setdiff(required_cols, colnames(contrast_table))
  
  if (length(missing_cols) > 0L) {
    stop(
      "`contrast_table` is missing required columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }
  
  ct <- contrast_table
  ct$case_group <- as.character(ct$case_group)
  ct$control_group <- as.character(ct$control_group)
  
  if (!"comparison_id" %in% colnames(ct)) {
    ct$comparison_id <- paste0(ct$case_group, "_vs_", ct$control_group)
  } else {
    ct$comparison_id <- as.character(ct$comparison_id)
    missing_id <- is.na(ct$comparison_id) | ct$comparison_id == ""
    ct$comparison_id[missing_id] <- paste0(
      ct$case_group[missing_id], "_vs_", ct$control_group[missing_id]
    )
  }
  
  if (anyDuplicated(ct$comparison_id) > 0L) {
    dup <- unique(ct$comparison_id[duplicated(ct$comparison_id)])
    stop(
      "`contrast_table` contains duplicated `comparison_id` values. Examples: ",
      paste(utils::head(dup, 10), collapse = ", "),
      call. = FALSE
    )
  }
  
  groups_present <- unique(sample_info[[group_col]])
  
  bad_case <- setdiff(unique(ct$case_group), groups_present)
  if (length(bad_case) > 0L) {
    stop(
      "The following `case_group` values from `contrast_table` are not present in `sample_info`: ",
      paste(bad_case, collapse = ", "),
      call. = FALSE
    )
  }
  
  bad_ctrl <- setdiff(unique(ct$control_group), groups_present)
  if (length(bad_ctrl) > 0L) {
    stop(
      "The following `control_group` values from `contrast_table` are not present in `sample_info`: ",
      paste(bad_ctrl, collapse = ", "),
      call. = FALSE
    )
  }
  
  out_list <- lapply(seq_len(nrow(ct)), function(i) {
    case_group <- ct$case_group[i]
    control_group <- ct$control_group[i]
    comparison_id <- ct$comparison_id[i]
    
    case_samples <- sample_info[[sample_col]][sample_info[[group_col]] == case_group]
    control_samples <- sample_info[[sample_col]][sample_info[[group_col]] == control_group]
    
    if (length(case_samples) == 0L) {
      stop("No samples found for case group: ", case_group, call. = FALSE)
    }
    
    if (length(control_samples) == 0L) {
      stop("No samples found for control group: ", control_group, call. = FALSE)
    }
    
    df <- data.frame(
      comparison_id = comparison_id,
      case_group = case_group,
      control_group = control_group,
      n_case = length(case_samples),
      n_control = length(control_samples),
      stringsAsFactors = FALSE
    )
    
    df$case_samples <- I(list(case_samples))
    df$control_samples <- I(list(control_samples))
    df
  })
  
  out <- do.call(rbind, out_list)
  rownames(out) <- NULL
  out
}

.compass_compare_build_compare_info_pairwise <- function(sample_info,
                                                         sample_col,
                                                         group_col,
                                                         contrast,
                                                         groups_present) {
  if (is.null(contrast) || length(contrast) != 2L || !is.character(contrast)) {
    stop(
      "For `mode = \"pairwise\"`, `contrast` must be a character vector of length 2: c(case, control).",
      call. = FALSE
    )
  }
  
  case_group <- contrast[1]
  control_group <- contrast[2]
  
  if (!case_group %in% groups_present) {
    stop("Case group not found in `sample_info`.", call. = FALSE)
  }
  
  if (!control_group %in% groups_present) {
    stop("Control group not found in `sample_info`.", call. = FALSE)
  }
  
  case_samples <- sample_info[[sample_col]][sample_info[[group_col]] == case_group]
  control_samples <- sample_info[[sample_col]][sample_info[[group_col]] == control_group]
  
  if (length(case_samples) == 0L) {
    stop("No samples found for case group: ", case_group, call. = FALSE)
  }
  
  if (length(control_samples) == 0L) {
    stop("No samples found for control group: ", control_group, call. = FALSE)
  }
  
  out <- data.frame(
    comparison_id = paste0(case_group, "_vs_", control_group),
    case_group = case_group,
    control_group = control_group,
    n_case = length(case_samples),
    n_control = length(control_samples),
    stringsAsFactors = FALSE
  )
  
  out$case_samples <- I(list(case_samples))
  out$control_samples <- I(list(control_samples))
  out
}

.compass_compare_build_compare_info_vs_control <- function(sample_info,
                                                           sample_col,
                                                           group_col,
                                                           control_group,
                                                           case_groups,
                                                           groups_present) {
  if (is.null(control_group) || length(control_group) != 1L || !is.character(control_group)) {
    stop(
      "For `mode = \"vs_control\"`, `control_group` must be a single character string.",
      call. = FALSE
    )
  }
  
  if (!control_group %in% groups_present) {
    stop("`control_group` not found in `sample_info`.", call. = FALSE)
  }
  
  if (is.null(case_groups)) {
    case_groups <- setdiff(groups_present, control_group)
  } else {
    if (!is.character(case_groups)) {
      stop("`case_groups` must be NULL or a character vector.", call. = FALSE)
    }
    bad_cases <- setdiff(case_groups, groups_present)
    if (length(bad_cases) > 0L) {
      stop(
        "The following `case_groups` are not present in `sample_info`: ",
        paste(bad_cases, collapse = ", "),
        call. = FALSE
      )
    }
    case_groups <- unique(case_groups)
  }
  
  if (length(case_groups) == 0L) {
    stop("No valid case groups available for `mode = \"vs_control\"`.", call. = FALSE)
  }
  
  out_list <- lapply(case_groups, function(case_group) {
    case_samples <- sample_info[[sample_col]][sample_info[[group_col]] == case_group]
    control_samples <- sample_info[[sample_col]][sample_info[[group_col]] == control_group]
    
    if (length(case_samples) == 0L) {
      stop("No samples found for case group: ", case_group, call. = FALSE)
    }
    
    if (length(control_samples) == 0L) {
      stop("No samples found for control group: ", control_group, call. = FALSE)
    }
    
    df <- data.frame(
      comparison_id = paste0(case_group, "_vs_", control_group),
      case_group = case_group,
      control_group = control_group,
      n_case = length(case_samples),
      n_control = length(control_samples),
      stringsAsFactors = FALSE
    )
    
    df$case_samples <- I(list(case_samples))
    df$control_samples <- I(list(control_samples))
    df
  })
  
  out <- do.call(rbind, out_list)
  rownames(out) <- NULL
  out
}

.compass_compare_signature_level <- function(compass_mat,
                                             pathway_vec,
                                             protein_vec,
                                             conf_vec,
                                             compare_info,
                                             include_effect_size = TRUE) {
  result_list <- vector("list", nrow(compare_info))
  
  for (i in seq_len(nrow(compare_info))) {
    case_group <- compare_info$case_group[i]
    control_group <- compare_info$control_group[i]
    comparison_id <- compare_info$comparison_id[i]
    
    case_samples <- compare_info$case_samples[[i]]
    control_samples <- compare_info$control_samples[[i]]
    
    case_mat <- compass_mat[, case_samples, drop = FALSE]
    control_mat <- compass_mat[, control_samples, drop = FALSE]
    
    mean_case <- rowMeans(case_mat, na.rm = TRUE)
    mean_control <- rowMeans(control_mat, na.rm = TRUE)
    
    sd_case <- if (ncol(case_mat) > 1L) {
      apply(case_mat, 1, stats::sd, na.rm = TRUE)
    } else {
      rep(NA_real_, nrow(case_mat))
    }
    
    sd_control <- if (ncol(control_mat) > 1L) {
      apply(control_mat, 1, stats::sd, na.rm = TRUE)
    } else {
      rep(NA_real_, nrow(control_mat))
    }
    
    delta <- mean_case - mean_control
    
    effect_size <- if (isTRUE(include_effect_size)) {
      .compass_compare_compute_effect_size(
        mean_case = mean_case,
        mean_control = mean_control,
        sd_case = sd_case,
        sd_control = sd_control,
        n_case = ncol(case_mat),
        n_control = ncol(control_mat)
      )
    } else {
      rep(NA_real_, length(delta))
    }
    
    df <- data.frame(
      comparison_id = comparison_id,
      case_group = case_group,
      control_group = control_group,
      pathway = pathway_vec,
      protein = protein_vec,
      conf_total = conf_vec,
      n_case = ncol(case_mat),
      n_control = ncol(control_mat),
      mean_case = as.numeric(mean_case),
      mean_control = as.numeric(mean_control),
      sd_case = as.numeric(sd_case),
      sd_control = as.numeric(sd_control),
      delta = as.numeric(delta),
      effect_size = as.numeric(effect_size),
      direction = ifelse(
        delta > 0, "higher_in_case",
        ifelse(delta < 0, "lower_in_case", "no_change")
      ),
      stringsAsFactors = FALSE
    )
    
    result_list[[i]] <- df
  }
  
  out <- do.call(rbind, result_list)
  out <- out[order(out$comparison_id, -abs(out$delta), out$pathway), , drop = FALSE]
  rownames(out) <- NULL
  out
}

.compass_compare_compute_effect_size <- function(mean_case,
                                                 mean_control,
                                                 sd_case,
                                                 sd_control,
                                                 n_case,
                                                 n_control) {
  if (n_case < 2L || n_control < 2L) {
    return(rep(NA_real_, length(mean_case)))
  }
  
  pooled_var <- (((n_case - 1) * sd_case^2) + ((n_control - 1) * sd_control^2)) /
    (n_case + n_control - 2)
  
  pooled_sd <- sqrt(pooled_var)
  
  out <- (mean_case - mean_control) / pooled_sd
  out[!is.finite(out)] <- NA_real_
  out
}

.compass_compare_protein_level <- function(signature_results,
                                           protein_aggregation = c("by_conf", "all_conf")) {
  protein_aggregation <- match.arg(protein_aggregation)
  
  if (protein_aggregation == "by_conf") {
    group_keys <- paste(
      signature_results$comparison_id,
      signature_results$protein,
      signature_results$conf_total,
      sep = "___"
    )
    
    split_df <- split(signature_results, group_keys)
    
    out_list <- lapply(split_df, function(df) {
      delta_vec <- df$delta
      
      data.frame(
        comparison_id = df$comparison_id[1],
        case_group = df$case_group[1],
        control_group = df$control_group[1],
        protein = df$protein[1],
        conf_total = df$conf_total[1],
        n_signatures = nrow(df),
        mean_delta = mean(delta_vec, na.rm = TRUE),
        median_delta = stats::median(delta_vec, na.rm = TRUE),
        min_delta = min(delta_vec, na.rm = TRUE),
        max_delta = max(delta_vec, na.rm = TRUE),
        n_pos_delta = sum(delta_vec > 0, na.rm = TRUE),
        n_neg_delta = sum(delta_vec < 0, na.rm = TRUE),
        consistency_index = abs(sum(delta_vec > 0, na.rm = TRUE) - sum(delta_vec < 0, na.rm = TRUE)) / nrow(df),
        stringsAsFactors = FALSE
      )
    })
    
    out <- do.call(rbind, out_list)
    out <- out[order(out$comparison_id, -abs(out$mean_delta), out$protein), , drop = FALSE]
    rownames(out) <- NULL
    return(out)
  }
  
  group_keys <- paste(
    signature_results$comparison_id,
    signature_results$protein,
    sep = "___"
  )
  
  split_df <- split(signature_results, group_keys)
  
  out_list <- lapply(split_df, function(df) {
    delta_vec <- df$delta
    conf_counts <- table(df$conf_total)
    conf_counts_str <- paste(
      paste0("c", names(conf_counts), ":", as.integer(conf_counts)),
      collapse = " | "
    )
    
    data.frame(
      comparison_id = df$comparison_id[1],
      case_group = df$case_group[1],
      control_group = df$control_group[1],
      protein = df$protein[1],
      conf_total = NA_integer_,
      conf_distribution = conf_counts_str,
      n_signatures = nrow(df),
      mean_delta = mean(delta_vec, na.rm = TRUE),
      median_delta = stats::median(delta_vec, na.rm = TRUE),
      min_delta = min(delta_vec, na.rm = TRUE),
      max_delta = max(delta_vec, na.rm = TRUE),
      n_pos_delta = sum(delta_vec > 0, na.rm = TRUE),
      n_neg_delta = sum(delta_vec < 0, na.rm = TRUE),
      consistency_index = abs(sum(delta_vec > 0, na.rm = TRUE) - sum(delta_vec < 0, na.rm = TRUE)) / nrow(df),
      stringsAsFactors = FALSE
    )
  })
  
  out <- do.call(rbind, out_list)
  out <- out[order(out$comparison_id, -abs(out$mean_delta), out$protein), , drop = FALSE]
  rownames(out) <- NULL
  out
}

.compass_compare_parse_conf <- function(pathways) {
  suppressWarnings(as.integer(sub(".*_c([0-9]+)$", "\\1", pathways)))
}

.compass_compare_derive_protein <- function(pathways) {
  out <- ifelse(
    grepl("_XPR", pathways),
    sub("_XPR.*$", "", pathways),
    ifelse(
      grepl("_HAHN", pathways),
      sub("_HAHN.*$", "", pathways),
      sub("_c[0-9]+$", "", pathways)
    )
  )
  
  out[is.na(out)] <- ""
  out
}

.compass_compare_print_summary <- function(compare_info,
                                           signature_results,
                                           protein_results = NULL,
                                           print_n = 10L) {
  cat("\n===========================\n")
  cat("      COMPASS compare      \n")
  cat("===========================\n")
  cat("Comparisons:", nrow(compare_info), "\n")
  cat("Signature rows:", nrow(signature_results), "\n")
  if (!is.null(protein_results)) {
    cat("Protein rows:", nrow(protein_results), "\n")
  }
  cat("print_n:", print_n, "\n")
  
  for (i in seq_len(nrow(compare_info))) {
    cmp_id <- compare_info$comparison_id[i]
    case_group <- compare_info$case_group[i]
    control_group <- compare_info$control_group[i]
    
    df_sig <- signature_results[signature_results$comparison_id == cmp_id, , drop = FALSE]
    
    cat("\n----------------------------------------\n")
    cat("Comparison:", cmp_id, "\n")
    cat("Case:    ", case_group, "\n", sep = "")
    cat("Control: ", control_group, "\n", sep = "")
    
    pos_df <- df_sig[is.finite(df_sig$delta) & df_sig$delta > 0, , drop = FALSE]
    pos_df <- pos_df[order(-pos_df$delta), , drop = FALSE]
    pos_df <- utils::head(pos_df, print_n)
    
    neg_df <- df_sig[is.finite(df_sig$delta) & df_sig$delta < 0, , drop = FALSE]
    neg_df <- neg_df[order(neg_df$delta), , drop = FALSE]
    neg_df <- utils::head(neg_df, print_n)
    
    cat("\nTop higher-in-case signatures\n")
    if (nrow(pos_df) == 0L) {
      cat("none\n")
    } else {
      print(
        data.frame(
          pathway = pos_df$pathway,
          protein = pos_df$protein,
          conf = paste0("c", pos_df$conf_total),
          delta = formatC(pos_df$delta, format = "f", digits = 3),
          stringsAsFactors = FALSE
        ),
        row.names = FALSE,
        right = FALSE
      )
    }
    
    cat("\nTop lower-in-case signatures\n")
    if (nrow(neg_df) == 0L) {
      cat("none\n")
    } else {
      print(
        data.frame(
          pathway = neg_df$pathway,
          protein = neg_df$protein,
          conf = paste0("c", neg_df$conf_total),
          delta = formatC(neg_df$delta, format = "f", digits = 3),
          stringsAsFactors = FALSE
        ),
        row.names = FALSE,
        right = FALSE
      )
    }
    
    if (!is.null(protein_results)) {
      df_prot <- protein_results[protein_results$comparison_id == cmp_id, , drop = FALSE]
      df_prot <- df_prot[order(-abs(df_prot$mean_delta)), , drop = FALSE]
      df_prot <- utils::head(df_prot, print_n)
      
      cat("\nTop protein-level rows (optional aggregation)\n")
      print(df_prot, row.names = FALSE, right = FALSE)
    }
  }
  
  invisible(NULL)
}
