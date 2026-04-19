#' Summarize COMPASS results
#'
#' Creates a compact summary table from a COMPASS result object and optionally
#' prints a readable summary to the console.
#'
#' Supported result types:
#' - fgsea-style data.frame returned by `compass(..., mode = "fgsea")`
#' - score matrix returned by `compass(..., mode = "gsva")`
#' - future AUCell score matrices should also work if they follow the same
#'   basic matrix structure (gene sets x samples/cells)
#'
#' @param compass_res A COMPASS result object:
#'   - fgsea-like result data.frame
#'   - matrix of COMPASS activity scores
#' @param top_n Integer. Number of top entries returned per category.
#'   Default: 100.
#' @param print_n Integer. Number of entries shown per section in the console.
#'   Default: 20.
#' @param sample Optional character scalar. For matrix-like results, if provided,
#'   positive/negative activity summaries are restricted to this sample/cell.
#'   RNA–activity discordance will also be printed for this sample only, but is
#'   still computed relative to the z-scored distribution across all overlapping
#'   samples.
#' @param expr_mat Optional expression matrix (genes x samples). Required only if
#'   `include_rna_concordance = TRUE` or `include_hypothesis_shortlist = TRUE`.
#' @param index Optional data.frame with pathway/protein metadata. If `NULL`,
#'   protein names are derived heuristically from pathway names.
#' @param include_rna_concordance Logical. If `TRUE`, compute RNA–activity
#'   discordance for matrix-like results. Default: `TRUE`.
#' @param include_variability Logical. If `TRUE`, include top variable signatures
#'   across all samples. Default: `FALSE`.
#' @param include_signature_consistency Logical. If `TRUE`, add an optional
#'   protein-level warning block for convergent vs divergent signatures within
#'   the selected sample. Currently supported for matrix-like COMPASS results.
#'   Default: `FALSE`.
#' @param include_hypothesis_shortlist Logical. If `TRUE`, add an optional
#'   hypothesis-generation shortlist based on signal strength, confidence,
#'   consistency, and support. Default: `FALSE`.
#' @param shortlist_n Integer. Number of shortlist entries per quadrant in the
#'   hypothesis matrix. Default: 3.
#' @param verbose Logical. If `TRUE`, print a compact console summary.
#'   Default: `TRUE`.
#'
#' @return A data.frame with one row per summary entry.
#'
#' @export
compass_summary <- function(compass_res,
                            top_n = 100L,
                            print_n = 20L,
                            sample = NULL,
                            expr_mat = NULL,
                            index = NULL,
                            include_rna_concordance = TRUE,
                            include_variability = FALSE,
                            include_signature_consistency = FALSE,
                            include_hypothesis_shortlist = FALSE,
                            shortlist_n = 3L,
                            verbose = TRUE) {
  if (!is.numeric(top_n) || length(top_n) != 1L || is.na(top_n) || top_n <= 0) {
    stop("`top_n` must be a single positive number.", call. = FALSE)
  }
  top_n <- as.integer(top_n)
  
  if (!is.numeric(print_n) || length(print_n) != 1L || is.na(print_n) || print_n <= 0) {
    stop("`print_n` must be a single positive number.", call. = FALSE)
  }
  print_n <- as.integer(print_n)
  
  if (!is.numeric(shortlist_n) || length(shortlist_n) != 1L || is.na(shortlist_n) || shortlist_n <= 0) {
    stop("`shortlist_n` must be a single positive number.", call. = FALSE)
  }
  shortlist_n <- as.integer(shortlist_n)
  
  if (!is.null(sample) && (!is.character(sample) || length(sample) != 1L || is.na(sample))) {
    stop("`sample` must be NULL or a single character string.", call. = FALSE)
  }
  
  if (!is.logical(include_rna_concordance) || length(include_rna_concordance) != 1L || is.na(include_rna_concordance)) {
    stop("`include_rna_concordance` must be TRUE or FALSE.", call. = FALSE)
  }
  
  if (!is.logical(include_variability) || length(include_variability) != 1L || is.na(include_variability)) {
    stop("`include_variability` must be TRUE or FALSE.", call. = FALSE)
  }
  
  if (!is.logical(include_signature_consistency) || length(include_signature_consistency) != 1L || is.na(include_signature_consistency)) {
    stop("`include_signature_consistency` must be TRUE or FALSE.", call. = FALSE)
  }
  
  if (!is.logical(include_hypothesis_shortlist) || length(include_hypothesis_shortlist) != 1L || is.na(include_hypothesis_shortlist)) {
    stop("`include_hypothesis_shortlist` must be TRUE or FALSE.", call. = FALSE)
  }
  
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("`verbose` must be TRUE or FALSE.", call. = FALSE)
  }
  
  is_fgsea_like <- is.data.frame(compass_res) &&
    all(c("pathway", "NES") %in% colnames(compass_res))
  
  if (is_fgsea_like) {
    summary_df <- .compass_summary_fgsea(
      compass_res = compass_res,
      top_n = top_n
    )
  } else if (is.matrix(compass_res) || is.data.frame(compass_res)) {
    summary_df <- .compass_summary_matrix(
      compass_res = as.matrix(compass_res),
      top_n = top_n,
      sample = sample,
      expr_mat = expr_mat,
      index = index,
      include_rna_concordance = include_rna_concordance,
      include_variability = include_variability,
      include_signature_consistency = include_signature_consistency,
      include_hypothesis_shortlist = include_hypothesis_shortlist,
      shortlist_n = shortlist_n
    )
  } else {
    stop(
      "`compass_res` must be either a fgsea-like data.frame or a score matrix.",
      call. = FALSE
    )
  }
  
  if (isTRUE(verbose)) {
    overview_info <- .compass_build_overview_info(
      compass_res = compass_res,
      summary_df = summary_df,
      sample = sample
    )
    
    .print_compass_summary(
      summary_df = summary_df,
      top_n = top_n,
      print_n = print_n,
      sample = sample,
      shortlist_n = shortlist_n,
      overview_info = overview_info
    )
  }
  
  summary_df
}

# Internal helpers -----------------------------------------------------------

.compass_summary_fgsea <- function(compass_res, top_n) {
  required_cols <- c("pathway", "NES")
  missing_cols <- setdiff(required_cols, colnames(compass_res))
  
  if (length(missing_cols) > 0L) {
    stop(
      "FGSEA-like result is missing required columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }
  
  df <- compass_res
  
  if (!"padj" %in% colnames(df)) df$padj <- NA_real_
  if (!"pval" %in% colnames(df)) df$pval <- NA_real_
  if (!"conf_total" %in% colnames(df)) {
    df$conf_total <- suppressWarnings(
      as.integer(sub(".*_c([0-9]+)$", "\\1", df$pathway))
    )
  }
  
  df$protein <- .derive_protein_from_pathway(df$pathway)
  df$sample <- NA_character_
  
  pos_df <- df[df$NES > 0, , drop = FALSE]
  pos_df <- pos_df[order(-pos_df$NES, pos_df$padj), , drop = FALSE]
  pos_df <- utils::head(pos_df, top_n)
  pos_df$summary_type <- "top_positive_nes"
  pos_df$direction <- "positive"
  pos_df$rank <- seq_len(nrow(pos_df))
  pos_df$score <- pos_df$NES
  
  neg_df <- df[df$NES < 0, , drop = FALSE]
  neg_df <- neg_df[order(neg_df$NES, neg_df$padj), , drop = FALSE]
  neg_df <- utils::head(neg_df, top_n)
  neg_df$summary_type <- "top_negative_nes"
  neg_df$direction <- "negative"
  neg_df$rank <- seq_len(nrow(neg_df))
  neg_df$score <- neg_df$NES
  
  out <- rbind(
    pos_df[, c("summary_type", "direction", "rank", "sample", "pathway", "protein", "score", "conf_total", "padj", "pval"), drop = FALSE],
    neg_df[, c("summary_type", "direction", "rank", "sample", "pathway", "protein", "score", "conf_total", "padj", "pval"), drop = FALSE]
  )
  
  rownames(out) <- NULL
  out
}

.compass_summary_matrix <- function(compass_res,
                                    top_n,
                                    sample = NULL,
                                    expr_mat = NULL,
                                    index = NULL,
                                    include_rna_concordance = FALSE,
                                    include_variability = FALSE,
                                    include_signature_consistency = FALSE,
                                    include_hypothesis_shortlist = FALSE,
                                    shortlist_n = 3L) {
  if (!is.numeric(compass_res)) {
    stop("Matrix-like COMPASS results must be numeric.", call. = FALSE)
  }
  
  if (is.null(rownames(compass_res))) {
    stop("Score matrix must have pathway names as rownames.", call. = FALSE)
  }
  
  if (is.null(colnames(compass_res))) {
    stop("Score matrix must have sample/cell names as colnames.", call. = FALSE)
  }
  
  sample_names <- colnames(compass_res)
  
  if (!is.null(sample)) {
    if (!sample %in% sample_names) {
      stop("Requested `sample` not found in `compass_res` columns.", call. = FALSE)
    }
    sample_names <- sample
  }
  
  protein_map <- .resolve_protein_map(
    pathways = rownames(compass_res),
    index = index
  )
  
  conf_total <- suppressWarnings(
    as.integer(sub(".*_c([0-9]+)$", "\\1", rownames(compass_res)))
  )
  
  summary_list <- list()
  
  # Top positive / negative scores per requested sample(s)
  for (s in sample_names) {
    vals <- compass_res[, s]
    
    df <- data.frame(
      summary_type = NA_character_,
      direction = NA_character_,
      rank = NA_integer_,
      sample = s,
      pathway = rownames(compass_res),
      protein = unname(protein_map[rownames(compass_res)]),
      score = as.numeric(vals),
      conf_total = conf_total,
      padj = NA_real_,
      pval = NA_real_,
      stringsAsFactors = FALSE
    )
    
    pos_df <- df[order(-df$score), , drop = FALSE]
    pos_df <- utils::head(pos_df, top_n)
    pos_df$summary_type <- "top_positive_activity"
    pos_df$direction <- "positive"
    pos_df$rank <- seq_len(nrow(pos_df))
    
    neg_df <- df[order(df$score), , drop = FALSE]
    neg_df <- utils::head(neg_df, top_n)
    neg_df$summary_type <- "top_negative_activity"
    neg_df$direction <- "negative"
    neg_df$rank <- seq_len(nrow(neg_df))
    
    summary_list[[paste0("pos_", s)]] <- pos_df
    summary_list[[paste0("neg_", s)]] <- neg_df
  }
  
  # Optional: top variable signatures across all samples
  if (isTRUE(include_variability)) {
    row_sd <- apply(compass_res, 1, stats::sd, na.rm = TRUE)
    
    var_df <- data.frame(
      summary_type = "top_variable_activity",
      direction = "variable",
      rank = seq_len(min(top_n, length(row_sd))),
      sample = NA_character_,
      pathway = rownames(compass_res),
      protein = unname(protein_map[rownames(compass_res)]),
      score = as.numeric(row_sd),
      conf_total = conf_total,
      padj = NA_real_,
      pval = NA_real_,
      stringsAsFactors = FALSE
    )
    
    var_df <- var_df[order(-var_df$score), , drop = FALSE]
    var_df <- utils::head(var_df, top_n)
    var_df$rank <- seq_len(nrow(var_df))
    
    summary_list[["variable"]] <- var_df
  }
  
  # Optional RNA–activity discordance
  if (isTRUE(include_rna_concordance)) {
    if (is.null(expr_mat)) {
      stop(
        "`expr_mat` must be provided when `include_rna_concordance = TRUE`.",
        call. = FALSE
      )
    }
    
    discordance_df <- .compute_rna_activity_discordance(
      compass_res = compass_res,
      expr_mat = expr_mat,
      protein_map = protein_map,
      top_n = top_n,
      sample = sample
    )
    
    if (nrow(discordance_df) > 0L) {
      summary_list[["discordance"]] <- discordance_df
    }
  }
  
  # Optional signature consistency / divergence
  if (isTRUE(include_signature_consistency)) {
    if (is.null(sample)) {
      stop(
        "`sample` must be provided when `include_signature_consistency = TRUE`.",
        call. = FALSE
      )
    }
    
    consistency_df <- .compute_signature_consistency(
      compass_res = compass_res,
      protein_map = protein_map,
      top_n = top_n,
      sample = sample
    )
    
    if (nrow(consistency_df) > 0L) {
      summary_list[["consistency"]] <- consistency_df
    }
  }
  
  # Optional hypothesis shortlist
  if (isTRUE(include_hypothesis_shortlist)) {
    if (is.null(sample)) {
      stop(
        "`sample` must be provided when `include_hypothesis_shortlist = TRUE`.",
        call. = FALSE
      )
    }
    
    if (is.null(expr_mat)) {
      stop(
        "`expr_mat` must be provided when `include_hypothesis_shortlist = TRUE`.",
        call. = FALSE
      )
    }
    
    shortlist_df <- .compute_hypothesis_shortlist(
      compass_res = compass_res,
      expr_mat = expr_mat,
      protein_map = protein_map,
      sample = sample,
      shortlist_n = shortlist_n
    )
    
    if (nrow(shortlist_df) > 0L) {
      summary_list[["hypothesis_shortlist"]] <- shortlist_df
    }
  }
  
  out <- .bind_summary_frames(summary_list)
  rownames(out) <- NULL
  out
}

.compute_rna_activity_discordance_full <- function(compass_res,
                                                   expr_mat,
                                                   protein_map,
                                                   sample = NULL) {
  if (!(is.matrix(expr_mat) || is.data.frame(expr_mat))) {
    stop("`expr_mat` must be a matrix or data.frame.", call. = FALSE)
  }
  
  expr_mat <- as.matrix(expr_mat)
  
  if (!is.numeric(expr_mat)) {
    stop("`expr_mat` must contain numeric values.", call. = FALSE)
  }
  
  if (is.null(rownames(expr_mat))) {
    stop("`expr_mat` must have gene names as rownames.", call. = FALSE)
  }
  
  common_samples <- intersect(colnames(compass_res), colnames(expr_mat))
  
  if (length(common_samples) < 2L) {
    stop(
      "RNA–activity concordance requires at least 2 overlapping samples ",
      "between `compass_res` and `expr_mat`.",
      call. = FALSE
    )
  }
  
  pathways <- rownames(compass_res)
  conf_total <- suppressWarnings(
    as.integer(sub(".*_c([0-9]+)$", "\\1", pathways))
  )
  
  rows <- vector("list", length = 0L)
  
  for (i in seq_along(pathways)) {
    pw <- pathways[i]
    protein <- protein_map[[pw]]
    
    if (is.na(protein) || !protein %in% rownames(expr_mat)) {
      next
    }
    
    activity_vals <- as.numeric(compass_res[pw, common_samples])
    rna_vals <- as.numeric(expr_mat[protein, common_samples])
    
    if (all(is.na(activity_vals)) || all(is.na(rna_vals))) {
      next
    }
    
    if (stats::sd(activity_vals, na.rm = TRUE) == 0 || stats::sd(rna_vals, na.rm = TRUE) == 0) {
      next
    }
    
    activity_z <- as.numeric(scale(activity_vals))
    rna_z <- as.numeric(scale(rna_vals))
    discordance <- activity_z - rna_z
    
    rows[[length(rows) + 1L]] <- data.frame(
      summary_type = NA_character_,
      direction = ifelse(discordance >= 0, "activity_gt_rna", "rna_gt_activity"),
      rank = NA_integer_,
      sample = common_samples,
      pathway = pw,
      protein = protein,
      score = discordance,
      conf_total = conf_total[i],
      padj = NA_real_,
      pval = NA_real_,
      stringsAsFactors = FALSE
    )
  }
  
  out <- .bind_summary_frames(rows)
  
  if (!is.null(sample) && nrow(out) > 0L) {
    out <- out[out$sample == sample, , drop = FALSE]
  }
  
  rownames(out) <- NULL
  out
}

.compute_rna_activity_discordance <- function(compass_res,
                                              expr_mat,
                                              protein_map,
                                              top_n,
                                              sample = NULL) {
  discordance_df <- .compute_rna_activity_discordance_full(
    compass_res = compass_res,
    expr_mat = expr_mat,
    protein_map = protein_map,
    sample = sample
  )
  
  if (nrow(discordance_df) == 0L) {
    return(discordance_df)
  }
  
  pos_df <- discordance_df[discordance_df$score > 0, , drop = FALSE]
  pos_df <- pos_df[order(-pos_df$score), , drop = FALSE]
  pos_df <- utils::head(pos_df, top_n)
  pos_df$summary_type <- "top_activity_gt_rna"
  pos_df$rank <- seq_len(nrow(pos_df))
  
  neg_df <- discordance_df[discordance_df$score < 0, , drop = FALSE]
  neg_df <- neg_df[order(neg_df$score), , drop = FALSE]
  neg_df <- utils::head(neg_df, top_n)
  neg_df$summary_type <- "top_rna_gt_activity"
  neg_df$rank <- seq_len(nrow(neg_df))
  
  out <- .bind_summary_frames(list(pos_df, neg_df))
  rownames(out) <- NULL
  out
}

.compute_signature_consistency <- function(compass_res,
                                           protein_map,
                                           top_n,
                                           sample) {
  if (is.null(sample) || length(sample) != 1L) {
    stop("`sample` must be a single sample name.", call. = FALSE)
  }
  
  if (!sample %in% colnames(compass_res)) {
    stop("Requested `sample` not found in `compass_res` columns.", call. = FALSE)
  }
  
  df <- data.frame(
    pathway = rownames(compass_res),
    protein = unname(protein_map[rownames(compass_res)]),
    score = as.numeric(compass_res[, sample]),
    stringsAsFactors = FALSE
  )
  
  df <- df[!is.na(df$protein) & df$protein != "" & is.finite(df$score), , drop = FALSE]
  protein_split <- split(df, df$protein)
  
  convergent_rows <- list()
  divergent_rows <- list()
  
  for (prot in names(protein_split)) {
    sub_df <- protein_split[[prot]]
    scores <- sub_df$score
    n_sig <- length(scores)
    
    if (n_sig < 2L) next
    
    n_pos <- sum(scores > 0)
    n_neg <- sum(scores < 0)
    
    mean_score <- mean(scores)
    min_score <- min(scores)
    max_score <- max(scores)
    range_score <- max_score - min_score
    log_factor <- log2(n_sig + 1)
    
    if (n_pos == 0L || n_neg == 0L) {
      convergence_score <- abs(mean_score) * log_factor
      
      convergent_rows[[length(convergent_rows) + 1L]] <- data.frame(
        summary_type = "top_convergent_proteins",
        direction = ifelse(mean_score >= 0, "positive_convergent", "negative_convergent"),
        rank = NA_integer_,
        sample = sample,
        pathway = NA_character_,
        protein = prot,
        score = convergence_score,
        conf_total = NA_integer_,
        padj = NA_real_,
        pval = NA_real_,
        n_signatures = n_sig,
        n_positive = n_pos,
        n_negative = n_neg,
        mean_score = mean_score,
        min_score = min_score,
        max_score = max_score,
        stringsAsFactors = FALSE
      )
    }
    
    if (n_pos > 0L && n_neg > 0L) {
      balance <- min(n_pos, n_neg) / n_sig
      divergence_score <- range_score * balance * log_factor
      
      divergent_rows[[length(divergent_rows) + 1L]] <- data.frame(
        summary_type = "top_divergent_proteins",
        direction = "mixed_sign",
        rank = NA_integer_,
        sample = sample,
        pathway = NA_character_,
        protein = prot,
        score = divergence_score,
        conf_total = NA_integer_,
        padj = NA_real_,
        pval = NA_real_,
        n_signatures = n_sig,
        n_positive = n_pos,
        n_negative = n_neg,
        mean_score = mean_score,
        min_score = min_score,
        max_score = max_score,
        stringsAsFactors = FALSE
      )
    }
  }
  
  conv_df <- .bind_summary_frames(convergent_rows)
  if (nrow(conv_df) > 0L) {
    conv_df <- conv_df[order(-conv_df$score, -abs(conv_df$mean_score)), , drop = FALSE]
    conv_df <- utils::head(conv_df, top_n)
    conv_df$rank <- seq_len(nrow(conv_df))
  }
  
  div_df <- .bind_summary_frames(divergent_rows)
  if (nrow(div_df) > 0L) {
    div_df <- div_df[order(-div_df$score, -div_df$n_signatures), , drop = FALSE]
    div_df <- utils::head(div_df, top_n)
    div_df$rank <- seq_len(nrow(div_df))
  }
  
  out <- .bind_summary_frames(list(conv_df, div_df))
  rownames(out) <- NULL
  out
}

.compute_protein_support_stats <- function(compass_res, protein_map, sample) {
  df <- data.frame(
    protein = unname(protein_map[rownames(compass_res)]),
    score = as.numeric(compass_res[, sample]),
    stringsAsFactors = FALSE
  )
  
  df <- df[!is.na(df$protein) & df$protein != "" & is.finite(df$score), , drop = FALSE]
  split_df <- split(df, df$protein)
  
  rows <- lapply(names(split_df), function(prot) {
    sub_df <- split_df[[prot]]
    scores <- sub_df$score
    n_sig <- length(scores)
    n_pos <- sum(scores > 0)
    n_neg <- sum(scores < 0)
    consistency_index <- abs(n_pos - n_neg) / n_sig
    support_weight <- 1 + 0.15 * log2(n_sig + 1)
    
    data.frame(
      protein = prot,
      n_signatures = n_sig,
      n_positive = n_pos,
      n_negative = n_neg,
      consistency_index = consistency_index,
      support_weight = support_weight,
      stringsAsFactors = FALSE
    )
  })
  
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

.compass_conf_weight <- function(conf_total) {
  ifelse(
    is.na(conf_total), 1.00,
    ifelse(conf_total >= 3, 1.40,
           ifelse(conf_total == 2, 1.20, 1.00)
    )
  )
}

.compute_hypothesis_shortlist <- function(compass_res,
                                          expr_mat,
                                          protein_map,
                                          sample,
                                          shortlist_n = 3L) {
  conf_total <- suppressWarnings(
    as.integer(sub(".*_c([0-9]+)$", "\\1", rownames(compass_res)))
  )
  
  activity_df <- data.frame(
    pathway = rownames(compass_res),
    protein = unname(protein_map[rownames(compass_res)]),
    score = as.numeric(compass_res[, sample]),
    conf_total = conf_total,
    stringsAsFactors = FALSE
  )
  
  activity_df <- activity_df[!is.na(activity_df$protein) & activity_df$protein != "" & is.finite(activity_df$score), , drop = FALSE]
  
  support_stats <- .compute_protein_support_stats(
    compass_res = compass_res,
    protein_map = protein_map,
    sample = sample
  )
  
  discordance_df <- .compute_rna_activity_discordance_full(
    compass_res = compass_res,
    expr_mat = expr_mat,
    protein_map = protein_map,
    sample = sample
  )
  
  .rank_category <- function(df, summary_type, direction_label) {
    if (nrow(df) == 0L) return(data.frame())
    
    idx <- match(df$protein, support_stats$protein)
    df$n_signatures <- support_stats$n_signatures[idx]
    df$n_positive <- support_stats$n_positive[idx]
    df$n_negative <- support_stats$n_negative[idx]
    df$consistency_index <- support_stats$consistency_index[idx]
    df$support_weight <- support_stats$support_weight[idx]
    
    df$conf_weight <- .compass_conf_weight(df$conf_total)
    df$consistency_weight <- 0.5 + 0.5 * df$consistency_index
    df$hypothesis_score <- abs(df$score) * df$conf_weight * df$consistency_weight * df$support_weight
    
    df <- df[order(-df$hypothesis_score, -abs(df$score), -df$conf_total), , drop = FALSE]
    df <- df[!duplicated(df$protein), , drop = FALSE]
    df <- utils::head(df, shortlist_n)
    
    if (nrow(df) == 0L) return(data.frame())
    
    df$summary_type <- summary_type
    df$direction <- direction_label
    df$rank <- seq_len(nrow(df))
    df$sample <- sample
    df$padj <- NA_real_
    df$pval <- NA_real_
    
    df[, c(
      "summary_type", "direction", "rank", "sample", "pathway", "protein",
      "score", "conf_total", "padj", "pval",
      "n_signatures", "n_positive", "n_negative",
      "consistency_index", "support_weight", "conf_weight", "hypothesis_score"
    ), drop = FALSE]
  }
  
  pos_act <- activity_df[activity_df$score > 0, , drop = FALSE]
  neg_act <- activity_df[activity_df$score < 0, , drop = FALSE]
  pos_dis <- discordance_df[discordance_df$score > 0, , drop = FALSE]
  neg_dis <- discordance_df[discordance_df$score < 0, , drop = FALSE]
  
  out <- .bind_summary_frames(list(
    .rank_category(pos_act, "shortlist_positive_activity", "positive"),
    .rank_category(neg_act, "shortlist_negative_activity", "negative"),
    .rank_category(pos_dis, "shortlist_activity_gt_rna", "activity_gt_rna"),
    .rank_category(neg_dis, "shortlist_rna_gt_activity", "rna_gt_activity")
  ))
  
  rownames(out) <- NULL
  out
}

.bind_summary_frames <- function(frames) {
  if (length(frames) == 0L) {
    return(data.frame())
  }
  
  is_empty <- vapply(frames, function(x) is.null(x) || nrow(x) == 0L, logical(1))
  frames <- frames[!is_empty]
  
  if (length(frames) == 0L) {
    return(data.frame())
  }
  
  all_cols <- unique(unlist(lapply(frames, names), use.names = FALSE))
  
  frames_aligned <- lapply(frames, function(df) {
    missing_cols <- setdiff(all_cols, names(df))
    if (length(missing_cols) > 0L) {
      for (mc in missing_cols) {
        df[[mc]] <- NA
      }
    }
    df[, all_cols, drop = FALSE]
  })
  
  out <- do.call(rbind, frames_aligned)
  rownames(out) <- NULL
  out
}

.resolve_protein_map <- function(pathways, index = NULL) {
  if (!is.null(index)) {
    required_cols <- c("pathway", "protein")
    missing_cols <- setdiff(required_cols, colnames(index))
    
    if (length(missing_cols) > 0L) {
      stop(
        "`index` is missing required columns: ",
        paste(missing_cols, collapse = ", "),
        call. = FALSE
      )
    }
    
    protein_map <- stats::setNames(as.character(index$protein), index$pathway)
    out <- protein_map[pathways]
    
    missing <- is.na(out) | out == ""
    if (any(missing)) {
      out[missing] <- .derive_protein_from_pathway(pathways[missing])
    }
    
    return(out)
  }
  
  out <- .derive_protein_from_pathway(pathways)
  names(out) <- pathways
  out
}

.derive_protein_from_pathway <- function(pathways) {
  ifelse(
    grepl("_XPR", pathways),
    sub("_XPR.*$", "", pathways),
    sub("_c[0-9]+$", "", pathways)
  )
}

.compass_build_overview_info <- function(compass_res, summary_df, sample = NULL) {
  is_fgsea_like <- is.data.frame(compass_res) &&
    all(c("pathway", "NES") %in% colnames(compass_res))
  
  if (is_fgsea_like) {
    score_vec <- compass_res$NES
    
    if ("conf_total" %in% colnames(compass_res)) {
      conf_all <- compass_res$conf_total
    } else {
      conf_all <- suppressWarnings(
        as.integer(sub(".*_c([0-9]+)$", "\\1", compass_res$pathway))
      )
    }
  } else {
    compass_mat <- as.matrix(compass_res)
    
    if (!is.null(sample)) {
      if (!sample %in% colnames(compass_mat)) {
        stop("Requested `sample` not found in `compass_res` columns.", call. = FALSE)
      }
      score_vec <- as.numeric(compass_mat[, sample])
    } else {
      score_vec <- as.numeric(compass_mat)
    }
    
    conf_all <- suppressWarnings(
      as.integer(sub(".*_c([0-9]+)$", "\\1", rownames(compass_mat)))
    )
  }
  
  score_vec <- score_vec[is.finite(score_vec) & !is.na(score_vec)]
  
  if (length(score_vec) == 0L) {
    stats_list <- list(
      median_score = NA_real_,
      mean_score = NA_real_,
      frac_positive = NA_real_,
      frac_negative = NA_real_,
      strongest_positive = NA_real_,
      strongest_negative = NA_real_,
      overall_shift = "not available"
    )
  } else {
    stats_list <- list(
      median_score = stats::median(score_vec),
      mean_score = mean(score_vec),
      frac_positive = mean(score_vec > 0),
      frac_negative = mean(score_vec < 0),
      strongest_positive = max(score_vec),
      strongest_negative = min(score_vec),
      overall_shift = .compass_describe_shift(stats::median(score_vec))
    )
  }
  
  section_conf <- lapply(split(summary_df$conf_total, summary_df$summary_type), .compass_extract_conf_counts)
  
  list(
    stats = stats_list,
    conf_all = .compass_extract_conf_counts(conf_all),
    conf_by_section = section_conf
  )
}

.compass_extract_conf_counts <- function(x) {
  x <- x[!is.na(x)]
  
  if (length(x) == 0L) {
    return(setNames(integer(), character()))
  }
  
  tbl <- table(x)
  ord <- order(as.numeric(names(tbl)))
  tbl <- tbl[ord]
  
  setNames(as.integer(tbl), paste0("c", names(tbl)))
}

.compass_describe_shift <- function(median_score) {
  if (is.na(median_score)) {
    return("not available")
  }
  
  if (median_score > 0.05) {
    return("positive-shifted")
  }
  
  if (median_score < -0.05) {
    return("negative-shifted")
  }
  
  "near balanced"
}

.format_conf_counts <- function(counts) {
  if (length(counts) == 0L) {
    return("none")
  }
  
  paste(paste0(names(counts), ": ", counts), collapse = " | ")
}

.print_overview_block <- function(overview_info) {
  s <- overview_info$stats
  
  cat("\nOverview\n")
  cat("- median score: ", formatC(s$median_score, format = "f", digits = 3), "\n", sep = "")
  cat("- mean score:   ", formatC(s$mean_score, format = "f", digits = 3), "\n", sep = "")
  cat("- positive:     ", formatC(100 * s$frac_positive, format = "f", digits = 1), "%\n", sep = "")
  cat("- negative:     ", formatC(100 * s$frac_negative, format = "f", digits = 1), "%\n", sep = "")
  cat("- strongest +:  ", formatC(s$strongest_positive, format = "f", digits = 3), "\n", sep = "")
  cat("- strongest -:  ", formatC(s$strongest_negative, format = "f", digits = 3), "\n", sep = "")
  cat("- overall shift: ", s$overall_shift, "\n", sep = "")
}

.print_confidence_block <- function(title, counts) {
  cat("\n", title, "\n", sep = "")
  cat("- ", .format_conf_counts(counts), "\n", sep = "")
}

.print_compass_summary <- function(summary_df,
                                   top_n,
                                   print_n = 5L,
                                   sample = NULL,
                                   shortlist_n = 3L,
                                   overview_info = NULL) {
  cat("\n===========================\n")
  cat("      COMPASS summary      \n")
  cat("===========================\n")
  
  if (!is.null(sample)) {
    cat("Sample: ", sample, "\n", sep = "")
  }
  
  cat("Returned rows: ", nrow(summary_df), "\n", sep = "")
  cat("top_n: ", top_n, "\n", sep = "")
  cat("print_n: ", print_n, "\n", sep = "")
  
  if (!is.null(overview_info)) {
    .print_overview_block(overview_info)
    .print_confidence_block(
      title = "Confidence overview (all signatures)",
      counts = overview_info$conf_all
    )
  }
  
  type_order <- c(
    "top_positive_activity",
    "top_negative_activity",
    "top_activity_gt_rna",
    "top_rna_gt_activity",
    "top_convergent_proteins",
    "top_divergent_proteins",
    "top_variable_activity",
    "top_positive_nes",
    "top_negative_nes"
  )
  
  present_types <- intersect(type_order, unique(summary_df$summary_type))
  
  for (st in present_types) {
    section_df <- summary_df[summary_df$summary_type == st, , drop = FALSE]
    if (nrow(section_df) == 0L) next
    
    title <- switch(
      st,
      top_positive_activity   = "[1] Top positive activity",
      top_negative_activity   = "[2] Top negative activity",
      top_activity_gt_rna     = "[3] Top activity > RNA discordance",
      top_rna_gt_activity     = "[4] Top RNA > activity discordance",
      top_convergent_proteins = "[5] Top convergent proteins",
      top_divergent_proteins  = "[6] Top divergent proteins (warning)",
      top_variable_activity   = "[7] Top variable activity",
      top_positive_nes        = "[1] Top positive NES",
      top_negative_nes        = "[2] Top negative NES",
      st
    )
    
    cat("\n", title, "\n", sep = "")
    
    if (!is.null(overview_info) &&
        !is.null(overview_info$conf_by_section[[st]]) &&
        length(overview_info$conf_by_section[[st]]) > 0L) {
      cat("Conf: ", .format_conf_counts(overview_info$conf_by_section[[st]]), "\n", sep = "")
    }
    
    .print_summary_section(
      df = section_df,
      print_n = print_n,
      sample = sample
    )
  }
  
  shortlist_types <- c(
    "shortlist_positive_activity",
    "shortlist_negative_activity",
    "shortlist_activity_gt_rna",
    "shortlist_rna_gt_activity"
  )
  
  if (any(shortlist_types %in% unique(summary_df$summary_type))) {
    .print_hypothesis_matrix(
      summary_df = summary_df,
      shortlist_n = shortlist_n
    )
  }
  
  invisible(NULL)
}

.print_summary_section <- function(df, print_n = 5L, sample = NULL) {
  show_df <- utils::head(df, min(print_n, nrow(df)))
  section_type <- unique(show_df$summary_type)[1]
  
  if (section_type %in% c("top_convergent_proteins", "top_divergent_proteins")) {
    disp <- data.frame(
      rank = show_df$rank,
      protein = show_df$protein,
      score = formatC(show_df$score, format = "f", digits = 3),
      n_sig = show_df$n_signatures,
      n_pos = show_df$n_positive,
      n_neg = show_df$n_negative,
      mean = formatC(show_df$mean_score, format = "f", digits = 3),
      min = formatC(show_df$min_score, format = "f", digits = 3),
      max = formatC(show_df$max_score, format = "f", digits = 3),
      stringsAsFactors = FALSE
    )
    
    print(disp, row.names = FALSE, right = FALSE)
    return(invisible(NULL))
  }
  
  disp <- data.frame(
    rank = show_df$rank,
    protein = show_df$protein,
    score = formatC(show_df$score, format = "f", digits = 3),
    conf = ifelse(is.na(show_df$conf_total), NA_character_, paste0("c", show_df$conf_total)),
    stringsAsFactors = FALSE
  )
  
  if ("padj" %in% colnames(show_df) && !all(is.na(show_df$padj))) {
    disp$padj <- formatC(show_df$padj, format = "e", digits = 2)
  }
  
  if ("pathway" %in% colnames(show_df) && !all(is.na(show_df$pathway))) {
    disp$pathway <- show_df$pathway
  }
  
  if (is.null(sample) && "sample" %in% colnames(show_df) && !all(is.na(show_df$sample))) {
    disp$sample <- show_df$sample
    
    if ("padj" %in% colnames(disp) && "pathway" %in% colnames(disp)) {
      disp <- disp[, c("rank", "protein", "score", "padj", "conf", "sample", "pathway"), drop = FALSE]
    } else if ("pathway" %in% colnames(disp)) {
      disp <- disp[, c("rank", "protein", "score", "conf", "sample", "pathway"), drop = FALSE]
    } else {
      disp <- disp[, c("rank", "protein", "score", "conf", "sample"), drop = FALSE]
    }
  } else {
    if ("padj" %in% colnames(disp) && "pathway" %in% colnames(disp)) {
      disp <- disp[, c("rank", "protein", "score", "padj", "conf", "pathway"), drop = FALSE]
    } else if ("pathway" %in% colnames(disp)) {
      disp <- disp[, c("rank", "protein", "score", "conf", "pathway"), drop = FALSE]
    } else {
      disp <- disp[, c("rank", "protein", "score", "conf"), drop = FALSE]
    }
  }
  
  print(disp, row.names = FALSE, right = FALSE)
}

.pad_or_trim <- function(x, width) {
  x <- ifelse(is.na(x), "", x)
  x <- enc2utf8(as.character(x))
  x <- ifelse(nchar(x, type = "width") > width, substr(x, 1, width), x)
  sprintf(paste0("%-", width, "s"), x)
}

.make_hypothesis_cell <- function(df, title, shortlist_n = 3L, width = 40L) {
  lines <- c(title)
  
  if (nrow(df) == 0L) {
    lines <- c(lines, "(none)")
  } else {
    df <- utils::head(df, shortlist_n)
    
    for (i in seq_len(nrow(df))) {
      prot <- substr(df$protein[i], 1, 10)
      conf <- ifelse(is.na(df$conf_total[i]), "-", paste0("c", df$conf_total[i]))
      n_sig <- ifelse(is.na(df$n_signatures[i]), 1L, df$n_signatures[i])
      hsc <- formatC(df$hypothesis_score[i], format = "f", digits = 2)
      
      line <- sprintf(
        "%d %-10s h=%s %s n=%d",
        i, prot, hsc, conf, n_sig
      )
      
      lines <- c(lines, line)
    }
  }
  
  target_len <- shortlist_n + 1L
  while (length(lines) < target_len) {
    lines <- c(lines, "")
  }
  
  vapply(lines, .pad_or_trim, character(1), width = width)
}

.print_hypothesis_matrix <- function(summary_df, shortlist_n = 3L) {
  cat("\n========================================\n")
  cat(" Hypothesis generation shortlist\n")
  cat(" h = |signal| x conf x consistency x support\n")
  cat(" conf weights: c1=1.00, c2=1.20, c3=1.40\n")
  cat("========================================\n")
  
  width <- 40L
  
  cell_tl <- .make_hypothesis_cell(
    summary_df[summary_df$summary_type == "shortlist_positive_activity", , drop = FALSE],
    title = "POS activity",
    shortlist_n = shortlist_n,
    width = width
  )
  
  cell_tr <- .make_hypothesis_cell(
    summary_df[summary_df$summary_type == "shortlist_negative_activity", , drop = FALSE],
    title = "NEG activity",
    shortlist_n = shortlist_n,
    width = width
  )
  
  cell_bl <- .make_hypothesis_cell(
    summary_df[summary_df$summary_type == "shortlist_activity_gt_rna", , drop = FALSE],
    title = "Act > RNA",
    shortlist_n = shortlist_n,
    width = width
  )
  
  cell_br <- .make_hypothesis_cell(
    summary_df[summary_df$summary_type == "shortlist_rna_gt_activity", , drop = FALSE],
    title = "RNA > Act",
    shortlist_n = shortlist_n,
    width = width
  )
  
  top_border <- paste0("┌", strrep("─", width), "┬", strrep("─", width), "┐")
  mid_border <- paste0("├", strrep("─", width), "┼", strrep("─", width), "┤")
  bot_border <- paste0("└", strrep("─", width), "┴", strrep("─", width), "┘")
  
  cat(top_border, "\n", sep = "")
  for (i in seq_along(cell_tl)) {
    cat("│", cell_tl[i], "│", cell_tr[i], "│\n", sep = "")
  }
  cat(mid_border, "\n", sep = "")
  for (i in seq_along(cell_bl)) {
    cat("│", cell_bl[i], "│", cell_br[i], "│\n", sep = "")
  }
  cat(bot_border, "\n", sep = "")
}