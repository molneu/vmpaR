# Internal gene-set building utilities for vmpaR
# These functions are not user-facing and should not be exported.

.vmpa_read_subset_gct <- function(subset_file) {
  if (!is.character(subset_file) || length(subset_file) != 1L || is.na(subset_file) || subset_file == "") {
    stop("`subset_file` must be a single non-empty character string.", call. = FALSE)
  }

  if (!file.exists(subset_file)) {
    stop("Subset file not found: ", subset_file, call. = FALSE)
  }

  if (!requireNamespace("cmapR", quietly = TRUE)) {
    stop("Package `cmapR` must be installed.", call. = FALSE)
  }

  gct <- tryCatch(
    readRDS(subset_file),
    error = function(e) {
      stop("Failed to read subset file: ", conditionMessage(e), call. = FALSE)
    }
  )

  if (!methods::is(gct, "GCT")) {
    stop("The loaded object is not a `GCT` object.", call. = FALSE)
  }

  required_rdesc_cols <- c("symbol")
  required_cdesc_cols <- c("id", "cmap_name", "cps_conf_total", "cancer_driver_summary")

  missing_rdesc <- setdiff(required_rdesc_cols, colnames(gct@rdesc))
  missing_cdesc <- setdiff(required_cdesc_cols, colnames(gct@cdesc))

  if (length(missing_rdesc) > 0L) {
    stop(
      "Missing required row metadata columns: ",
      paste(missing_rdesc, collapse = ", "),
      call. = FALSE
    )
  }

  if (length(missing_cdesc) > 0L) {
    stop(
      "Missing required column metadata columns: ",
      paste(missing_cdesc, collapse = ", "),
      call. = FALSE
    )
  }

  if (nrow(gct@mat) != nrow(gct@rdesc)) {
    stop(
      "Row metadata and expression matrix row count do not match.",
      call. = FALSE
    )
  }

  if (ncol(gct@mat) != nrow(gct@cdesc)) {
    stop(
      "Column metadata and expression matrix column count do not match.",
      call. = FALSE
    )
  }

  gct
}

.vmpa_build_gene_sets <- function(gct,
                                     n = 250L,
                                     min_conf = 1L,
                                     targets = NULL,
                                     driver_filter = FALSE) {
  if (!methods::is(gct, "GCT")) {
    stop("`gct` must be a `GCT` object.", call. = FALSE)
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

  keep_idx <- which(gct@cdesc$cps_conf_total >= min_conf)

  if (driver_filter) {
    # Broad filter: retain any signature with non-"None" cancer_driver_summary.
    # This does not restrict to canonical drivers or specific driver roles.
    cds <- gct@cdesc$cancer_driver_summary
    keep_idx <- keep_idx[cds[keep_idx] != "None"]
  }

  if (!is.null(targets)) {
    keep_idx <- keep_idx[gct@cdesc$cmap_name[keep_idx] %in% targets]
  }

  if (length(keep_idx) == 0L) {
    out <- list()
    attr(out, "signature_metadata") <- data.frame()
    return(out)
  }

  mat <- gct@mat[, keep_idx, drop = FALSE]
  rownames(mat) <- as.character(gct@rdesc$symbol)

  sig_meta <- gct@cdesc[keep_idx, , drop = FALSE]
  gene_set_names <- paste0(sig_meta$id, "_c", sig_meta$cps_conf_total)
  sig_meta$gene_set_name <- gene_set_names

  if (anyDuplicated(gene_set_names) > 0L) {
    dup_names <- unique(gene_set_names[duplicated(gene_set_names)])
    stop(
      "Duplicate gene-set names detected. Examples: ",
      paste(utils::head(dup_names, 10), collapse = ", "),
      call. = FALSE
    )
  }

  # `n` refers to the final number of unique genes per set.
  # For duplicated gene symbols, the most downregulated occurrence is retained
  # because duplicates are removed after ordering and before truncation to `n`.

  gene_set_list <- lapply(seq_len(ncol(mat)), function(j) {
    vals <- as.numeric(mat[, j])
    genes <- trimws(rownames(mat))
    keep <- is.finite(vals) & !is.na(genes) & nzchar(genes)

    ord <- order(vals[keep], decreasing = FALSE)
    ordered_genes <- genes[keep][ord]
    ordered_unique_genes <- unique(ordered_genes)

    utils::head(ordered_unique_genes, n)
  })

  names(gene_set_list) <- gene_set_names

  keep_nonempty <- lengths(gene_set_list) > 0L
  gene_set_list <- gene_set_list[keep_nonempty]
  sig_meta <- sig_meta[keep_nonempty, , drop = FALSE]
  rownames(sig_meta) <- sig_meta$gene_set_name

  if (!identical(sig_meta$gene_set_name, names(gene_set_list))) {
    stop(
      "Internal error: signature metadata are not aligned with gene_set_list.",
      call. = FALSE
    )
  }

  attr(gene_set_list, "signature_metadata") <- sig_meta

  gene_set_list
}
