# Internal query-input and gene-set-overlap utilities for vmpaR.

.vmpa_validate_gene_ids <- function(gene_ids, location) {
  if (is.null(gene_ids)) {
    stop("`", location, "` must contain gene symbols.", call. = FALSE)
  }

  if (anyNA(gene_ids)) {
    stop("`", location, "` contains missing gene symbols.", call. = FALSE)
  }

  trimmed_ids <- trimws(gene_ids)

  if (any(!nzchar(trimmed_ids))) {
    stop("`", location, "` contains empty or whitespace-only gene symbols.", call. = FALSE)
  }

  if (any(gene_ids != trimmed_ids)) {
    stop(
      "`", location, "` contains gene symbols with leading or trailing whitespace.",
      call. = FALSE
    )
  }

  if (anyDuplicated(gene_ids) > 0L) {
    duplicate_ids <- unique(gene_ids[duplicated(gene_ids)])
    stop(
      "`", location, "` contains duplicated gene symbols. Examples: ",
      paste(utils::head(duplicate_ids, 10L), collapse = ", "),
      call. = FALSE
    )
  }

  invisible(gene_ids)
}

.vmpa_prepare_gsva_input <- function(expr_mat) {
  if (!(is.matrix(expr_mat) || is.data.frame(expr_mat))) {
    stop(
      "`algorithm = \"gsva\"` requires `input` to be a matrix or data.frame ",
      "with genes in rows and samples in columns.",
      call. = FALSE
    )
  }

  expr_mat <- as.matrix(expr_mat)

  if (!is.numeric(expr_mat)) {
    stop("`input` must contain numeric expression values.", call. = FALSE)
  }

  .vmpa_validate_gene_ids(rownames(expr_mat), "rownames(input)")

  if (is.null(colnames(expr_mat))) {
    stop("`input` must have sample names as colnames.", call. = FALSE)
  }

  if (anyNA(expr_mat) || any(!is.finite(expr_mat))) {
    stop("`input` contains NA, NaN, or infinite values.", call. = FALSE)
  }

  expr_mat
}

.vmpa_prepare_fgsea_input <- function(stats_vec) {
  if (!is.numeric(stats_vec)) {
    stop(
      "`algorithm = \"fgsea\"` requires `input` to be a named numeric vector.",
      call. = FALSE
    )
  }

  .vmpa_validate_gene_ids(names(stats_vec), "names(input)")

  stats_vec <- stats_vec[is.finite(stats_vec) & !is.na(stats_vec)]

  if (length(stats_vec) == 0L) {
    stop("`input` contains no finite numeric values after filtering.", call. = FALSE)
  }

  sort(stats_vec, decreasing = TRUE)
}

.vmpa_validate_min_size <- function(min_size, argument) {
  if (!is.numeric(min_size) || length(min_size) != 1L || is.na(min_size) ||
      !is.finite(min_size) || min_size <= 0 || min_size != as.integer(min_size)) {
    stop("`", argument, "` must be a single positive integer.", call. = FALSE)
  }

  as.integer(min_size)
}

.vmpa_filter_gene_sets_by_overlap <- function(gene_set_list,
                                                input_genes,
                                                min_size,
                                                backend,
                                                verbose = TRUE) {
  min_size <- .vmpa_validate_min_size(min_size, paste0(backend, "_min_size"))
  .vmpa_validate_gene_ids(input_genes, "input gene identifiers")

  if (!is.list(gene_set_list) || length(gene_set_list) == 0L) {
    stop("`gene_set_list` must be a non-empty list.", call. = FALSE)
  }

  overlap_sizes <- vapply(
    gene_set_list,
    function(gene_set) length(intersect(input_genes, unique(gene_set))),
    integer(1)
  )

  reference_genes <- unique(unlist(gene_set_list, use.names = FALSE))
  total_overlap <- length(intersect(input_genes, reference_genes))

  if (total_overlap == 0L) {
    stop(
      "No input gene symbols match the selected VMPA gene sets. ",
      "VMPA expects case-sensitive human gene symbols matching the VMPA/LINCS reference; ",
      "Ensembl and Entrez identifiers are not converted automatically.",
      call. = FALSE
    )
  }

  keep <- overlap_sizes >= min_size
  total_sets <- length(gene_set_list)
  retained_sets <- sum(keep)

  .vmpa_msg(
    verbose,
    total_overlap, " unique input genes match the selected VMPA gene sets."
  )
  .vmpa_msg(
    verbose,
    retained_sets, " of ", total_sets, " selected VMPA gene sets have at least ",
    min_size, if (min_size == 1L) " represented gene" else " represented genes",
    "; ", total_sets - retained_sets, " excluded."
  )

  if (retained_sets == 0L) {
    stop(
      "No selected VMPA gene set has at least ", min_size,
      if (min_size == 1L) " represented input gene." else " represented input genes.",
      call. = FALSE
    )
  }

  filtered <- gene_set_list[keep]
  kept_overlap_sizes <- overlap_sizes[keep]
  signature_metadata <- attr(gene_set_list, "signature_metadata")

  if (is.data.frame(signature_metadata) &&
      "gene_set_name" %in% colnames(signature_metadata)) {
    metadata_index <- match(names(filtered), signature_metadata$gene_set_name)
    signature_metadata <- signature_metadata[metadata_index, , drop = FALSE]
    signature_metadata$input_overlap <- unname(kept_overlap_sizes)
    rownames(signature_metadata) <- signature_metadata$gene_set_name
    attr(filtered, "signature_metadata") <- signature_metadata
  }

  unique_selection <- attr(gene_set_list, "unique_selection")
  if (is.data.frame(unique_selection) &&
      "gene_set_name" %in% colnames(unique_selection)) {
    unique_selection <- unique_selection[
      unique_selection$gene_set_name %in% names(filtered),
      ,
      drop = FALSE
    ]
  }
  if (!is.null(unique_selection)) {
    attr(filtered, "unique_selection") <- unique_selection
  }

  attr(filtered, "input_overlap") <- kept_overlap_sizes
  attr(filtered, "total_input_overlap") <- total_overlap
  filtered
}

.vmpa_subset_gene_sets_to_backend <- function(gene_set_list, used_names) {
  used_names <- unique(as.character(used_names))

  if (length(used_names) == 0L) {
    stop("The selected backend did not use any VMPA gene sets.", call. = FALSE)
  }

  missing_names <- setdiff(used_names, names(gene_set_list))
  if (length(missing_names) > 0L) {
    stop("Internal error: backend returned unknown VMPA gene-set names.", call. = FALSE)
  }

  filtered <- gene_set_list[used_names]

  signature_metadata <- attr(gene_set_list, "signature_metadata")
  if (is.data.frame(signature_metadata) &&
      "gene_set_name" %in% colnames(signature_metadata)) {
    metadata_index <- match(used_names, signature_metadata$gene_set_name)
    signature_metadata <- signature_metadata[metadata_index, , drop = FALSE]
    rownames(signature_metadata) <- signature_metadata$gene_set_name
    attr(filtered, "signature_metadata") <- signature_metadata
  }

  overlap_sizes <- attr(gene_set_list, "input_overlap")
  if (!is.null(overlap_sizes)) {
    attr(filtered, "input_overlap") <- overlap_sizes[used_names]
  }

  attr(filtered, "total_input_overlap") <- attr(gene_set_list, "total_input_overlap")
  unique_selection <- attr(gene_set_list, "unique_selection")
  if (is.data.frame(unique_selection) &&
      "gene_set_name" %in% colnames(unique_selection)) {
    unique_selection <- unique_selection[
      unique_selection$gene_set_name %in% used_names,
      ,
      drop = FALSE
    ]
  }
  attr(filtered, "unique_selection") <- unique_selection
  filtered
}
