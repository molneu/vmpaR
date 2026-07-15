#' Run VMPA
#'
#' `vmpa()` applies context-matched VMPA reference signatures to a user
#' query using either a GSVA-based or FGSEA-based workflow.
#'
#' With `algorithm = "gsva"`, the input is a gene-by-sample expression matrix.
#' VMPA returns a `vmpa_result` data frame containing VMPA metadata
#' and one score column per sample.
#'
#' With `algorithm = "fgsea"`, the input is a named numeric vector of gene-level
#' statistics or rankings. VMPA returns a `vmpa_result` data frame
#' containing FGSEA enrichment statistics together with VMPA metadata.
#'
#' @details
#' The selected cancer context determines which VMPA reference signatures are
#' used. Required reference data are resolved automatically by vmpaR.
#'
#' VMPA reference data may contain multiple perturbation signatures for the
#' same target/protein.
#'
#' With `unique = FALSE`, all signatures passing the selected filters are kept
#' as separate features.
#'
#' With `unique = TRUE`, signatures are grouped by their target/protein
#' annotation (`cmap_name`). If validation metadata are available, validated
#' signatures are prioritized. Among the remaining candidates, signatures with
#' the highest `cps_conf_total` are retained.
#'
#' For `algorithm = "gsva"`, all equally prioritized signatures for the same
#' target/protein are scored separately. If `gsva_score_scaling` is not
#' `"none"`, scaling is applied before target-level averaging. The final output
#' contains one row per target/protein.
#'
#' For `algorithm = "fgsea"`, retained signatures are tested separately.
#' If multiple equally prioritized signatures remain for the same target/protein,
#' they are returned as separate FGSEA rows rather than being averaged or
#' arbitrarily collapsed.
#'
#' The `seed` is applied locally around the complete FGSEA backend call. The
#' user's global random-number state is restored afterwards, including when the
#' backend raises an error. Set `seed = NULL` to disable package-controlled
#' seeding. The parameter has no effect when `algorithm = "gsva"`.
#'
#' Input gene identifiers must be unique, case-sensitive human gene symbols
#' matching those used by the VMPA/LINCS reference (for example `TP53`, `EGFR`,
#' or `AKT1`). Ensembl and Entrez identifiers are not converted automatically.
#' Missing, blank, whitespace-padded, or duplicated input symbols are rejected.
#'
#' Overlap is evaluated separately for every selected VMPA gene set after
#' non-finite FGSEA statistics have been removed. GSVA retains gene sets with at
#' least `gsva_min_size` represented genes; FGSEA uses `fgsea_min_size`.
#' If no input symbol matches or no gene set reaches the relevant threshold,
#' VMPA stops before calling the backend. With `verbose = TRUE`, VMPA reports
#' the global overlap, threshold exclusions, and the number of sets actually
#' used by the backend.
#'
#' Duplicate symbols in the VMPA reference are handled during gene-set
#' construction: the occurrence with the lowest reference-signature value is
#' retained, and a symbol counts at most once toward `n`.
#'
#' @param input Query input.
#'   For `algorithm = "gsva"`, provide a numeric matrix or data frame with genes
#'   in rows and samples in columns. Unique human gene symbols must be stored in
#'   row names, and sample identifiers must be stored in column names.
#'   For `algorithm = "fgsea"`, provide a named numeric vector of gene-level
#'   statistics or rankings. Names must contain unique human gene symbols.
#'   Gene-symbol matching is case-sensitive. Ensembl and Entrez identifiers are
#'   not converted automatically.
#' @param context Character scalar. Cancer context to use. One of:
#'   `"glioma"`, `"melanoma"`, `"nsclc"`, `"gastric"`, `"ovarian"`,
#'   `"crc"`, `"breast"`, `"prostate"`, `"pdac"`, or `"headneck"`.
#' @param algorithm Character scalar. Analysis workflow. Either `"gsva"` or
#'   `"fgsea"`. Default: `"gsva"`.
#' @param gsva_score_scaling Character scalar. Optional post-processing applied
#'   to GSVA-derived VMPA scores. One of `"none"`, `"sample_z"`,
#'   `"sample_pop_sd"`, or `"signature_z"`. Only used when
#'   `algorithm = "gsva"`. Default: `"none"`.
#' @param gsva_min_size Integer scalar. Minimum number of input genes represented
#'   in a VMPA gene set for the GSVA backend. Default: `1L`, preserving the GSVA
#'   default and the original VMPA GSVA workflow.
#' @param fgsea_min_size Integer scalar. Minimum number of finite ranked input
#'   genes represented in a VMPA gene set for the FGSEA backend. Default: `10L`,
#'   preserving the original VMPA FGSEA workflow.
#' @param unique Logical scalar. If `TRUE`, reduce results to one
#'   target/protein-level output per target. If `FALSE`, keep all VMPA
#'   signatures separately. Default: `TRUE`.
#' @param n Integer scalar. Number of bottom-ranked unique genes to include per
#'   VMPA signature. Default: `250L`.
#' @param min_conf Integer scalar. Minimum confidence score required for a
#'   VMPA signature to be included. Must be one of `1`, `2`, or `3`.
#'   Default: `1L`.
#' @param targets Optional character vector. If provided, only signatures for
#'   matching targets are retained. Default: `NULL`.
#' @param driver_filter Logical scalar. If `TRUE`, retain only signatures with
#'   a non-`"None"` cancer-driver annotation. This is a broad annotation filter,
#'   not a strict canonical-driver filter. Default: `FALSE`.
#' @param return_gene_sets Logical scalar. If `TRUE`, return a list containing
#'   the VMPA result, the gene sets used for the analysis, and selection
#'   metadata. Default: `FALSE`.
#' @param seed `NULL` or a single integer. Random seed applied locally to the
#'   complete FGSEA backend call. Set to `NULL` to use the current global random
#'   state without package-controlled seeding. Not used by GSVA. Default:
#'   `123L`.
#' @param verbose Logical scalar. If `TRUE`, print progress messages. Default:
#'   `TRUE`.
#'
#' @return
#' If `return_gene_sets = FALSE`, a `vmpa_result` data frame.
#'
#' For `algorithm = "gsva"`, rows correspond to targets/proteins when
#' `unique = TRUE` and to individual VMPA signatures when `unique = FALSE`.
#' Score columns correspond to samples.
#'
#' For `algorithm = "fgsea"`, rows correspond to tested VMPA signatures and
#' columns contain FGSEA enrichment statistics together with VMPA metadata.
#'
#' If `return_gene_sets = TRUE`, a list with:
#'
#' - `vmpa_result`: the `vmpa_result` data frame.
#' - `gene_sets`: gene sets used for the analysis.
#' - `signature_metadata`: metadata for the gene sets used.
#' - `context`: selected cancer context.
#' - `algorithm`: selected analysis workflow.
#' - `unique`: whether signature reduction was used.
#' - `unique_selection`: metadata describing selected signatures when
#'   `unique = TRUE`; otherwise `NULL`.
#'
#' @examples
#' \dontrun{
#' if (requireNamespace("Biobase", quietly = TRUE)) {
#'   data(kebir_gb, package = "vmpaR")
#'
#'   # GSVA workflow
#'   expr_mat <- Biobase::exprs(kebir_gb)
#'
#'   gsva_result <- vmpa(
#'     input = expr_mat,
#'     context = "glioma",
#'     algorithm = "gsva"
#'   )
#'
#'   # FGSEA workflow
#'   sample_metadata <- Biobase::pData(kebir_gb)
#'
#'   relapsed <- sample_metadata$relapse_TYPE != "n"
#'   naive <- sample_metadata$relapse_TYPE == "n"
#'
#'   stats_vec <- rowMeans(expr_mat[, relapsed, drop = FALSE]) -
#'     rowMeans(expr_mat[, naive, drop = FALSE])
#'
#'   stats_vec <- stats_vec[is.finite(stats_vec) & !is.na(stats_vec)]
#'   stats_vec <- sort(stats_vec, decreasing = TRUE)
#'
#'   fgsea_result <- vmpa(
#'     input = stats_vec,
#'     context = "glioma",
#'     algorithm = "fgsea"
#'   )
#' }
#' }
#'
#' @export
vmpa <- function(input,
                    context,
                    algorithm = c("gsva", "fgsea"),
                    gsva_score_scaling = c(
                      "none",
                      "sample_z",
                      "sample_pop_sd",
                      "signature_z"
                    ),
                    unique = TRUE,
                    n = 250L,
                    min_conf = 1L,
                    targets = NULL,
                    driver_filter = FALSE,
                    return_gene_sets = FALSE,
                    verbose = TRUE,
                    seed = 123L,
                    gsva_min_size = 1L,
                    fgsea_min_size = 10L) {
  algorithm <- match.arg(algorithm)
  gsva_score_scaling <- match.arg(gsva_score_scaling)
  context <- .vmpa_validate_context(context)
  seed <- .vmpa_validate_seed(seed)
  gsva_min_size <- .vmpa_validate_min_size(gsva_min_size, "gsva_min_size")
  fgsea_min_size <- .vmpa_validate_min_size(fgsea_min_size, "fgsea_min_size")

  input <- if (algorithm == "gsva") {
    .vmpa_prepare_gsva_input(input)
  } else {
    .vmpa_prepare_fgsea_input(input)
  }

  if (!is.logical(unique) || length(unique) != 1L || is.na(unique)) {
    stop("`unique` must be TRUE or FALSE.", call. = FALSE)
  }

  if (!is.logical(return_gene_sets) || length(return_gene_sets) != 1L || is.na(return_gene_sets)) {
    stop("`return_gene_sets` must be TRUE or FALSE.", call. = FALSE)
  }

  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("`verbose` must be TRUE or FALSE.", call. = FALSE)
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

  if (algorithm != "gsva" && gsva_score_scaling != "none") {
    stop(
      "`gsva_score_scaling` is only available when `algorithm = \"gsva\"`.",
      call. = FALSE
    )
  }

  # 1) Resolve subset file from cache or download it automatically
  subset_file <- .vmpa_resolve_subset_file(
    context = context,
    verbose = verbose
  )

  # 2) Read subset GCT and build VMPA gene sets internally
  gct <- .vmpa_read_subset_gct(subset_file)

  gene_set_list <- .vmpa_build_gene_sets(
    gct = gct,
    n = n,
    min_conf = min_conf,
    targets = targets,
    driver_filter = driver_filter
  )

  if (length(gene_set_list) == 0L) {
    stop("No VMPA gene sets available after filtering.", call. = FALSE)
  }

  # 3) Optionally prioritize signatures per target/protein
  if (isTRUE(unique)) {
    gene_set_list <- .vmpa_select_unique_candidates(
      gene_set_list = gene_set_list
    )
  }

  if (length(gene_set_list) == 0L) {
    stop("No VMPA gene sets available after unique reduction.", call. = FALSE)
  }

  # Metadata for the gene sets that are actually used for scoring
  signature_metadata <- attr(gene_set_list, "signature_metadata")

  backend_min_size <- if (algorithm == "gsva") gsva_min_size else fgsea_min_size
  input_genes <- if (algorithm == "gsva") rownames(input) else names(input)
  gene_set_list <- .vmpa_filter_gene_sets_by_overlap(
    gene_set_list = gene_set_list,
    input_genes = input_genes,
    min_size = backend_min_size,
    backend = algorithm,
    verbose = verbose
  )
  signature_metadata <- attr(gene_set_list, "signature_metadata")

  # 4) Run analysis with the selected algorithm
  if (algorithm == "gsva") {
    vmpa_result <- .vmpa_run_gsva(
      expr_mat = input,
      gene_set_list = gene_set_list,
      min_size = gsva_min_size,
      score_scaling = gsva_score_scaling,
      verbose = verbose
    )

    .vmpa_msg(
      verbose,
      nrow(vmpa_result), " VMPA gene sets were used by the GSVA backend."
    )

    gene_set_list <- .vmpa_subset_gene_sets_to_backend(
      gene_set_list,
      rownames(vmpa_result)
    )
    signature_metadata <- attr(gene_set_list, "signature_metadata")

    if (isTRUE(unique)) {
      vmpa_result <- .vmpa_reduce_unique_gsva_scores(
        score_mat = vmpa_result,
        signature_metadata = signature_metadata
      )
    }

    vmpa_result <- .vmpa_format_gsva_result(
      score_mat = vmpa_result,
      signature_metadata = signature_metadata,
      context = context,
      algorithm = algorithm,
      unique = unique
    )

  } else if (algorithm == "fgsea") {
    vmpa_result <- .vmpa_run_fgsea(
      stats_vec = input,
      gene_set_list = gene_set_list,
      context = context,
      gs_size = n,
      conf = min_conf,
      min_size = fgsea_min_size,
      max_size = 500L,
      n_perm_simple = 5000L,
      seed = seed
    )

    .vmpa_msg(
      verbose,
      nrow(vmpa_result), " VMPA gene sets were tested by the FGSEA backend."
    )

    gene_set_list <- .vmpa_subset_gene_sets_to_backend(
      gene_set_list,
      vmpa_result$pathway
    )
    signature_metadata <- attr(gene_set_list, "signature_metadata")

    if (!is.null(signature_metadata) &&
        is.data.frame(signature_metadata) &&
        nrow(signature_metadata) > 0L &&
        "gene_set_name" %in% colnames(signature_metadata)) {
      meta_idx <- match(vmpa_result$pathway, signature_metadata$gene_set_name)

      if ("cmap_name" %in% colnames(signature_metadata)) {
        vmpa_result$target <- signature_metadata$cmap_name[meta_idx]
      }

      if ("id" %in% colnames(signature_metadata)) {
        vmpa_result$signature_id <- signature_metadata$id[meta_idx]
      }

      if ("cps_conf_total" %in% colnames(signature_metadata)) {
        vmpa_result$cps_conf_total <- signature_metadata$cps_conf_total[meta_idx]
      }
    }

    vmpa_result <- .vmpa_format_fgsea_result(
      vmpa_result = vmpa_result,
      context = context,
      algorithm = algorithm,
      unique = unique
    )
  }

  # 5) Return result
  if (isTRUE(return_gene_sets)) {
    return(list(
      vmpa_result = vmpa_result,
      gene_sets = gene_set_list,
      signature_metadata = attr(gene_set_list, "signature_metadata"),
      context = context,
      algorithm = algorithm,
      unique = unique,
      unique_selection = attr(gene_set_list, "unique_selection")
    ))
  }

  vmpa_result
}
