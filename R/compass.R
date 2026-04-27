#' Run COMPASS for one context and one analysis algorithm
#'
#' Main user-facing function of the protivity package.
#'
#' `compass()` resolves the subset database for the selected cancer context,
#' builds COMPASS gene sets internally, and then applies either GSVA or FGSEA
#' to the user-provided query input.
#'
#' If the required context-specific subset is already available in the package
#' cache, it is reused. Otherwise, it is downloaded automatically and cached
#' for future runs.
#'
#' @details
#' COMPASS reference subsets may contain multiple perturbation signatures for
#' the same target/protein. Internally, these signatures are identified by
#' `gene_set_name`, while their target/protein is stored in `cmap_name`.
#'
#' With `unique = FALSE`, all COMPASS signatures passing the selected filters
#' are retained as separate features.
#'
#' With `unique = TRUE`, `compass()` reduces the output to target/protein level.
#' Signatures are grouped by `cmap_name`. If validation metadata are available
#' in the signature metadata, validated signatures are prioritized. If no
#' validation metadata are available, signatures with the highest
#' `cps_conf_total` are used.
#'
#' For `algorithm = "gsva"`, all equally prioritized signatures for the same
#' target are scored and then averaged sample-wise. The returned score matrix
#' therefore has one row per target/protein, with row names such as `"AKT1"`
#' rather than signature names such as `"AKT1_..._c3"`.
#'
#' For `algorithm = "fgsea"`, one representative signature per target is
#' selected before enrichment analysis, because FGSEA p-values, adjusted
#' p-values, and leading-edge genes are signature-level results and are not
#' averaged.
#'
#' If `return_gene_sets = TRUE`, the returned list contains the gene sets used
#' for the analysis and, when `unique = TRUE`, `unique_selection`, a metadata
#' table describing which signatures were selected for target-level output.
#'
#' @param input Main query input.
#'   - For `algorithm = "gsva"`: numeric matrix or data.frame with genes in rows
#'     and samples in columns.
#'   - For `algorithm = "fgsea"`: named numeric vector of gene-level statistics
#'     or rankings.
#' @param context Character; preferred context, one of:
#'   `"glioma"`, `"melanoma"`, `"nsclc"`, `"gastric"`, `"ovarian"`,
#'   `"crc"`, `"breast"`, `"prostate"`, `"pdac"`, `"headneck"`.
#' @param algorithm Character scalar. Either `"gsva"` or `"fgsea"`.
#'   Default: `"gsva"`.
#' @param gsva_score_scaling Character scalar. Post-processing applied to GSVA
#'   score matrices. One of `"none"`, `"sample_z"`, `"sample_pop_sd"`, or
#'   `"signature_z"`. `"sample_z"` applies z-score scaling within each sample.
#'   `"sample_pop_sd"` applies sample-wise scaling using the population standard
#'   deviation. `"signature_z"` applies z-score scaling across samples for each
#'   signature. Default: `"none"`. Only used when `algorithm = "gsva"`.
#' @param unique Logical. If `TRUE`, reduce COMPASS output to one
#'   target/protein-level result per target. Multiple COMPASS signatures for
#'   the same target are grouped by `cmap_name`. Validation metadata are used
#'   for prioritization if available; otherwise signatures with the highest
#'   `cps_conf_total` are used. For `algorithm = "gsva"`, equally prioritized
#'   signatures are averaged at the score level. For `algorithm = "fgsea"`,
#'   one representative signature per target is selected before enrichment
#'   analysis. If `FALSE`, all COMPASS signatures are retained separately.
#'   Default: `TRUE`.
#' @param n Integer. Number of bottom-ranked genes to include per reference
#'   signature when building COMPASS gene sets. Default: `250L`.
#' @param min_conf Integer. Minimum `cps_conf_total` required for a
#'   reference signature to be included. Must be one of `1`, `2`, or `3`.
#'   Default: `1L`.
#' @param targets Optional character vector. If provided, only signatures with
#'   matching targets are retained. Default: `NULL`.
#' @param driver_filter Logical. If `TRUE`, only signatures with
#'   `cancer_driver_summary != "None"` are retained. This is a broad filter:
#'   any non-`"None"` driver-related annotation is kept, not only canonical
#'   drivers. Default: `FALSE`.
#' @param return_gene_sets Logical. If `TRUE`, return a list containing the
#'   COMPASS result, the gene sets used for the analysis, and signature
#'   selection metadata. If `unique = FALSE`, `gene_sets` contains all gene sets
#'   passing the filters. If `unique = TRUE`, `gene_sets` contains the
#'   prioritized gene sets used to compute the target-level output. For GSVA,
#'   several gene sets may be retained for one target if they share the same
#'   highest priority; their scores are averaged in `compass_result`.
#'   Default: `FALSE`.
#' @param verbose Logical. If `TRUE`, print progress messages for subset
#'   resolution, download steps, and optional GSVA score scaling. Default: `TRUE`.
#'
#' @return
#' If `return_gene_sets = FALSE`:
#' - for `algorithm = "gsva"`: a numeric matrix of COMPASS scores. With
#'   `unique = TRUE`, rows correspond to target/protein names. With
#'   `unique = FALSE`, rows correspond to individual COMPASS signature names.
#' - for `algorithm = "fgsea"`: a data.frame of FGSEA results, annotated with
#'   target and signature metadata when available.
#'
#' If `return_gene_sets = TRUE`, a list with:
#' - `compass_result`: COMPASS score matrix or FGSEA result table
#' - `gene_sets`: gene sets used for the analysis
#' - `context`: selected cancer context
#' - `algorithm`: selected analysis algorithm
#' - `unique`: whether target-level reduction was used
#' - `unique_selection`: metadata describing the selected signatures when
#'   `unique = TRUE`; otherwise `NULL`
#'
#' @examples
#' \dontrun{
#' if (requireNamespace("Biobase", quietly = TRUE)) {
#'   data(kebir_gb, package = "protivity")
#'
#'   # Example 1: GSVA workflow
#'   df <- Biobase::exprs(kebir_gb)
#'
#'   gsva_result <- compass(
#'     input = df,
#'     context = "glioma",
#'     algorithm = "gsva"
#'   )
#'
#'   # Sample metadata are available via pData()
#'   sample_metadata <- Biobase::pData(kebir_gb)
#'
#'   # Example 2: FGSEA workflow
#'   # Build a simple ranked vector contrasting relapse vs treatment-naive
#'   # samples. This is a minimal example for demonstrating the required
#'   # input format for algorithm = "fgsea".
#'   relapsed <- sample_metadata$relapse_TYPE != "n"
#'   naive <- sample_metadata$relapse_TYPE == "n"
#'
#'   stats_vec <- rowMeans(df[, relapsed, drop = FALSE]) -
#'     rowMeans(df[, naive, drop = FALSE])
#'
#'   stats_vec <- stats_vec[is.finite(stats_vec) & !is.na(stats_vec)]
#'   stats_vec <- sort(stats_vec, decreasing = TRUE)
#'
#'   fgsea_result <- compass(
#'     input = stats_vec,
#'     context = "glioma",
#'     algorithm = "fgsea"
#'   )
#' }
#' }
#'
#' @export
compass <- function(input,
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
                    verbose = TRUE) {
  algorithm <- match.arg(algorithm)
  gsva_score_scaling <- match.arg(gsva_score_scaling)
  context <- .compass_validate_context(context)

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
  subset_file <- .compass_resolve_subset_file(
    context = context,
    verbose = verbose
  )

  # 2) Read subset GCT and build COMPASS gene sets internally
  gct <- .compass_read_subset_gct(subset_file)

  gene_set_list <- .compass_build_gene_sets(
    gct = gct,
    n = n,
    min_conf = min_conf,
    targets = targets,
    driver_filter = driver_filter
  )

  if (length(gene_set_list) == 0L) {
    stop("No COMPASS gene sets available after filtering.", call. = FALSE)
  }

  # 3) Optionally reduce multiple signatures per target/protein
  if (isTRUE(unique)) {
    if (algorithm == "gsva") {
      gene_set_list <- .compass_select_unique_candidates(
        gene_set_list = gene_set_list,
        mode = "aggregate"
      )
    } else if (algorithm == "fgsea") {
      gene_set_list <- .compass_select_unique_candidates(
        gene_set_list = gene_set_list,
        mode = "representative"
      )
    }
  }

  if (length(gene_set_list) == 0L) {
    stop("No COMPASS gene sets available after unique reduction.", call. = FALSE)
  }

  # 4) Run analysis with the selected algorithm
  if (algorithm == "gsva") {
    compass_result <- .compass_run_gsva(
      expr_mat = input,
      gene_set_list = gene_set_list,
      score_scaling = gsva_score_scaling,
      verbose = verbose
    )

    if (isTRUE(unique)) {
      compass_result <- .compass_reduce_unique_gsva_scores(
        score_mat = compass_result,
        signature_metadata = attr(gene_set_list, "signature_metadata")
      )
    }

  } else if (algorithm == "fgsea") {
    compass_result <- .compass_run_fgsea(
      stats_vec = input,
      gene_set_list = gene_set_list,
      context = context,
      gs_size = n,
      conf = min_conf,
      min_size = 10L,
      max_size = 500L,
      n_perm_simple = 5000L
    )

    signature_metadata <- attr(gene_set_list, "signature_metadata")

    if (!is.null(signature_metadata) &&
        is.data.frame(signature_metadata) &&
        nrow(signature_metadata) > 0L &&
        "gene_set_name" %in% colnames(signature_metadata)) {
      meta_idx <- match(compass_result$pathway, signature_metadata$gene_set_name)

      if ("cmap_name" %in% colnames(signature_metadata)) {
        compass_result$target <- signature_metadata$cmap_name[meta_idx]
      }

      if ("id" %in% colnames(signature_metadata)) {
        compass_result$signature_id <- signature_metadata$id[meta_idx]
      }

      if ("cps_conf_total" %in% colnames(signature_metadata)) {
        compass_result$cps_conf_total <- signature_metadata$cps_conf_total[meta_idx]
      }
    }
  }

  # 5) Return result
  if (isTRUE(return_gene_sets)) {
    return(list(
      compass_result = compass_result,
      gene_sets = gene_set_list,
      context = context,
      algorithm = algorithm,
      unique = unique,
      unique_selection = attr(gene_set_list, "unique_selection")
    ))
  }

  compass_result
}
