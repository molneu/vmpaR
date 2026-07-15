#' Retrieve VMPA gene sets
#'
#' `vmpa_gsc()` returns VMPA gene sets for a selected cancer context
#' without running GSVA or FGSEA.
#'
#' Use this function to inspect, export, or reuse the reference gene sets that
#' are used by [vmpa()].
#'
#' @details
#' The selected cancer context determines which VMPA reference signatures are
#' used. Required reference data are resolved automatically by vmpaR.
#'
#' By default, `vmpa_gsc()` returns all signatures passing the selected
#' filters. If `unique = TRUE`, signatures are prioritized within each
#' target/protein: validated signatures are preferred when available, and among
#' the remaining candidates signatures with the highest `cps_conf_total` are
#' retained. Tied candidates are kept.
#'
#' @param context Character scalar. Cancer context to use. One of:
#'   `"glioma"`, `"melanoma"`, `"nsclc"`, `"gastric"`, `"ovarian"`,
#'   `"crc"`, `"breast"`, `"prostate"`, `"pdac"`, or `"headneck"`.
#' @param n Integer scalar. Number of bottom-ranked unique genes to include per
#'   VMPA signature. Default: `250L`.
#' @param min_conf Integer scalar. Minimum VMPA confidence level required for
#'   a signature to be included. Must be one of `1`, `2`, or `3`.
#'   Default: `1L`.
#' @param targets Optional character vector. If provided, only signatures for
#'   matching targets are retained. Default: `NULL`.
#' @param driver_filter Logical scalar. If `TRUE`, retain only signatures with
#'   a non-`"None"` cancer-driver annotation. This is a broad annotation filter,
#'   not a strict canonical-driver filter. Default: `FALSE`.
#' @param unique Logical scalar. If `TRUE`, retain the highest-priority gene
#'   sets per target/protein. If `FALSE`, keep all VMPA signatures passing
#'   the selected filters. Default: `FALSE`.
#' @param output Character scalar. Output format. One of `"list"`, `"df"`, or
#'   `"gsc"`. Default: `"list"`.
#' @param verbose Logical scalar. If `TRUE`, print progress messages. Default:
#'   `TRUE`.
#'
#' @return
#' Depending on `output`:
#'
#' - `"list"`: a named list of character vectors.
#' - `"df"`: a data frame with one column per gene set. Shorter gene sets are
#'   padded with `NA`.
#' - `"gsc"`: a `GSEABase::GeneSetCollection`.
#'
#' Signature-level metadata are attached as `attr(x, "signature_metadata")`.
#' If `unique = TRUE`, selection metadata are attached as
#' `attr(x, "unique_selection")` when available.
#'
#' @examples
#' \dontrun{
#' # Return all glioma VMPA gene sets as a named list
#' gene_sets <- vmpa_gsc("glioma")
#'
#' # Return highest-priority gene sets per target/protein
#' # Tied candidates are retained.
#' gene_sets_unique <- vmpa_gsc("glioma", unique = TRUE)
#'
#' # Inspect signature metadata
#' sig_meta <- attr(gene_sets, "signature_metadata")
#' head(sig_meta)
#'
#' # Return a GeneSetCollection
#' if (requireNamespace("GSEABase", quietly = TRUE)) {
#'   gsc <- vmpa_gsc("glioma", output = "gsc")
#' }
#' }
#'
#' @export
vmpa_gsc <- function(context,
                        n = 250L,
                        min_conf = 1L,
                        targets = NULL,
                        driver_filter = FALSE,
                        unique = FALSE,
                        output = c("list", "df", "gsc"),
                        verbose = TRUE) {
  output <- match.arg(output)
  context <- .vmpa_validate_context(context)

  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("`verbose` must be TRUE or FALSE.", call. = FALSE)
  }

  if (!is.numeric(n) || length(n) != 1L || is.na(n) ||
      n <= 0 || n != as.integer(n)) {
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

  if (!is.logical(driver_filter) || length(driver_filter) != 1L ||
      is.na(driver_filter)) {
    stop("`driver_filter` must be TRUE or FALSE.", call. = FALSE)
  }

  if (!is.logical(unique) || length(unique) != 1L || is.na(unique)) {
    stop("`unique` must be TRUE or FALSE.", call. = FALSE)
  }

  subset_file <- .vmpa_resolve_subset_file(
    context = context,
    verbose = verbose
  )

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

  if (isTRUE(unique)) {
    gene_set_list <- .vmpa_select_unique_candidates(
      gene_set_list = gene_set_list
    )
  }

  if (length(gene_set_list) == 0L) {
    stop("No VMPA gene sets available after unique reduction.", call. = FALSE)
  }

  signature_metadata <- attr(gene_set_list, "signature_metadata")
  unique_selection <- attr(gene_set_list, "unique_selection")

  if (is.data.frame(signature_metadata) &&
      nrow(signature_metadata) > 0L &&
      "gene_set_name" %in% colnames(signature_metadata)) {
    rownames(signature_metadata) <- signature_metadata$gene_set_name
  }

  if (output == "list") {
    attr(gene_set_list, "signature_metadata") <- signature_metadata
    attr(gene_set_list, "unique_selection") <- unique_selection
    return(gene_set_list)
  }

  if (output == "df") {
    max_len <- max(lengths(gene_set_list))

    out_df <- as.data.frame(
      do.call(cbind, lapply(gene_set_list, function(x) {
        length(x) <- max_len
        x
      })),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

    attr(out_df, "signature_metadata") <- signature_metadata
    attr(out_df, "unique_selection") <- unique_selection

    return(out_df)
  }

  if (!requireNamespace("GSEABase", quietly = TRUE)) {
    stop(
      "Package `GSEABase` must be installed for `output = \"gsc\"`.",
      call. = FALSE
    )
  }

  gsc <- GSEABase::GeneSetCollection(
    lapply(seq_along(gene_set_list), function(i) {
      GSEABase::GeneSet(
        geneIds = gene_set_list[[i]],
        setName = names(gene_set_list)[i]
      )
    })
  )

  attr(gsc, "signature_metadata") <- signature_metadata
  attr(gsc, "unique_selection") <- unique_selection

  gsc
}
