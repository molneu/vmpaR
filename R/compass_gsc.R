#' Build COMPASS gene sets for one context
#'
#' Public helper to extract COMPASS reference gene sets without running
#' downstream scoring. This function reuses the same internal subset-resolution
#' and gene-set-building logic as `compass()`.
#'
#' @param context Character scalar. Cancer context used to resolve the
#'   corresponding subset database.
#' @param subset_dir Optional character scalar. Local directory containing
#'   subset files such as `glioma_subset.rds`. Checked before cache/download.
#' @param cache_dir Optional character scalar. Override for the package cache
#'   directory. If `NULL`, a package-specific user cache directory is used.
#' @param download_if_missing Logical. If `TRUE`, missing subset files are
#'   downloaded and cached automatically. Default: `TRUE`.
#' @param n Integer. Number of bottom-ranked unique genes to include per
#'   reference signature. Default: `250L`.
#' @param min_conf Integer. Minimum `cps_conf_total` required for a
#'   reference signature to be included. Must be one of `1`, `2`, or `3`.
#'   Default: `1L`.
#' @param targets Optional character vector. If provided, only signatures with
#'   matching targets are retained. Default: `NULL`.
#' @param driver_filter Logical. If `TRUE`, only signatures with
#'   `cancer_driver_summary != "None"` are retained. This is a broad filter:
#'   any non-`"None"` driver-related annotation is kept, not only canonical
#'   drivers. Default: `FALSE`.
#' @param output Character scalar. Output format. One of `"list"`, `"df"`,
#'   or `"gsc"`. Default: `"list"`.
#' @param verbose Logical. If `TRUE`, print progress messages for subset
#'   resolution and download steps. Default: `TRUE`.
#'
#' @return
#' Depending on `output`:
#' - `"list"`: named list of character vectors
#' - `"df"`: data frame with one column per gene set
#' - `"gsc"`: `GSEABase::GeneSetCollection`
#'
#' For all outputs, signature-level metadata are attached as
#' `attr(x, "signature_metadata")`. These metadata describe the reference
#' signatures underlying the returned gene sets and may include fields such as
#' `gene_set_name`, `id`, `cmap_name`, `cps_conf_total`, and
#' `cancer_driver_summary`. Row names correspond to gene-set names.
#'
#' @examples
#' gsc <- compass_gsc("glioma", output = "gsc")
#' sig_meta <- attr(gsc, "signature_metadata")
#' head(sig_meta)
#'
#' @export
compass_gsc <- function(context,
                        subset_dir = NULL,
                        cache_dir = NULL,
                        download_if_missing = TRUE,
                        n = 250L,
                        min_conf = 1L,
                        targets = NULL,
                        driver_filter = FALSE,
                        output = c("list", "df", "gsc"),
                        verbose = TRUE) {
  output <- match.arg(output)
  context <- .compass_validate_context(context)

  if (!is.logical(download_if_missing) || length(download_if_missing) != 1L || is.na(download_if_missing)) {
    stop("`download_if_missing` must be TRUE or FALSE.", call. = FALSE)
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

  subset_file <- .compass_resolve_subset_file(
    context = context,
    subset_dir = subset_dir,
    cache_dir = cache_dir,
    download_if_missing = download_if_missing,
    verbose = verbose
  )

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

  signature_metadata <- attr(gene_set_list, "signature_metadata")

  if (!is.null(signature_metadata) && nrow(signature_metadata) > 0L) {
    rownames(signature_metadata) <- signature_metadata$gene_set_name
  }

  if (output == "list") {
    attr(gene_set_list, "signature_metadata") <- signature_metadata
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
    return(out_df)
  }

  if (!requireNamespace("GSEABase", quietly = TRUE)) {
    stop("Package `GSEABase` must be installed for `output = \"gsc\"`.", call. = FALSE)
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

  gsc
}
