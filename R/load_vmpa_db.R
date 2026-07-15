# Internal subset-resolution utilities for vmpaR
# These functions are not user-facing and should not be exported.

.vmpa_resolve_subset_file <- function(context,
                                         cache_dir = NULL,
                                         verbose = TRUE) {
  context <- .vmpa_validate_context(context)
  
  if (!is.null(cache_dir) &&
      (!is.character(cache_dir) || length(cache_dir) != 1L || is.na(cache_dir) || cache_dir == "")) {
    stop("`cache_dir` must be NULL or a single non-empty character string.", call. = FALSE)
  }
  
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("`verbose` must be TRUE or FALSE.", call. = FALSE)
  }
  
  subset_file <- .vmpa_subset_filename(context)
  
  # 1) package cache ----------------------------------------------------------
  cache_dir <- .vmpa_cache_dir(cache_dir = cache_dir)
  cache_path <- file.path(cache_dir, subset_file)
  
  .vmpa_msg(verbose, "Checking cache: ", cache_path)
  
  if (file.exists(cache_path)) {
    .vmpa_msg(verbose, "Found cached subset file.")
    return(normalizePath(cache_path, winslash = "/", mustWork = TRUE))
  }
  
  # 2) download if not cached -------------------------------------------------
  subset_url <- .vmpa_subset_url(context)
  
  if (is.na(subset_url) || subset_url == "") {
    stop(
      "No download URL is configured yet for context `", context, "`.",
      call. = FALSE
    )
  }
  
  .vmpa_msg(verbose, "Subset not found in cache.")
  .vmpa_msg(verbose, "Downloading subset for context `", context, "` ...")
  
  tmp <- tempfile(pattern = paste0(context, "_subset_"), fileext = ".rds")
  on.exit(unlink(tmp), add = TRUE)
  
  download_method <- if (isTRUE(capabilities("libcurl"))) {
    "libcurl"
  } else {
    "auto"
  }
  
  download_status <- tryCatch(
    utils::download.file(
      url = subset_url,
      destfile = tmp,
      mode = "wb",
      method = download_method,
      quiet = !isTRUE(verbose)
    ),
    error = function(e) e
  )
  
  if (inherits(download_status, "error")) {
    stop(
      "Failed to download subset for context `", context, "`: ",
      conditionMessage(download_status),
      call. = FALSE
    )
  }
  
  if (!is.numeric(download_status) || length(download_status) != 1L ||
      is.na(download_status) || download_status != 0 || !file.exists(tmp)) {
    stop(
      "Failed to download subset for context `", context, "`.",
      call. = FALSE
    )
  }
  
  if (file.info(tmp)$size <= 0) {
    stop(
      "Downloaded subset file for context `", context, "` is empty.",
      call. = FALSE
    )
  }
  
  ok_copy <- file.copy(from = tmp, to = cache_path, overwrite = TRUE)
  
  if (!isTRUE(ok_copy) || !file.exists(cache_path)) {
    stop(
      "Download succeeded but writing to cache failed: ", cache_path,
      call. = FALSE
    )
  }
  
  .vmpa_msg(verbose, "Cached subset file at: ", cache_path)
  
  normalizePath(cache_path, winslash = "/", mustWork = TRUE)
}

.vmpa_cache_dir <- function(cache_dir = NULL) {
  out <- cache_dir
  
  if (is.null(out)) {
    out <- file.path(tools::R_user_dir("vmpaR", which = "cache"), "subsets")
  }
  
  if (!dir.exists(out)) {
    dir.create(out, recursive = TRUE, showWarnings = FALSE)
  }
  
  if (!dir.exists(out)) {
    stop("Could not create cache directory: ", out, call. = FALSE)
  }
  
  normalizePath(out, winslash = "/", mustWork = TRUE)
}

.vmpa_context_registry <- function() {
  data.frame(
    context = c(
      "breast",
      "crc",
      "gastric",
      "glioma",
      "headneck",
      "melanoma",
      "nsclc",
      "ovarian",
      "pdac",
      "prostate"
    ),
    subset_file = c(
      "breast_subset.rds",
      "crc_subset.rds",
      "gastric_subset.rds",
      "glioma_subset.rds",
      "headneck_subset.rds",
      "melanoma_subset.rds",
      "nsclc_subset.rds",
      "ovarian_subset.rds",
      "pdac_subset.rds",
      "prostate_subset.rds"
    ),
    subset_url = c(
      "https://ndownloader.figshare.com/files/63887454",
      "https://ndownloader.figshare.com/files/63887457",
      "https://ndownloader.figshare.com/files/63887466",
      "https://ndownloader.figshare.com/files/63887475",
      "https://ndownloader.figshare.com/files/63887463",
      "https://ndownloader.figshare.com/files/63887472",
      "https://ndownloader.figshare.com/files/63887481",
      "https://ndownloader.figshare.com/files/63887478",
      "https://ndownloader.figshare.com/files/63887469",
      "https://ndownloader.figshare.com/files/63887460"
    ),
    stringsAsFactors = FALSE
  )
}

.vmpa_subset_filename <- function(context) {
  registry <- .vmpa_context_registry()
  registry$subset_file[match(context, registry$context)]
}

.vmpa_subset_url <- function(context) {
  registry <- .vmpa_context_registry()
  registry$subset_url[match(context, registry$context)]
}

.vmpa_validate_context <- function(context) {
  valid_contexts <- .vmpa_context_registry()$context
  
  if (!is.character(context) || length(context) != 1L || is.na(context) || context == "") {
    stop("`context` must be a single non-empty character string.", call. = FALSE)
  }
  
  if (!context %in% valid_contexts) {
    stop(
      "`context` must be one of: ",
      paste(valid_contexts, collapse = ", "),
      call. = FALSE
    )
  }
  
  context
}

.vmpa_msg <- function(verbose, ...) {
  if (isTRUE(verbose)) {
    message("[vmpaR] ", ...)
  }
}