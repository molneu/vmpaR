# Internal subset-resolution utilities for protivity
# These functions are not user-facing and should not be exported.

.compass_resolve_subset_file <- function(context,
                                         subset_dir = NULL,
                                         cache_dir = NULL,
                                         download_if_missing = TRUE,
                                         verbose = TRUE) {
  context <- .compass_validate_context(context)
  
  if (!is.null(subset_dir) &&
      (!is.character(subset_dir) || length(subset_dir) != 1L || is.na(subset_dir) || subset_dir == "")) {
    stop("`subset_dir` must be NULL or a single non-empty character string.", call. = FALSE)
  }
  
  if (!is.null(cache_dir) &&
      (!is.character(cache_dir) || length(cache_dir) != 1L || is.na(cache_dir) || cache_dir == "")) {
    stop("`cache_dir` must be NULL or a single non-empty character string.", call. = FALSE)
  }
  
  if (!is.logical(download_if_missing) || length(download_if_missing) != 1L || is.na(download_if_missing)) {
    stop("`download_if_missing` must be TRUE or FALSE.", call. = FALSE)
  }
  
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("`verbose` must be TRUE or FALSE.", call. = FALSE)
  }
  
  subset_file <- .compass_subset_filename(context)
  
  # 1) explicit local directory ----------------------------------------------
  if (!is.null(subset_dir)) {
    candidate <- file.path(subset_dir, subset_file)
    .compass_msg(verbose, "Checking subset_dir: ", candidate)
    
    if (file.exists(candidate)) {
      .compass_msg(verbose, "Found local subset file.")
      return(normalizePath(candidate, winslash = "/", mustWork = TRUE))
    }
  }
  
  # 2) package cache ----------------------------------------------------------
  cache_dir <- .compass_cache_dir(cache_dir = cache_dir)
  cache_path <- file.path(cache_dir, subset_file)
  
  .compass_msg(verbose, "Checking cache: ", cache_path)
  
  if (file.exists(cache_path)) {
    .compass_msg(verbose, "Found cached subset file.")
    return(normalizePath(cache_path, winslash = "/", mustWork = TRUE))
  }
  
  # 3) optional download ------------------------------------------------------
  if (!isTRUE(download_if_missing)) {
    stop(
      "Subset file for context `", context, "` was not found.\n",
      "Checked:\n",
      if (!is.null(subset_dir)) paste0("- subset_dir: ", file.path(subset_dir, subset_file), "\n") else "",
      "- cache: ", cache_path, "\n",
      "Set `download_if_missing = TRUE` or provide a valid `subset_dir`.",
      call. = FALSE
    )
  }
  
  subset_url <- .compass_subset_url(context)
  
  if (is.na(subset_url) || subset_url == "") {
    stop(
      "No download URL is configured yet for context `", context, "`.",
      call. = FALSE
    )
  }
  
  .compass_msg(verbose, "Downloading subset for context `", context, "` ...")
  
  tmp <- tempfile(pattern = paste0(context, "_subset_"), fileext = ".rds")
  on.exit(unlink(tmp), add = TRUE)
  
  download_status <- tryCatch(
    utils::download.file(
      url = subset_url,
      destfile = tmp,
      mode = "wb",
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
  
  ok_copy <- file.copy(from = tmp, to = cache_path, overwrite = TRUE)
  
  if (!isTRUE(ok_copy) || !file.exists(cache_path)) {
    stop(
      "Download succeeded but writing to cache failed: ", cache_path,
      call. = FALSE
    )
  }
  
  .compass_msg(verbose, "Cached subset file at: ", cache_path)
  
  normalizePath(cache_path, winslash = "/", mustWork = TRUE)
}

.compass_cache_dir <- function(cache_dir = NULL) {
  out <- cache_dir
  
  if (is.null(out)) {
    out <- file.path(tools::R_user_dir("protivity", which = "cache"), "subsets")
  }
  
  if (!dir.exists(out)) {
    dir.create(out, recursive = TRUE, showWarnings = FALSE)
  }
  
  if (!dir.exists(out)) {
    stop("Could not create cache directory: ", out, call. = FALSE)
  }
  
  normalizePath(out, winslash = "/", mustWork = TRUE)
}

.compass_context_registry <- function() {
  contexts <- c(
    "breast", "crc", "gastric", "glioma", "headneck",
    "melanoma", "nsclc", "ovarian", "pdac", "prostate"
  )
  
  data.frame(
    context = contexts,
    subset_file = paste0(contexts, "_subset.rds"),
    subset_url = NA_character_,
    stringsAsFactors = FALSE
  )
}

.compass_subset_filename <- function(context) {
  registry <- .compass_context_registry()
  registry$subset_file[match(context, registry$context)]
}

.compass_subset_url <- function(context) {
  registry <- .compass_context_registry()
  registry$subset_url[match(context, registry$context)]
}

.compass_validate_context <- function(context) {
  valid_contexts <- .compass_context_registry()$context
  
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

.compass_msg <- function(verbose, ...) {
  if (isTRUE(verbose)) {
    message("[protivity] ", ...)
  }
}