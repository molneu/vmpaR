# Create small FGSEA-ready example dataset for protivity
#
# Source:
#   COMPASS paper reproducibility / Figure 4e context
#   Source object: drugs_U251MG_withLINCsignatures.cleaned.rds
#
# Selected example:
#   Drug: AZD-8055
#   Dose: 0.12 uM
#   Target context: MTOR
#
# Output object written to data/:
#   - example_u251_azd8055_rank
#
# Notes:
#   - This script is intended to be run from the package root
#     (the directory containing DESCRIPTION).
#   - By default, it uses a local Windows path.
#   - Optionally, you can set an environment variable
#     PROTIVITY_U251_DRUGS_GCT to point to the source RDS file.

# ---- input path ----------------------------------------------------------

input_path <- Sys.getenv(
  "PROTIVITY_U251_DRUGS_GCT",
  unset = "C:/Users/ACER/Desktop/Uni Essen/SHK Neurochirurgie/Compass reproducibility/compass_reproducibility/reproducibility/figure4/data/drugs_U251MG_withLINCsignatures.cleaned.rds"
)

# ---- project checks ------------------------------------------------------

if (!file.exists("DESCRIPTION")) {
  stop(
    "Please run this script from the package root (directory containing DESCRIPTION).",
    call. = FALSE
  )
}

# ---- package checks ------------------------------------------------------

required_pkgs <- c("cmapR", "usethis")

missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs) > 0) {
  stop(
    "Missing required package(s): ",
    paste(missing_pkgs, collapse = ", "),
    call. = FALSE
  )
}

# ---- helpers -------------------------------------------------------------

.normalize_dose <- function(x) {
  x <- as.character(x)
  x <- gsub("\u00B5", "u", x, fixed = TRUE)
  x <- gsub("\\s+", " ", x)
  trimws(x)
}

# ---- read source object --------------------------------------------------

if (!file.exists(input_path)) {
  stop("Input file does not exist: ", input_path, call. = FALSE)
}

gct <- readRDS(input_path)

if (!methods::is(gct, "GCT")) {
  stop("Input object is not a `GCT` object.", call. = FALSE)
}

# ---- basic validation ----------------------------------------------------

required_rdesc_cols <- c("symbol")
required_cdesc_cols <- c("sig_id", "cmap_name", "pert_idose", "tas")

missing_rdesc <- setdiff(required_rdesc_cols, colnames(gct@rdesc))
missing_cdesc <- setdiff(required_cdesc_cols, colnames(gct@cdesc))

if (length(missing_rdesc) > 0) {
  stop(
    "Missing required row metadata column(s): ",
    paste(missing_rdesc, collapse = ", "),
    call. = FALSE
  )
}

if (length(missing_cdesc) > 0) {
  stop(
    "Missing required column metadata column(s): ",
    paste(missing_cdesc, collapse = ", "),
    call. = FALSE
  )
}

if (!is.matrix(gct@mat) || !is.numeric(gct@mat)) {
  stop("`gct@mat` is not a numeric matrix.", call. = FALSE)
}

if (nrow(gct@mat) != nrow(gct@rdesc)) {
  stop("Row metadata and matrix row count do not match.", call. = FALSE)
}

if (ncol(gct@mat) != nrow(gct@cdesc)) {
  stop("Column metadata and matrix column count do not match.", call. = FALSE)
}

# ---- select AZD-8055 example signature -----------------------------------

cdesc <- as.data.frame(gct@cdesc, stringsAsFactors = FALSE)
cdesc$pert_idose_norm <- .normalize_dose(cdesc$pert_idose)

keep_idx <- which(
  cdesc$cmap_name == "AZD-8055" &
    cdesc$pert_idose_norm == "0.12 uM"
)

if (length(keep_idx) == 0L) {
  stop(
    "No signature found for cmap_name = 'AZD-8055' and pert_idose = '0.12 uM'.",
    call. = FALSE
  )
}

if (length(keep_idx) > 1L) {
  matched_ids <- cdesc$sig_id[keep_idx]
  stop(
    "Expected exactly one signature for AZD-8055 at 0.12 uM, found ",
    length(keep_idx), ". Matching sig_id values: ",
    paste(matched_ids, collapse = ", "),
    call. = FALSE
  )
}

sel_idx <- keep_idx[[1]]

# ---- build ranking vector ------------------------------------------------

stats_vec <- as.numeric(gct@mat[, sel_idx])
gene_symbols <- as.character(gct@rdesc$symbol)

if (length(stats_vec) != length(gene_symbols)) {
  stop("Length mismatch between selected signature values and gene symbols.", call. = FALSE)
}

if (anyNA(gene_symbols) || any(gene_symbols == "")) {
  stop("Missing or empty gene symbols found in `gct@rdesc$symbol`.", call. = FALSE)
}

dup_symbols <- unique(gene_symbols[duplicated(gene_symbols)])

if (length(dup_symbols) > 0L) {
  warning(
    "Duplicated gene symbols were excluded from example_u251_azd8055_rank. ",
    "Examples: ",
    paste(utils::head(dup_symbols, 10), collapse = ", "),
    if (length(dup_symbols) > 10) " ..." else "",
    call. = FALSE
  )

  keep_unique <- !(gene_symbols %in% dup_symbols)
  stats_vec <- stats_vec[keep_unique]
  gene_symbols <- gene_symbols[keep_unique]
}

names(stats_vec) <- gene_symbols

example_u251_azd8055_rank <- stats_vec[is.finite(stats_vec) & !is.na(stats_vec)]
example_u251_azd8055_rank <- sort(example_u251_azd8055_rank, decreasing = TRUE)

# ---- final checks --------------------------------------------------------

if (!is.numeric(example_u251_azd8055_rank)) {
  stop("`example_u251_azd8055_rank` is not numeric.", call. = FALSE)
}

if (is.null(names(example_u251_azd8055_rank))) {
  stop("`example_u251_azd8055_rank` has no gene names.", call. = FALSE)
}

if (anyNA(names(example_u251_azd8055_rank)) || any(names(example_u251_azd8055_rank) == "")) {
  stop("`example_u251_azd8055_rank` contains missing or empty gene names.", call. = FALSE)
}

if (anyDuplicated(names(example_u251_azd8055_rank)) > 0L) {
  stop("`example_u251_azd8055_rank` contains duplicated gene names.", call. = FALSE)
}

if (length(example_u251_azd8055_rank) == 0L) {
  stop("`example_u251_azd8055_rank` is empty after filtering.", call. = FALSE)
}

if (anyNA(example_u251_azd8055_rank) || any(!is.finite(example_u251_azd8055_rank))) {
  stop("`example_u251_azd8055_rank` contains NA, NaN, or infinite values.", call. = FALSE)
}

# ---- informative messages ------------------------------------------------

message(
  "Selected signature: ",
  cdesc$sig_id[sel_idx],
  " | drug = ", cdesc$cmap_name[sel_idx],
  " | dose = ", cdesc$pert_idose_norm[sel_idx],
  " | tas = ", cdesc$tas[sel_idx]
)

message(
  "Created example_u251_azd8055_rank with ",
  length(example_u251_azd8055_rank),
  " ranked genes."
)

# ---- save into package data/ ---------------------------------------------

usethis::use_data(
  example_u251_azd8055_rank,
  overwrite = TRUE,
  compress = "xz"
)

message("Saved package dataset: example_u251_azd8055_rank")
