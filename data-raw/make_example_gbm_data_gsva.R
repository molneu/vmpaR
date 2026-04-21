# Create small GSVA-ready example datasets for protivity
#
# Source:
#   GEO accession GSE145128
#   Source object: ExpressionSet_GSE145128.rds
#   Context: COMPASS paper / Figure-6-like relapse comparison
#
# Filtering applied:
#   1) exclude relapse_TYPE == "eR"
#   2) exclude patients B046 and B078
#
# Output objects written to data/:
#   - example_gbm_expr
#   - example_gbm_sample_info
#
# Notes:
#   - This script is intended to be run from the package root
#     (the directory containing DESCRIPTION).
#   - By default, it uses a local Windows path.
#   - Optionally, you can set an environment variable
#     PROTIVITY_GSE145128_ESET to point to the source RDS file.

# ---- input path ----------------------------------------------------------

input_path <- Sys.getenv(
  "PROTIVITY_GSE145128_ESET",
  unset = "C:/Users/ACER/Desktop/Uni Essen/SHK Neurochirurgie/Compass reproducibility/compass_reproducibility/reproducibility/figure6/data/ExpressionSet_GSE145128.rds"
)

# ---- project checks ------------------------------------------------------

if (!file.exists("DESCRIPTION")) {
  stop(
    "Please run this script from the package root (directory containing DESCRIPTION).",
    call. = FALSE
  )
}

# ---- package checks ------------------------------------------------------

required_pkgs <- c("Biobase", "usethis")

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

# ---- read source object --------------------------------------------------

if (!file.exists(input_path)) {
  stop("Input file does not exist: ", input_path, call. = FALSE)
}

eset <- readRDS(input_path)

if (!methods::is(eset, "ExpressionSet")) {
  stop("Input object is not an ExpressionSet.", call. = FALSE)
}

expr <- Biobase::exprs(eset)
pheno <- Biobase::pData(eset)

# ---- basic validation ----------------------------------------------------

required_cols <- c(
  "patient_ID",
  "relapse_TYPE",
  "MGMT_STATUS",
  "geo_accession",
  "sample type:ch1"
)

missing_cols <- setdiff(required_cols, colnames(pheno))

if (length(missing_cols) > 0) {
  stop(
    "Missing required phenotype column(s): ",
    paste(missing_cols, collapse = ", "),
    call. = FALSE
  )
}

if (!is.matrix(expr) || !is.numeric(expr)) {
  stop("Expression data extracted from ExpressionSet is not a numeric matrix.", call. = FALSE)
}

if (is.null(rownames(expr))) {
  stop("Expression matrix has no gene names (rownames).", call. = FALSE)
}

if (is.null(colnames(expr))) {
  stop("Expression matrix has no sample names (colnames).", call. = FALSE)
}

if (anyDuplicated(rownames(expr)) > 0) {
  stop("Duplicate gene names found in source expression matrix.", call. = FALSE)
}

if (anyDuplicated(colnames(expr)) > 0) {
  stop("Duplicate sample names found in source expression matrix.", call. = FALSE)
}

if (!all(colnames(expr) %in% rownames(pheno))) {
  stop(
    "Not all expression-matrix sample names are present in pData rownames.",
    call. = FALSE
  )
}

# Reorder phenotype data to expression-matrix column order
pheno <- pheno[colnames(expr), , drop = FALSE]

# Explicitly validate filter-critical columns
filter_cols <- c("patient_ID", "relapse_TYPE")

if (anyNA(pheno[, filter_cols, drop = FALSE])) {
  stop(
    "Missing values found in filter-critical phenotype columns: ",
    paste(filter_cols, collapse = ", "),
    call. = FALSE
  )
}

# ---- Figure-6-like filtering ---------------------------------------------

excluded_patients <- c("B046", "B078")

keep_idx <- pheno[["relapse_TYPE"]] != "eR" &
  !pheno[["patient_ID"]] %in% excluded_patients

if (anyNA(keep_idx)) {
  stop("Filtering produced NA values in keep_idx.", call. = FALSE)
}

example_gbm_expr <- expr[, keep_idx, drop = FALSE]
example_gbm_sample_info <- pheno[keep_idx, required_cols, drop = FALSE]

# ---- tidy sample metadata ------------------------------------------------

example_gbm_sample_info <- data.frame(
  sample_id = colnames(example_gbm_expr),
  patient_ID = example_gbm_sample_info[["patient_ID"]],
  relapse_TYPE = example_gbm_sample_info[["relapse_TYPE"]],
  MGMT_STATUS = example_gbm_sample_info[["MGMT_STATUS"]],
  geo_accession = example_gbm_sample_info[["geo_accession"]],
  sample_type = example_gbm_sample_info[["sample type:ch1"]],
  stringsAsFactors = FALSE,
  row.names = NULL
)

# Validate expected relapse labels before ordering
valid_relapse_types <- c("n", "cR")

unexpected_relapse_types <- setdiff(
  unique(example_gbm_sample_info$relapse_TYPE),
  valid_relapse_types
)

if (length(unexpected_relapse_types) > 0) {
  stop(
    "Unexpected relapse_TYPE value(s) after filtering: ",
    paste(unexpected_relapse_types, collapse = ", "),
    call. = FALSE
  )
}

# Stable ordering: patient_ID, then n before cR
example_gbm_sample_info$relapse_TYPE <- factor(
  example_gbm_sample_info$relapse_TYPE,
  levels = valid_relapse_types
)

ord <- order(
  example_gbm_sample_info$patient_ID,
  example_gbm_sample_info$relapse_TYPE
)

example_gbm_sample_info <- example_gbm_sample_info[ord, , drop = FALSE]
example_gbm_expr <- example_gbm_expr[, example_gbm_sample_info$sample_id, drop = FALSE]

# ---- semantic structure checks -------------------------------------------

# Expected filtered structure from your Figure-6-like summary:
#   - 10 samples total
#   - 5 patients total
#   - each patient represented exactly once in n and once in cR

if (ncol(example_gbm_expr) != 10) {
  stop(
    "Expected 10 samples after filtering, found ",
    ncol(example_gbm_expr), ".",
    call. = FALSE
  )
}

patient_counts <- table(example_gbm_sample_info$patient_ID)

if (length(patient_counts) != 5) {
  stop(
    "Expected 5 patients after filtering, found ",
    length(patient_counts), ".",
    call. = FALSE
  )
}

if (!all(patient_counts == 2)) {
  stop(
    "Expected exactly 2 samples per patient after filtering.",
    call. = FALSE
  )
}

pair_tab <- table(
  example_gbm_sample_info$patient_ID,
  example_gbm_sample_info$relapse_TYPE
)

if (!all(pair_tab[, valid_relapse_types, drop = FALSE] == 1)) {
  stop(
    "Expected exactly one 'n' and one 'cR' sample per patient after filtering.",
    call. = FALSE
  )
}

# Convert factor back to character for cleaner package dataset behavior
example_gbm_sample_info$relapse_TYPE <- as.character(example_gbm_sample_info$relapse_TYPE)

# ---- final checks --------------------------------------------------------

if (!is.matrix(example_gbm_expr) || !is.numeric(example_gbm_expr)) {
  stop("example_gbm_expr is not a numeric matrix.", call. = FALSE)
}

if (ncol(example_gbm_expr) != nrow(example_gbm_sample_info)) {
  stop(
    "Number of samples in expression matrix and sample info do not match.",
    call. = FALSE
  )
}

if (!identical(colnames(example_gbm_expr), example_gbm_sample_info$sample_id)) {
  stop(
    "Sample order mismatch between expression matrix and sample info.",
    call. = FALSE
  )
}

if (anyDuplicated(rownames(example_gbm_expr)) > 0) {
  stop("Duplicate gene names found in example_gbm_expr.", call. = FALSE)
}

if (anyDuplicated(colnames(example_gbm_expr)) > 0) {
  stop("Duplicate sample names found in example_gbm_expr.", call. = FALSE)
}

if (anyNA(example_gbm_expr)) {
  stop("Missing values found in example_gbm_expr.", call. = FALSE)
}

if (anyNA(example_gbm_sample_info)) {
  stop("Missing values found in example_gbm_sample_info.", call. = FALSE)
}

# ---- informative messages ------------------------------------------------

message(
  "Created example_gbm_expr with dimensions: ",
  nrow(example_gbm_expr), " genes x ", ncol(example_gbm_expr), " samples"
)

message(
  "Patients retained: ",
  paste(unique(example_gbm_sample_info$patient_ID), collapse = ", ")
)

message(
  "Relapse groups: ",
  paste(
    names(table(example_gbm_sample_info$relapse_TYPE)),
    table(example_gbm_sample_info$relapse_TYPE),
    sep = "=",
    collapse = ", "
  )
)

# ---- save into package data/ ---------------------------------------------

usethis::use_data(
  example_gbm_expr,
  example_gbm_sample_info,
  overwrite = TRUE,
  compress = "xz"
)

message("Saved package datasets: example_gbm_expr, example_gbm_sample_info")