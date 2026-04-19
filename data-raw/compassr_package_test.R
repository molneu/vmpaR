# ============================================================================
# test_compass_workflow.R
# ----------------------------------------------------------------------------
# Purpose:
#   End-to-end walkthrough for the current COMPASS workflow:
#   1) compass_gsc()
#   2) compass()
#   3) compass_summary()
#   4) optional compass_annotate()
#
# This script is intentionally verbose and modifiable.
# It is meant as:
#   - a runnable test script
#   - a usage example
#   - a mini internal guide for what each function expects
#
# IMPORTANT:
#   - GSVA mode works directly with an expression matrix
#   - FGSEA mode requires a named numeric ranking vector
#   - AUCell is not implemented yet
# ============================================================================


# ============================================================================
# 0) USER INPUT: set your project root here
# ============================================================================

project_root <- "C:/Users/ACER/Desktop/Uni Essen/SHK Neurochirurgie/Compensatory Mechanisms/COMPASS R Package"

project_root <- normalizePath(project_root, winslash = "/", mustWork = TRUE)


# ============================================================================
# 1) Paths and options
# ============================================================================

subset_dir <- file.path(project_root, "subsets")
if (!file.exists(subset_dir)) {
  subset_dir <- file.path(project_root, "Data", "subsets")
}

if (!file.exists(subset_dir)) {
  stop("Could not find subset directory.", call. = FALSE)
}

go_annotation_path <- file.path(project_root, "annotation", "compass_go_annotations.rds")
expr_xlsx_path <- file.path(project_root, "Data", "input", "expr_logcpm.xlsx")

results_dir <- file.path(project_root, "results")
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

# Set package-style options so functions can find local resources automatically
options(compass.subset_dir = subset_dir)

if (file.exists(go_annotation_path)) {
  options(compass.go_annotation_path = go_annotation_path)
}


# ============================================================================
# 2) Source current development functions
#    (adjust this if/when the package becomes formally installed)
# ============================================================================

source(file.path(project_root, "R", "compass_gsc.R"))
source(file.path(project_root, "R", "compass.R"))
source(file.path(project_root, "R", "compass_summary.R"))

if (file.exists(file.path(project_root, "R", "compass_annotate.R"))) {
  source(file.path(project_root, "R", "compass_annotate.R"))
}


# ============================================================================
# 3) Helper: load expression matrix from xlsx
# ----------------------------------------------------------------------------
# Assumption:
#   - first column = gene symbols
#   - remaining columns = numeric expression values
#
# If your xlsx structure changes, this is the only part you need to adapt.
# ============================================================================

load_expr_matrix_from_xlsx <- function(xlsx_path) {
  if (!requireNamespace("readxl", quietly = TRUE)) {
    stop("Package `readxl` must be installed to read the xlsx input.", call. = FALSE)
  }
  
  if (!file.exists(xlsx_path)) {
    stop("Expression xlsx file not found: ", xlsx_path, call. = FALSE)
  }
  
  df <- readxl::read_xlsx(xlsx_path)
  
  if (ncol(df) < 2L) {
    stop("Input xlsx must contain at least one gene column and one sample column.", call. = FALSE)
  }
  
  gene_col <- names(df)[1]
  
  genes <- as.character(df[[gene_col]])
  expr_df <- df[, -1, drop = FALSE]
  
  expr_mat <- as.matrix(expr_df)
  storage.mode(expr_mat) <- "numeric"
  rownames(expr_mat) <- genes
  
  # Basic checks
  if (anyNA(rownames(expr_mat)) || any(rownames(expr_mat) == "")) {
    stop("Gene column contains missing or empty gene names.", call. = FALSE)
  }
  
  if (anyDuplicated(rownames(expr_mat)) > 0L) {
    dup <- unique(rownames(expr_mat)[duplicated(rownames(expr_mat))])
    stop(
      "Duplicated gene names found in expression matrix. Examples: ",
      paste(utils::head(dup, 10), collapse = ", "),
      call. = FALSE
    )
  }
  
  if (is.null(colnames(expr_mat)) || any(colnames(expr_mat) == "")) {
    stop("Expression matrix must have valid sample names as column names.", call. = FALSE)
  }
  
  if (any(!is.finite(expr_mat), na.rm = TRUE)) {
    stop("Expression matrix contains non-finite values.", call. = FALSE)
  }
  
  expr_mat
}


# ============================================================================
# 4) Load example expression input
# ----------------------------------------------------------------------------
# INPUT TYPE for GSVA mode:
#   numeric matrix
#   rows = genes
#   columns = samples
# ============================================================================

expr_mat <- load_expr_matrix_from_xlsx(expr_xlsx_path)

cat("\n=== Expression matrix loaded ===\n")
cat("Dimensions:", nrow(expr_mat), "genes x", ncol(expr_mat), "samples\n")
cat("First genes:\n")
print(utils::head(rownames(expr_mat), 10))
cat("First samples:\n")
print(utils::head(colnames(expr_mat), 10))


# ============================================================================
# 5) Optional: define sample groups for demo comparisons
# ----------------------------------------------------------------------------
# This section is ONLY needed if you want to create a demo ranking vector
# for FGSEA mode.
#
# You should adapt the grep patterns to your own sample naming scheme.
# ============================================================================

dmso_samples <- grep("DMSO", colnames(expr_mat), value = TRUE)
temsi_50_samples <- grep("Temsi_50", colnames(expr_mat), value = TRUE)
temsi_500_samples <- grep("Temsi_500", colnames(expr_mat), value = TRUE)

cat("\n=== Sample groups detected ===\n")
cat("DMSO samples:", length(dmso_samples), "\n")
cat("Temsi 50 samples:", length(temsi_50_samples), "\n")
cat("Temsi 500 samples:", length(temsi_500_samples), "\n")

print(dmso_samples)
print(temsi_50_samples)
print(temsi_500_samples)


# ============================================================================
# 6) compass_gsc() examples
# ----------------------------------------------------------------------------
# PURPOSE:
#   Build context-specific reference gene sets from local subset files.
#
# IMPORTANT:
#   This is an upstream/reference-building function.
#   It does NOT score the user matrix yet.
# ============================================================================

cat("\n=== compass_gsc(): basic list output ===\n")

gsc_list_glioma <- compass_gsc(
  context = "glioma",
  n = 250,
  min_conf = 1,
  output = "list"
)

cat("Number of glioma gene sets:", length(gsc_list_glioma), "\n")
cat("First set names:\n")
print(utils::head(names(gsc_list_glioma), 5))


cat("\n=== compass_gsc(): higher-confidence only ===\n")

gsc_list_glioma_c23 <- compass_gsc(
  context = "glioma",
  n = 250,
  min_conf = 2,
  output = "list"
)

cat("Number of glioma gene sets with min_conf >= 2:", length(gsc_list_glioma_c23), "\n")


cat("\n=== compass_gsc(): target-filtered example ===\n")

gsc_targets_example <- compass_gsc(
  context = "glioma",
  n = 250,
  min_conf = 1,
  targets = c("AKT1", "BRAF", "BRD4", "CCNA2"),
  output = "list"
)

cat("Target-filtered gene sets:\n")
print(names(gsc_targets_example))


cat("\n=== compass_gsc(): data.frame output ===\n")

gsc_df_example <- compass_gsc(
  context = "glioma",
  n = 250,
  min_conf = 1,
  targets = c("AKT1", "BRAF"),
  output = "df"
)

print(dim(gsc_df_example))
print(utils::head(gsc_df_example[, 1:min(2, ncol(gsc_df_example)), drop = FALSE], 10))


cat("\n=== compass_gsc(): GeneSetCollection output ===\n")

gsc_gsc_example <- compass_gsc(
  context = "glioma",
  n = 250,
  min_conf = 1,
  targets = c("AKT1", "BRAF"),
  output = "gsc",
  meta_cols = "all"
)

print(gsc_gsc_example)


# ============================================================================
# 7) compass() in GSVA mode
# ----------------------------------------------------------------------------
# INPUT TYPE:
#   numeric expression matrix
#   rows = genes
#   columns = samples
#
# This is the most user-friendly current mode.
# ============================================================================

cat("\n=== compass(): GSVA mode, basic run ===\n")

compass_res_gsva <- compass(
  input = expr_mat,
  context = "glioma",
  mode = "gsva",
  conf = 1,
  gs.size = 250
)

cat("GSVA result dimensions:", dim(compass_res_gsva)[1], "x", dim(compass_res_gsva)[2], "\n")
print(compass_res_gsva[1:5, 1:min(5, ncol(compass_res_gsva)), drop = FALSE])


cat("\n=== compass(): GSVA mode, higher-confidence only ===\n")

compass_res_gsva_c23 <- compass(
  input = expr_mat,
  context = "glioma",
  mode = "gsva",
  conf = 2,
  gs.size = 250
)

cat("GSVA c2/c3 result dimensions:", dim(compass_res_gsva_c23)[1], "x", dim(compass_res_gsva_c23)[2], "\n")


cat("\n=== compass(): GSVA mode, scaled scores ===\n")

compass_res_gsva_scaled <- compass(
  input = expr_mat,
  context = "glioma",
  mode = "gsva",
  conf = 1,
  gs.size = 250,
  scale = TRUE
)

cat("Scaled GSVA result dimensions:", dim(compass_res_gsva_scaled)[1], "x", dim(compass_res_gsva_scaled)[2], "\n")


cat("\n=== compass(): GSVA mode, return gene sets used ===\n")

compass_res_gsva_with_sets <- compass(
  input = expr_mat,
  context = "glioma",
  mode = "gsva",
  conf = 1,
  gs.size = 250,
  return_gene_sets = TRUE
)

cat("Returned object names:\n")
print(names(compass_res_gsva_with_sets))
cat("Gene sets used:", length(compass_res_gsva_with_sets$gene_sets), "\n")


# ============================================================================
# 8) compass() in FGSEA mode
# ----------------------------------------------------------------------------
# INPUT TYPE:
#   named numeric vector
#   names = gene symbols
#   values = signed gene-level ranking statistic
#
# IMPORTANT:
#   This section uses a SIMPLE DEMO ranking:
#     rowMeans(Temsi_500) - rowMeans(DMSO)
#
#   This is only a placeholder for testing the API.
#   It is NOT a formal differential expression pipeline.
# ============================================================================

if (length(dmso_samples) > 0L && length(temsi_500_samples) > 0L) {
  cat("\n=== compass(): FGSEA mode, demo ranking vector ===\n")
  
  stats_vec_demo <- rowMeans(expr_mat[, temsi_500_samples, drop = FALSE]) -
    rowMeans(expr_mat[, dmso_samples, drop = FALSE])
  
  stats_vec_demo <- stats_vec_demo[is.finite(stats_vec_demo) & !is.na(stats_vec_demo)]
  stats_vec_demo <- sort(stats_vec_demo, decreasing = TRUE)
  
  cat("Length of demo ranking vector:", length(stats_vec_demo), "\n")
  print(utils::head(stats_vec_demo, 10))
  
  compass_res_fgsea <- compass(
    input = stats_vec_demo,
    context = "glioma",
    mode = "fgsea",
    conf = 1,
    gs.size = 250,
    min_size = 10,
    max_size = 500,
    n_perm_simple = 5000
  )
  
  cat("FGSEA result dimensions:", dim(compass_res_fgsea)[1], "x", dim(compass_res_fgsea)[2], "\n")
  print(utils::head(compass_res_fgsea, 10))
} else {
  cat("\n=== FGSEA demo skipped ===\n")
  cat("Reason: could not detect both DMSO and Temsi_500 samples from the current column names.\n")
}


# ============================================================================
# 9) Choose one sample of interest for summary examples
# ----------------------------------------------------------------------------
# Change this manually if you want to focus on another sample.
# ============================================================================

sample_of_interest <- "1_BN91_DMSO_S1_L001.dedup"

if (!sample_of_interest %in% colnames(compass_res_gsva)) {
  stop(
    "Selected `sample_of_interest` not found in GSVA result columns.\n",
    "Please adapt `sample_of_interest`.",
    call. = FALSE
  )
}

cat("\n=== Sample of interest ===\n")
cat(sample_of_interest, "\n")


# ============================================================================
# 10) compass_summary(): basic sample summary
# ----------------------------------------------------------------------------
# PURPOSE:
#   Organize COMPASS results into readable blocks.
#
# This is the simplest meaningful summary for one sample.
# ============================================================================

cat("\n=== compass_summary(): basic sample summary ===\n")

summary_basic <- compass_summary(
  compass_res = compass_res_gsva,
  top_n = 30,
  print_n = 10,
  sample = sample_of_interest,
  include_rna_concordance = FALSE,
  include_variability = FALSE,
  include_signature_consistency = FALSE,
  include_hypothesis_shortlist = FALSE,
  verbose = TRUE
)

cat("summary_basic dimensions:", dim(summary_basic)[1], "x", dim(summary_basic)[2], "\n")


# ============================================================================
# 11) compass_summary(): with RNA–activity discordance
# ----------------------------------------------------------------------------
# REQUIRES:
#   expr_mat
#
# This compares activity-like behavior with RNA behavior across samples.
# ============================================================================

cat("\n=== compass_summary(): with RNA–activity discordance ===\n")

summary_discordance <- compass_summary(
  compass_res = compass_res_gsva,
  top_n = 50,
  print_n = 10,
  sample = sample_of_interest,
  expr_mat = expr_mat,
  include_rna_concordance = TRUE,
  include_variability = FALSE,
  include_signature_consistency = FALSE,
  include_hypothesis_shortlist = FALSE,
  verbose = TRUE
)

cat("summary_discordance dimensions:", dim(summary_discordance)[1], "x", dim(summary_discordance)[2], "\n")


# ============================================================================
# 12) compass_summary(): with consistency / divergence warnings
# ----------------------------------------------------------------------------
# REQUIRES:
#   sample to be specified
#
# This adds:
#   - top convergent proteins
#   - top divergent proteins (warning)
# ============================================================================

cat("\n=== compass_summary(): with signature consistency ===\n")

summary_consistency <- compass_summary(
  compass_res = compass_res_gsva,
  top_n = 50,
  print_n = 10,
  sample = sample_of_interest,
  expr_mat = expr_mat,
  include_rna_concordance = TRUE,
  include_variability = FALSE,
  include_signature_consistency = TRUE,
  include_hypothesis_shortlist = FALSE,
  verbose = TRUE
)

cat("summary_consistency dimensions:", dim(summary_consistency)[1], "x", dim(summary_consistency)[2], "\n")


# ============================================================================
# 13) compass_summary(): full summary with four-field hypothesis shortlist
# ----------------------------------------------------------------------------
# This is the current most feature-rich summary mode.
#
# The final four-field hypothesis matrix uses:
#   hypothesis_score = |signal| * conf_weight * consistency_weight * support_weight
#
# It is NOT a truth statement.
# It is a transparent prioritization heuristic for manual follow-up.
# ============================================================================

cat("\n=== compass_summary(): full summary with hypothesis shortlist ===\n")

summary_full <- compass_summary(
  compass_res = compass_res_gsva,
  top_n = 100,
  print_n = 20,
  sample = sample_of_interest,
  expr_mat = expr_mat,
  include_rna_concordance = TRUE,
  include_variability = FALSE,
  include_signature_consistency = TRUE,
  include_hypothesis_shortlist = TRUE,
  shortlist_n = 3,
  verbose = TRUE
)

cat("summary_full dimensions:", dim(summary_full)[1], "x", dim(summary_full)[2], "\n")


# ============================================================================
# 14) compass_summary(): matrix-level view instead of one sample
# ----------------------------------------------------------------------------
# If sample = NULL, the function summarizes the whole matrix object.
#
# Useful especially for:
#   - variability overview
#   - non-sample-specific inspection
#
# NOTE:
#   For the console printout, a single sample is often more readable.
# ============================================================================

cat("\n=== compass_summary(): matrix-level / all-sample view ===\n")

summary_matrix_level <- compass_summary(
  compass_res = compass_res_gsva,
  top_n = 30,
  print_n = 10,
  sample = NULL,
  include_rna_concordance = FALSE,
  include_variability = TRUE,
  include_signature_consistency = FALSE,
  include_hypothesis_shortlist = FALSE,
  verbose = TRUE
)

cat("summary_matrix_level dimensions:", dim(summary_matrix_level)[1], "x", dim(summary_matrix_level)[2], "\n")


# ============================================================================
# 15) Optional GO annotation
# ----------------------------------------------------------------------------
# This uses the prebuilt annotation file:
#   annotation/compass_go_annotations.rds
#
# COMMENT:
#   This is currently better treated as a downstream utility than as the core
#   interpretation engine of COMPASS summary.
# ============================================================================

if (exists("compass_annotate") && file.exists(go_annotation_path)) {
  cat("\n=== compass_annotate(): annotate full summary ===\n")
  
  summary_full_annot <- compass_annotate(summary_full)
  
  cat("Annotated summary dimensions:", dim(summary_full_annot)[1], "x", dim(summary_full_annot)[2], "\n")
  
  print(
    utils::head(
      summary_full_annot[, c(
        "protein",
        "n_go_total",
        "n_go_bp",
        "n_go_mf",
        "n_go_cc",
        "go_bp_terms"
      )],
      10
    )
  )
  
  cat("\n=== compass_annotate(): BP-only example ===\n")
  
  summary_full_annot_bp <- compass_annotate(
    summary_full,
    namespaces = "biological_process"
  )
  
  print(
    utils::head(
      summary_full_annot_bp[, c("protein", "n_go_total", "n_go_bp", "go_bp_terms")],
      10
    )
  )
  
  cat("\n=== compass_annotate(): stricter evidence example ===\n")
  
  summary_full_annot_strict <- compass_annotate(
    summary_full,
    evidence_codes = c("IDA", "IMP", "IGI", "IPI", "TAS")
  )
  
  print(
    utils::head(
      summary_full_annot_strict[, c("protein", "n_go_total", "go_bp_terms")],
      10
    )
  )
} else {
  cat("\n=== GO annotation section skipped ===\n")
  cat("Reason: `compass_annotate()` or `compass_go_annotations.rds` not available.\n")
}


# ============================================================================
# 16) Save selected outputs
# ----------------------------------------------------------------------------
# Adjust this section as needed.
# ============================================================================

saveRDS(compass_res_gsva, file.path(results_dir, "compass_res_gsva_glioma.rds"))
saveRDS(summary_full, file.path(results_dir, "summary_full_glioma.rds"))

if (exists("summary_full_annot")) {
  saveRDS(summary_full_annot, file.path(results_dir, "summary_full_glioma_annotated.rds"))
}

if (exists("compass_res_fgsea")) {
  saveRDS(compass_res_fgsea, file.path(results_dir, "compass_res_fgsea_glioma_demo.rds"))
}

cat("\n=== Saved selected outputs to ===\n")
cat(results_dir, "\n")


# ============================================================================
# 17) Notes / caveats
# ----------------------------------------------------------------------------
# - GSVA mode:
#     input = expression matrix
#     output = signature x sample score matrix
#
# - FGSEA mode:
#     input = named numeric ranking vector
#     output = NES-style data.frame
#
# - AUCell mode:
#     not implemented yet
#
# - RNA–activity discordance:
#     requires overlapping samples between compass_res and expr_mat
#
# - hypothesis shortlist:
#     heuristic prioritization layer, not final truth
#
# - GO annotation:
#     currently useful as infrastructure / optional downstream utility,
#     not yet the main interpretation engine
# ============================================================================

cat("\n=== Workflow test completed successfully ===\n")