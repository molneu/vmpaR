make_overlap_gene_sets <- function() {
  gene_sets <- list(
    ONE_c1 = "GENE1",
    TEN_c1 = paste0("GENE", 1:10),
    ZERO_c1 = paste0("OTHER", 1:10)
  )
  attr(gene_sets, "signature_metadata") <- data.frame(
    gene_set_name = names(gene_sets),
    stringsAsFactors = FALSE
  )
  gene_sets
}

test_that("backend-specific minimum-size defaults preserve prior behavior", {
  expect_identical(formals(vmpa)$gsva_min_size, 1L)
  expect_identical(formals(vmpa)$fgsea_min_size, 10L)
})

test_that("valid unique human gene symbols are accepted", {
  expr <- matrix(
    1:6,
    nrow = 3,
    dimnames = list(c("TP53", "EGFR", "AKT1"), c("S1", "S2"))
  )
  stats <- stats::setNames(c(3, 2, 1), c("TP53", "EGFR", "AKT1"))

  expect_identical(vmpaR:::.vmpa_prepare_gsva_input(expr), expr)
  expect_identical(names(vmpaR:::.vmpa_prepare_fgsea_input(stats)), names(stats))
})

test_that("missing, blank, padded, and duplicated input symbols are rejected", {
  invalid_ids <- list(
    c("TP53", NA_character_),
    c("TP53", ""),
    c("TP53", "   "),
    c("TP53", " EGFR"),
    c("TP53", "TP53")
  )

  for (ids in invalid_ids) {
    expect_error(vmpaR:::.vmpa_validate_gene_ids(ids, "test input"))
  }

  expect_error(
    vmpaR:::.vmpa_prepare_gsva_input(matrix(1:4, nrow = 2)),
    "must contain gene symbols"
  )
  expect_error(
    vmpaR:::.vmpa_prepare_fgsea_input(c(2, 1)),
    "must contain gene symbols"
  )
})

test_that("input identifiers are rejected before reference resolution", {
  invalid_expr <- matrix(
    1:4,
    nrow = 2,
    dimnames = list(c("TP53", "TP53"), c("S1", "S2"))
  )

  expect_error(
    testthat::with_mocked_bindings(
      vmpa(
        input = invalid_expr,
        context = "glioma",
        algorithm = "gsva",
        verbose = FALSE
      ),
      .vmpa_resolve_subset_file = function(...) {
        stop("reference resolution should not be reached")
      },
      .package = "vmpaR"
    ),
    "duplicated gene symbols"
  )
})

test_that("GSVA and FGSEA use separate overlap thresholds", {
  gene_sets <- make_overlap_gene_sets()
  input_genes <- paste0("GENE", 1:10)

  gsva_message <- capture_messages(
    gsva_sets <- vmpaR:::.vmpa_filter_gene_sets_by_overlap(
      gene_sets, input_genes, min_size = 1L, backend = "gsva", verbose = TRUE
    )
  )
  fgsea_message <- capture_messages(
    fgsea_sets <- vmpaR:::.vmpa_filter_gene_sets_by_overlap(
      gene_sets, input_genes, min_size = 10L, backend = "fgsea", verbose = TRUE
    )
  )

  expect_named(gsva_sets, c("ONE_c1", "TEN_c1"))
  expect_named(fgsea_sets, "TEN_c1")
  expect_true(any(grepl("2 of 3 selected VMPA gene sets", gsva_message)))
  expect_true(any(grepl("1 of 3 selected VMPA gene sets", fgsea_message)))
  expect_identical(
    attr(gsva_sets, "signature_metadata")$input_overlap,
    c(1L, 10L)
  )
})

test_that("a complete identifier mismatch stops before the backend", {
  expect_error(
    vmpaR:::.vmpa_filter_gene_sets_by_overlap(
      make_overlap_gene_sets(),
      c("ENSG00000141510", "ENSG00000146648"),
      min_size = 1L,
      backend = "gsva",
      verbose = FALSE
    ),
    "Ensembl and Entrez identifiers are not converted automatically",
    fixed = TRUE
  )
})

test_that("insufficient per-set overlap stops before the backend", {
  expect_error(
    vmpaR:::.vmpa_filter_gene_sets_by_overlap(
      make_overlap_gene_sets(),
      paste0("GENE", 1:9),
      min_size = 10L,
      backend = "fgsea",
      verbose = FALSE
    ),
    "No selected VMPA gene set has at least 10 represented input genes",
    fixed = TRUE
  )
})

test_that("FGSEA overlap is calculated after non-finite statistics are removed", {
  stats <- stats::setNames(seq_len(10), paste0("GENE", 1:10))
  stats[["GENE1"]] <- NA_real_
  prepared <- vmpaR:::.vmpa_prepare_fgsea_input(stats)

  expect_false("GENE1" %in% names(prepared))
  expect_error(
    vmpaR:::.vmpa_filter_gene_sets_by_overlap(
      list(TEN_c1 = paste0("GENE", 1:10)),
      names(prepared),
      min_size = 10L,
      backend = "fgsea",
      verbose = FALSE
    ),
    "No selected VMPA gene set has at least 10 represented input genes",
    fixed = TRUE
  )
})

test_that("GSVA backend uses its mapped gene sets at minSize 1", {
  expr <- outer(1:5, 1:4, function(gene, sample) gene * sample + sample^2)
  dimnames(expr) <- list(paste0("GENE", 1:5), paste0("S", 1:4))
  gene_sets <- vmpaR:::.vmpa_filter_gene_sets_by_overlap(
    make_overlap_gene_sets(),
    rownames(expr),
    min_size = 1L,
    backend = "gsva",
    verbose = FALSE
  )

  result <- suppressWarnings(vmpaR:::.vmpa_run_gsva(
    expr_mat = expr,
    gene_set_list = gene_sets,
    min_size = 1L,
    score_scaling = "none",
    verbose = FALSE
  ))

  expect_setequal(rownames(result), c("ONE_c1", "TEN_c1"))
  used_gene_sets <- vmpaR:::.vmpa_subset_gene_sets_to_backend(
    gene_sets,
    rownames(result)
  )
  expect_equal(length(used_gene_sets), nrow(result))
  expect_identical(
    attr(used_gene_sets, "signature_metadata")$gene_set_name,
    rownames(result)
  )
})

test_that("FGSEA excludes sets below ten represented genes", {
  stats <- seq(2, -2, length.out = 20L)
  names(stats) <- paste0("GENE", 1:20)
  gene_sets <- list(
    NINE_c1 = paste0("GENE", 1:9),
    TEN_c1 = paste0("GENE", 1:10)
  )

  result <- vmpaR:::.vmpa_run_fgsea(
    stats_vec = stats,
    gene_set_list = gene_sets,
    context = "glioma",
    gs_size = 10L,
    conf = 1L,
    min_size = 10L,
    max_size = 20L,
    n_perm_simple = 100L,
    seed = 123L
  )

  expect_identical(result$pathway, "TEN_c1")
  expect_identical(result$size, 10L)
})
