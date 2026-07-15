make_test_gct <- function(mat, symbols, ids = colnames(mat)) {
  methods::new(
    "GCT",
    mat = mat,
    rid = paste0("row", seq_len(nrow(mat))),
    cid = ids,
    rdesc = data.frame(symbol = symbols, stringsAsFactors = FALSE),
    cdesc = data.frame(
      id = ids,
      cmap_name = ids,
      cps_conf_total = rep(1L, length(ids)),
      cancer_driver_summary = rep("None", length(ids)),
      stringsAsFactors = FALSE
    )
  )
}

test_that("gene sets contain the lowest-valued finite unique genes", {
  mat <- matrix(
    c(2, 5, -2, -10, -9, NA, NaN, Inf, -5, -Inf),
    ncol = 1,
    dimnames = list(NULL, "SIG1")
  )
  symbols <- c("B", "A", "A", NA, "", "C", "D", "E", "F", "G")
  gct <- make_test_gct(mat, symbols)

  result <- vmpaR:::.vmpa_build_gene_sets(gct, n = 3L)
  repeated_result <- vmpaR:::.vmpa_build_gene_sets(gct, n = 3L)

  expect_identical(unname(result[[1]]), c("F", "A", "B"))
  expect_identical(repeated_result, result)
  expect_length(result[[1]], 3L)
  expect_false(anyNA(result[[1]]))
  expect_false(anyDuplicated(result[[1]]) > 0L)
})

test_that("invalid values and blank symbols cannot fill a gene set", {
  mat <- matrix(
    c(3, 1, 0, NA, NaN, Inf),
    ncol = 2,
    dimnames = list(NULL, c("VALID", "EMPTY"))
  )
  symbols <- c(" GENE2 ", "GENE1", " ")
  gct <- make_test_gct(mat, symbols)

  result <- vmpaR:::.vmpa_build_gene_sets(gct, n = 10L)

  expect_named(result, "VALID_c1")
  expect_identical(unname(result[[1]]), c("GENE1", "GENE2"))
  expect_identical(
    attr(result, "signature_metadata")$gene_set_name,
    "VALID_c1"
  )
})
