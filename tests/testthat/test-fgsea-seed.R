make_fgsea_seed_input <- function() {
  stats <- seq(2, -2, length.out = 100L)
  names(stats) <- paste0("GENE", seq_along(stats))

  list(
    stats = stats,
    pathways = list(
      SET_A_c1 = paste0("GENE", 1:25),
      SET_B_c1 = paste0("GENE", 26:50),
      SET_C_c1 = paste0("GENE", 51:75)
    )
  )
}

run_seed_test_fgsea <- function(seed) {
  x <- make_fgsea_seed_input()

  vmpaR:::.vmpa_run_fgsea(
    stats_vec = x$stats,
    gene_set_list = x$pathways,
    context = "glioma",
    gs_size = 25L,
    conf = 1L,
    min_size = 10L,
    max_size = 50L,
    n_perm_simple = 100L,
    seed = seed
  )
}

test_that("FGSEA is exactly reproducible with a fixed seed", {
  expect_identical(formals(vmpa)$seed, 123L)

  first <- run_seed_test_fgsea(123L)
  second <- run_seed_test_fgsea(123L)

  expect_identical(first, second)
})

test_that("FGSEA restores the user's global random-number state", {
  set.seed(987L)
  state_before <- .Random.seed

  invisible(run_seed_test_fgsea(123L))

  expect_identical(.Random.seed, state_before)
})

test_that("FGSEA does not leave a global seed when none existed", {
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    old_state <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    on.exit(assign(".Random.seed", old_state, envir = .GlobalEnv), add = TRUE)
    rm(".Random.seed", envir = .GlobalEnv)
  }

  invisible(run_seed_test_fgsea(123L))

  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
})

test_that("seed validation rejects invalid values", {
  invalid <- list(c(1, 2), NA_integer_, "123", 1.5, Inf, NaN)

  for (seed in invalid) {
    expect_error(
      vmpaR:::.vmpa_validate_seed(seed),
      "must be NULL or a single non-missing integer",
      fixed = TRUE
    )
  }
})

test_that("seed NULL runs FGSEA without package-controlled seeding", {
  expect_s3_class(run_seed_test_fgsea(NULL), "data.frame")
})

test_that("GSVA results are independent of the seed argument", {
  gene_sets <- list(SET_A_c1 = c("GENE1", "GENE2"))
  attr(gene_sets, "signature_metadata") <- data.frame(
    gene_set_name = "SET_A_c1",
    stringsAsFactors = FALSE
  )

  results <- testthat::with_mocked_bindings(
    {
      list(
        seed_1 = vmpa(
          input = matrix(1, nrow = 1, dimnames = list("GENE1", "SAMPLE1")),
          context = "glioma",
          algorithm = "gsva",
          unique = FALSE,
          seed = 1L,
          verbose = FALSE
        ),
        seed_2 = vmpa(
          input = matrix(1, nrow = 1, dimnames = list("GENE1", "SAMPLE1")),
          context = "glioma",
          algorithm = "gsva",
          unique = FALSE,
          seed = 999L,
          verbose = FALSE
        )
      )
    },
    .vmpa_resolve_subset_file = function(...) "synthetic.rds",
    .vmpa_read_subset_gct = function(...) NULL,
    .vmpa_build_gene_sets = function(...) gene_sets,
    .vmpa_run_gsva = function(...) {
      matrix(0.5, nrow = 1, dimnames = list("SET_A_c1", "SAMPLE1"))
    },
    .vmpa_format_gsva_result = function(score_mat, ...) score_mat,
    .package = "vmpaR"
  )

  expect_identical(results$seed_1, results$seed_2)
})
