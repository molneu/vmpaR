test_that("the public API uses the VMPA names", {
  exports <- getNamespaceExports("vmpaR")
  retired_exports <- paste0("com", "pass", c("", "_gsc"))

  expect_true(all(c("vmpa", "vmpa_gsc") %in% exports))
  expect_false(any(retired_exports %in% exports))
})

test_that("VMPA results use the new S3 class", {
  result <- vmpaR:::.vmpa_new_result(
    data.frame(score = 1),
    context = "glioma",
    algorithm = "gsva",
    unique = TRUE
  )

  expect_s3_class(result, "vmpa_result")
  expect_match(capture.output(print(result))[1], "^VMPA result$")
})

test_that("the context registry points to the VMPA reference subsets", {
  registry <- vmpaR:::.vmpa_context_registry()

  expect_equal(nrow(registry), 10L)
  expect_setequal(
    registry$context,
    c(
      "breast", "crc", "gastric", "glioma", "headneck",
      "melanoma", "nsclc", "ovarian", "pdac", "prostate"
    )
  )
  expect_true(all(grepl(
    "^https://ndownloader[.]figshare[.]com/files/[0-9]+$",
    registry$subset_url
  )))
})
