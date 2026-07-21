test_that("temporary output paths are safe for CRAN", {
  out_dir <- file.path(tempdir(), "biofetchR_test_output")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  expect_true(dir.exists(out_dir))
  expect_true(startsWith(
    normalizePath(out_dir, winslash = "/", mustWork = FALSE),
    normalizePath(tempdir(), winslash = "/", mustWork = FALSE)
  ))
})
