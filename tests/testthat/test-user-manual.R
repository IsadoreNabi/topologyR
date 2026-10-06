# ======================================================================
# The user manual is executable (0.4.0)
#
# Up to 0.3.0 the manual shipped in inst/manuals was never run, and it
# drifted from the package: it still documented a field removed in
# 0.3.0 and printed outputs written by hand. Its code is now extracted
# and run against the installed package at every check.
# ======================================================================

test_that("the code of the user manual runs with the installed package", {
  skip_if_not_installed("knitr")
  rmd <- system.file("manuals", "user_manual.Rmd", package = "topologyR")
  skip_if(!nzchar(rmd), "the manual is not installed")
  code <- tempfile(fileext = ".R")
  on.exit(unlink(code), add = TRUE)
  knitr::purl(rmd, output = code, quiet = TRUE, documentation = 0L)
  env <- new.env(parent = globalenv())
  expect_no_error(utils::capture.output(sys.source(code, envir = env)))
  expect_identical(env$d, -env$d_rev)
  expect_identical(env$inv$forward_base_size, env$inv$backward_base_size)
})
