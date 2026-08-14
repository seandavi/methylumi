## Plot functions: smoke tests only. These exist because a defunct lattice call
## sat in parplot() undetected for years and eventually broke installation of
## the whole package (#28). Rendering is not compared -- the point is that every
## plotting entry point is executed by something.

## Draw to a null device for the whole file. setup()/teardown() are deprecated
## in testthat's 3rd edition, and every test here needs the same device, so open
## it once at file scope.
pdf(NULL)

test_that("parplot works (regression test for #28)", {
  ## parplot() called lattice::parallel, which went deprecated -> defunct ->
  ## unexported. Because the import resolved at lazy-load time, this took the
  ## whole package down, not just this function.
  data(mldat, package = "methylumi", envir = environment())
  p <- methylumi:::parplot(mldat)
  expect_s3_class(p, "trellis")
  expect_silent(print(p))
})

test_that("parplot accepts the documented quantiles and what arguments", {
  data(mldat, package = "methylumi", envir = environment())
  expect_s3_class(methylumi:::parplot(mldat, quantiles = seq(0, 1, 0.5)), "trellis")
  expect_s3_class(methylumi:::parplot(mldat, what = "pvals"), "trellis")
})

test_that("qcplot runs on a MethyLumiSet and on its QC data", {
  data(mldat, package = "methylumi", envir = environment())
  expect_no_error(qcplot(mldat, "NEGATIVE"))
  expect_no_error(qcplot(QCdata(mldat), "NEGATIVE"))
})

test_that("corplot and plotSampleIntensities run", {
  data(mldat, package = "methylumi", envir = environment())
  expect_no_error(corplot(mldat))
  expect_no_error(plotSampleIntensities(mldat))
})

test_that("qc.probe.plot runs on IDAT-derived data", {
  expect_no_error(print(qc.probe.plot(example_idats())))
})

test_that("the ggplot2 plots build without deprecation warnings (#40)", {
  ## These used qplot(), deprecated in ggplot2 3.4.0. lifecycle only warns once
  ## per session by default, so force it on for the duration of the test.
  old <- options(lifecycle_verbosity = "warning")
  on.exit(options(old), add = TRUE)
  data(mldat, package = "methylumi", envir = environment())
  for (p in list(qc.probe.plot(example_idats()),
                 plotNAs(mldat),
                 plotProbeNAs(mldat))) {
    expect_no_warning(expect_no_error(ggplot2::ggplot_build(p)))
  }
})

test_that("plotNegOob runs (regression test for #36)", {
  ## scale_y_continuous(breaks = NA) errored under current ggplot2, which wants
  ## NULL to mean "draw no breaks".
  expect_no_error(plotNegOob(example_idats()))
})
