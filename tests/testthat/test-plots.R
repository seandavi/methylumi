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

test_that("plotNegOob is broken by modern ggplot2 (BUG, see #36)", {
  ## Characterization test. R/plotNegOob.R:48 calls
  ##
  ##     scale_y_continuous(breaks = NA)
  ##
  ## ggplot2 used to tolerate NA there as "draw no breaks"; it now rejects it
  ## and asks for NULL. The function is exported, so this is a user-visible
  ## breakage, not an internal wart. One-word fix, deliberately not made here --
  ## #27 pins current behaviour and the fix lands separately.
  ##
  ## Invert this when #36 is fixed.
  expect_error(plotNegOob(example_idats()), "breaks")
})
