## Raw IDAT import: methylumIDAT() on the three bundled HumanMethylation27
## IDAT pairs, plus background correction and the minfi coercion.

test_that("barcodes are discovered from the bundled IDATs", {
  expect_equal(idat_barcodes(),
               c("5318317007_A", "5318317007_B", "5318317007_C"))
})

test_that("methylumIDAT imports the bundled 27k IDATs", {
  mi <- example_idats()
  expect_s4_class(mi, "MethyLumiSet")
  expect_equal(unname(dim(mi)), c(27578L, 3L))
  expect_equal(annotation(mi), "IlluminaHumanMethylation27k")
  expect_equal(sampleNames(mi), idat_barcodes())
  expect_equal(head(featureNames(mi), 2), c("cg00000292", "cg00002426"))
})

test_that("methylumIDAT assay values are unchanged", {
  mi <- example_idats()
  expect_equal(sum(betas(mi), na.rm = TRUE),        25691.821526, tolerance = TOL)
  expect_equal(sum(methylated(mi)),                    458162070, tolerance = TOL)
  expect_equal(sum(unmethylated(mi)),                 1171217850, tolerance = TOL)
  expect_equal(sum(pvals(mi)),                              0.25, tolerance = TOL)

  ## Exactly one probe has an undefined beta (zero total intensity).
  expect_equal(sum(is.na(betas(mi))), 1L)
})

test_that("methylumIDAT attaches QC and out-of-band data", {
  mi <- example_idats()
  expect_equal(unname(dim(QCdata(mi))), c(144L, 3L))
  expect_equal(sum(unlist(intensities.OOB(mi)), na.rm = TRUE),
               167458050, tolerance = TOL)
})

test_that("the parallel path gives the same answer as the serial path", {
  ## This is the guard for replacing the .mclapply shim in #30: whatever the
  ## parallel backend becomes, it must not change the numbers.
  mi <- example_idats()
  par <- suppressMessages(
    methylumIDAT(idat_barcodes(), idatPath = idat_path(), parallel = TRUE)
  )
  expect_equal(betas(par),        betas(mi))
  expect_equal(methylated(par),   methylated(mi))
  expect_equal(unmethylated(par), unmethylated(mi))
})

test_that("noob background correction is unchanged", {
  bg <- suppressMessages(methylumi.bgcorr(example_idats()))
  expect_s4_class(bg, "MethyLumiSet")
  expect_equal(sum(betas(bg), na.rm = TRUE), 25441.561351, tolerance = TOL)
  expect_equal(sum(is.na(betas(bg))), 3L)
})

test_that("stripOOB drops the out-of-band assays without touching betas", {
  mi <- example_idats()
  so <- stripOOB(mi)
  expect_equal(sort(assayDataElementNames(so)),
               c("betas", "methylated", "pvals", "unmethylated"))
  expect_equal(betas(so), betas(mi))
})

test_that("MethyLumiSet coerces to a minfi MethylSet", {
  ms <- as(example_idats(), "MethylSet")
  expect_s4_class(ms, "MethylSet")
  expect_equal(unname(dim(ms)), c(27578L, 3L))
  expect_equal(sum(minfi::getBeta(ms), na.rm = TRUE), 25692.398717, tolerance = TOL)
})

test_that("IDATsToMatrices reads both channels", {
  mats <- suppressMessages(
    IDATsToMatrices(idat_barcodes(), idatPath = idat_path())
  )
  expect_length(mats, 3L)
  expect_named(mats, idat_barcodes())
})
