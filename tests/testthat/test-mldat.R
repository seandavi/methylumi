## The bundled GoldenGate MethyLumiSet: accessors, subsetting, normalization,
## coercion and filtering.

data(mldat, package = "methylumi", envir = environment())

test_that("mldat loads with the expected shape", {
  expect_s4_class(mldat, "MethyLumiSet")
  expect_equal(unname(dim(mldat)), c(1536L, 10L))
  expect_equal(head(sampleNames(mldat), 3), c("M_1", "M_2", "M_3"))
  expect_equal(head(featureNames(mldat), 3),
               c("AATK_E63_R", "AATK_P519_R", "AATK_P709_R"))
})

test_that("assay accessors return the expected values", {
  expect_equal(sum(betas(mldat)),        5691.890947, tolerance = TOL)
  expect_equal(sum(pvals(mldat)),         473.826186, tolerance = TOL)
  expect_equal(sum(methylated(mldat)),   74208765.64, tolerance = TOL)
  expect_equal(sum(unmethylated(mldat)), 95192882.17, tolerance = TOL)

  for (f in list(betas, pvals, methylated, unmethylated)) {
    expect_equal(unname(dim(f(mldat))), c(1536L, 10L))
    expect_false(anyNA(f(mldat)))
  }
})

test_that("betas are bounded and pvals are probabilities", {
  b <- betas(mldat)
  expect_true(all(b >= 0 & b <= 1))
  p <- pvals(mldat)
  expect_true(all(p >= 0 & p <= 1))
})

test_that("mldat copies are independent (regression test for #35)", {
  ## The shipped mldat.rda used to have storageMode == "environment" rather
  ## than Biobase's default "lockedEnvironment", so its assayData was a bare
  ## environment shared between every copy and modifying a copy silently
  ## mutated the original.
  expect_equal(storageMode(mldat), "lockedEnvironment")

  x <- mldat
  before <- sum(methylated(mldat))
  methylated(x) <- methylated(x) + 1
  expect_equal(sum(methylated(mldat)), before)     # original untouched
  expect_false(sum(methylated(x)) == before)       # the copy did change
})

test_that("replacement methods round-trip", {
  x <- mldat
  b <- betas(x)
  betas(x) <- b * 0 + 0.5
  expect_true(all(betas(x) == 0.5))
  betas(x) <- b
  expect_equal(betas(x), b)

  m <- methylated(x)
  methylated(x) <- m + 1
  expect_equal(methylated(x), m + 1)
  expect_equal(sum(methylated(mldat)), sum(m))     # original untouched
})

test_that("QC data is present and has the documented control types", {
  qc <- QCdata(mldat)
  expect_s4_class(qc, "MethyLumiQC")
  expect_equal(unname(dim(qc)), c(44L, 10L))
  expect_equal(controlTypes(mldat),
               c("ALLELE SPECIFIC EXTENSION", "EXTENSION GAP",
                 "FIRST HYBRIDIZATION", "GENDER", "NEGATIVE",
                 "PCR CONTAMINATION", "SECOND HYBRIDIZATION"))
})

test_that("subsetting selects the right features and samples", {
  s <- mldat[1:100, 1:3]
  expect_equal(unname(dim(s)), c(100L, 3L))
  expect_equal(sum(betas(s)), 117.900646, tolerance = TOL)
  expect_equal(betas(s), betas(mldat)[1:100, 1:3])
  expect_equal(sampleNames(s), head(sampleNames(mldat), 3))

  expect_equal(unname(dim(mldat[, 1:3])), c(1536L, 3L))
  expect_equal(unname(dim(mldat[1:100, ])), c(100L, 10L))
})

test_that("subsetting keeps QC data (regression test for #21)", {
  ## eSet's "[" drops the QC slot, so the method saves it before calling
  ## callNextMethod() and reattaches it afterwards. QC probes are control
  ## probes in their own feature space, so feature indices must not touch them
  ## -- only sample selection applies.
  both <- QCdata(mldat[1:100, 1:3])
  expect_s4_class(both, "MethyLumiQC")
  expect_equal(unname(dim(both)), c(44L, 3L))
  expect_equal(sampleNames(both), head(sampleNames(mldat), 3))

  samples <- QCdata(mldat[, 1:3])
  expect_s4_class(samples, "MethyLumiQC")
  expect_equal(unname(dim(samples)), c(44L, 3L))

  ## Feature-only subsetting leaves QC entirely alone.
  feats <- QCdata(mldat[1:100, ])
  expect_s4_class(feats, "MethyLumiQC")
  expect_equal(unname(dim(feats)), c(44L, 10L))

  ## And the QC values themselves are the right columns, not just the right
  ## shape. (MethyLumiQC has no betas method; use an assay it does carry.)
  expect_equal(methylated(QCdata(mldat))[, 1:3], methylated(both))
})

test_that("subsetting records history", {
  h <- getHistory(mldat[1:100, 1:3])
  expect_true(any(grepl("Subset of 100 features & 3 samples", h$command)))
})

test_that("normalizeMethyLumiSet is unchanged", {
  n <- normalizeMethyLumiSet(mldat)
  expect_s4_class(n, "MethyLumiSet")
  expect_equal(unname(dim(n)), c(1536L, 10L))
  expect_equal(sum(betas(n), na.rm = TRUE), 5517.872753, tolerance = TOL)
})

test_that("noob background correction refuses non-IDAT input (#24)", {
  ## mldat has no methylated.OOB/unmethylated.OOB -- like any set built from a
  ## GEO series matrix or GenomeStudio output. It used to die inside
  ## intensities.OOB() ("unable to find an inherited method ... for signature
  ## MethyLumiQC, missing"); now it says what is missing and what to do.
  expect_error(methylumi.bgcorr(mldat), "out-of-band")
  expect_error(methylumi.bgcorr(mldat), "methylumIDAT")
  expect_error(methylumi.bgcorr(mldat, method = "mode"), "out-of-band")

  ## The signature reported in #24: a MethyLumiQC reaching methylumi.bgcorr had
  ## no intensities.OOB method at all, so dispatch failed before any check.
  expect_error(methylumi.bgcorr(QCdata(mldat)), "out-of-band")
})

test_that("MethyLumiSet coerces to MethyLumiM", {
  mm <- as(mldat, "MethyLumiM")
  expect_s4_class(mm, "MethyLumiM")
  expect_equal(unname(dim(mm)), c(1536L, 10L))
  expect_equal(sum(exprs(mm), na.rm = TRUE), -14514.904905, tolerance = TOL)
})

test_that("the MethyLumiM round-trip is lossy, by design", {
  ## Betas -> M-values -> betas does not return the original values, because
  ## estimateM applies an offset. Pinned so the size of the discrepancy cannot
  ## drift unnoticed.
  back <- as(as(mldat, "MethyLumiM"), "MethyLumiSet")
  expect_equal(sum(betas(back), na.rm = TRUE), 5761.247359, tolerance = TOL)
})

test_that("varFilter keeps half the features", {
  vf <- varFilter(mldat)
  expect_equal(unname(nrow(vf$eset)), 768L)
  expect_true(all(featureNames(vf$eset) %in% featureNames(mldat)))
})

test_that("combine of disjoint sample sets restores the whole", {
  a <- mldat[, 1:5]
  b <- mldat[, 6:10]
  ab <- combine(a, b)
  expect_equal(unname(dim(ab)), unname(dim(mldat)))
  expect_equal(betas(ab)[, sampleNames(mldat)], betas(mldat))
})

test_that("combine keeps control probes now that subsetting does (#21)", {
  ## combine() used to emit "Dropped control probes: any(is.null(QCdata(x),
  ## QCdata(y))) == TRUE" for any combine of subsets, because subsetting had
  ## already thrown the QC data away.
  ab <- combine(mldat[, 1:5], mldat[, 6:10])
  qc <- QCdata(ab)
  expect_s4_class(qc, "MethyLumiQC")
  expect_equal(unname(dim(qc)), c(44L, 10L))
})
