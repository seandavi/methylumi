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

test_that("mldat shares assayData by reference (BUG, see #35)", {
  ## Characterization test. The shipped mldat.rda has
  ## storageMode == "environment" rather than Biobase's default
  ## "lockedEnvironment", so its assayData is a bare environment shared between
  ## every copy of the object. R's copy-on-modify does not protect the payload:
  ##
  ##     x <- mldat
  ##     methylated(x) <- methylated(x) + 1   # also modifies mldat
  ##
  ## methylumiR() and methylumIDAT() both produce "lockedEnvironment" objects,
  ## so this is confined to the saved example dataset -- but that dataset is
  ## what every reader of the vignette experiments with.
  ##
  ## Invert this when #35 is fixed.
  expect_equal(storageMode(mldat), "environment")

  x <- mldat
  before <- sum(methylated(mldat))
  methylated(x) <- methylated(x) + 1
  expect_false(sum(methylated(mldat)) == before)   # the original was mutated

  methylated(x) <- methylated(x) - 1               # undo, for later tests
  expect_equal(sum(methylated(mldat)), before)
})

test_that("replacement methods round-trip", {
  ## Work on a properly isolated copy so this test cannot leak into the others,
  ## which is precisely what #35 causes.
  x <- mldat
  storageMode(x) <- "lockedEnvironment"

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

test_that("subsetting drops QC data (BUG, see #21)", {
  ## Characterization test: this documents current behaviour, it does not
  ## endorse it. eSet's "[" method drops the QC slot in callNextMethod(), and
  ## the guard added for #21 then runs
  ##
  ##     x@QC <- x@QC[, j, drop = FALSE]
  ##
  ## against an already-NULL slot. NULL[, j, drop = FALSE] silently returns
  ## NULL rather than erroring, so the fix is a no-op and QC data is lost on
  ## every form of subsetting.
  ##
  ## When #21 is actually fixed, these expectations must be inverted to assert
  ## that QC data survives, subset to the selected samples.
  expect_null(QCdata(mldat[1:100, 1:3]))
  expect_null(QCdata(mldat[, 1:3]))
  expect_null(QCdata(mldat[1:100, ]))
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
