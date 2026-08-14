## Text import: methylumiR() on the bundled GoldenGate sample and control files.

test_that("methylumiR imports the bundled sample file", {
  mr <- example_text()
  expect_s4_class(mr, "MethyLumiSet")
  expect_equal(unname(dim(mr)), c(1536L, 10L))
  expect_equal(head(sampleNames(mr), 3),
               c("1632405013_R006_C012", "1632405013_R007_C001",
                 "1632405013_R007_C002"))
})

test_that("methylumiR assay values are unchanged", {
  mr <- example_text()
  expect_equal(sum(betas(mr)),        5691.890947, tolerance = TOL)
  expect_equal(sum(pvals(mr)),         473.826186, tolerance = TOL)
  expect_equal(sum(methylated(mr)),   74208765.64, tolerance = TOL)
  expect_equal(sum(unmethylated(mr)), 95192882.17, tolerance = TOL)
  expect_false(anyNA(betas(mr)))
})

test_that("methylumiR reproduces the shipped mldat values", {
  ## mldat is the saved result of importing these same files, so the numbers
  ## must agree. Only the sample names differ: mldat carries the friendly
  ## M_1..M_10 labels, the fresh import carries the sentrix barcodes.
  data(mldat, package = "methylumi", envir = environment())
  mr <- example_text()
  expect_equal(unname(betas(mr)), unname(betas(mldat)), tolerance = TOL)
  expect_equal(unname(pvals(mr)), unname(pvals(mldat)), tolerance = TOL)
  expect_false(identical(sampleNames(mr), sampleNames(mldat)))
})

test_that("methylumiR attaches QC data from the controls file", {
  qc <- QCdata(example_text())
  expect_s4_class(qc, "MethyLumiQC")
  expect_equal(unname(ncol(qc)), 10L)
})

test_that("methylumiR works without a controls file", {
  mr <- suppressMessages(methylumiR(
    system.file("extdata", "exampledata.samples.txt", package = "methylumi")
  ))
  expect_s4_class(mr, "MethyLumiSet")
  expect_equal(unname(dim(mr)), c(1536L, 10L))
  expect_equal(sum(betas(mr)), 5691.890947, tolerance = TOL)
})

test_that("extractBarcodeAndPosition parses sentrix identifiers", {
  got <- extractBarcodeAndPosition("1632405013_R006_C012")
  expect_s3_class(got, "data.frame")
  expect_equal(names(got),
               c("sentrix", "row", "column", "rowNumber", "columnNumber"))
  expect_equal(as.character(got$sentrix), "1632405013")
  expect_equal(as.character(got$row), "R006")
  expect_equal(as.character(got$column), "C012")
  expect_equal(as.integer(got$rowNumber), 6L)
  expect_equal(as.integer(got$columnNumber), 12L)
})

test_that("extractBarcodeAndPosition is vectorised", {
  got <- extractBarcodeAndPosition(
    c("1632405013_R006_C012", "1632405013_R007_C001")
  )
  expect_equal(nrow(got), 2L)
  expect_equal(as.integer(got$rowNumber), c(6L, 7L))
  expect_equal(as.integer(got$columnNumber), c(12L, 1L))
})
