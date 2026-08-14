## Golden values in these tests were captured from methylumi 2.59.1 -- the first
## version that installs on Bioconductor 3.24 -- and verified identical to the
## behaviour of the last release that installed anywhere. They exist to make
## behavioural drift visible during the modernization work in #30 and #34.
##
## Comparisons use relative tolerance rather than identity so that a BLAS or
## platform change does not turn the suite red for the wrong reason.

TOL <- 1e-6

## methylumIDAT() on the three bundled 27k IDATs takes several seconds, so
## import once and reuse across test files.
.cache <- new.env(parent = emptyenv())

idat_path <- function() system.file("extdata", package = "methylumi")

idat_barcodes <- function() {
  unique(sub("_(Grn|Red)\\.idat$", "",
             basename(list.files(idat_path(), pattern = "\\.idat$"))))
}

## Cached methylumIDAT() import of the bundled 27k IDATs.
example_idats <- function() {
  if (is.null(.cache$idats)) {
    .cache$idats <- suppressMessages(
      methylumIDAT(idat_barcodes(), idatPath = idat_path())
    )
  }
  .cache$idats
}

## Cached methylumiR() import of the bundled GoldenGate text files.
example_text <- function() {
  if (is.null(.cache$text)) {
    .cache$text <- suppressMessages(methylumiR(
      system.file("extdata", "exampledata.samples.txt", package = "methylumi"),
      qcfile = system.file("extdata", "exampledata.controls.txt", package = "methylumi")
    ))
  }
  .cache$text
}
