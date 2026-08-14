# methylumi 2.59.2

## Bug fixes

* The package installs again. `lattice` removed `parallel()` from its exports
  (it was renamed `parallelplot()` years ago, then deprecated, then made
  defunct). Because the import resolved at lazy-load time, this broke
  `R CMD INSTALL` for the whole package rather than just `parplot()`, and took
  several downstream packages with it.

* `QCdata()` is no longer silently lost when a `MethyLumiSet` is subset.
  `eSet`'s `[` method drops the `QC` slot, and the previous guard tried to
  subset it *after* that had already happened — `NULL[, j, drop = FALSE]`
  returns `NULL` without erroring, so the guard was a no-op. QC data is now
  captured before dispatch and reattached. Feature indices no longer apply to
  QC probes, which live in their own feature space; only sample selection does.
  As a consequence `combine()` no longer reports "Dropped control probes" when
  combining subsets.

* The example dataset `mldat` no longer shares its `assayData` between copies.
  It had been saved with `storageMode = "environment"` instead of
  `"lockedEnvironment"`, so `x <- mldat; methylated(x) <- ...` silently mutated
  the original. The assay values themselves are unchanged.

* `plotNegOob()` works again. It had two independent ggplot2 breakages:
  `scale_y_continuous(breaks = NA)` (now `NULL`) and `opts()` (removed from
  ggplot2 years ago, now `labs()`).

## Internal changes

* `.mclapply()` is now a thin wrapper around `BiocParallel::bplapply()` instead
  of `parallel::mclapply()` behind a `require()`. `mclapply()` forks, so it
  degraded to serial on Windows regardless. The signature is unchanged, and
  parallel results are verified identical to serial results. Configure it the
  standard Bioconductor way, with `BiocParallel::register()`.

* `require()` calls inside package code have been replaced with declared
  imports, or with `requireNamespace()` guards that produce an actionable error
  for genuinely optional packages (MASS, Biostrings, lumi). `require()` returns
  `FALSE` rather than erroring, so a missing package used to change behaviour
  silently.

* Bare `T`/`F` are now `TRUE`/`FALSE` throughout (105 occurrences). These are
  ordinary variables, not reserved words, so a user with a variable named `T`
  in scope could change the meaning of a default argument.

* `1:length(x)`, `1:nrow(x)` and `1:ncol(x)` are now `seq_along()`/`seq_len()`,
  which behave correctly when the length is zero.

* Base-package functions used by the package are now properly imported, rather
  than resolving by accident because the packages happened to be attached.

## Documentation

* Most man pages are generated from roxygen blocks now. The S4 class and
  generic documentation remains hand-written on purpose; see the note in
  `R/data.R`.

* `extractBarcodeAndPosition()` was documented as returning three columns with
  numeric row and column. It returns five, and `row`/`column` are strings such
  as `"R006"`; the numeric ones are `rowNumber`/`columnNumber`.

## Testing and infrastructure

* The package has a test suite for the first time: 36 tests covering import
  from IDAT and text, normalization, background correction, coercions,
  filtering, subsetting and every plotting entry point.

* CI runs `R CMD check` against Bioconductor devel on GitHub Actions,
  replacing a `.travis.yml` that had not run in years.

## Known issues

* `methylumi.bgcorr()` with `method = "gamma"` or `method = "mode"` does not
  work. It calls `gamma.mle()`, `gamma.mode()` and `gamma.integral()` from
  `rGammaGamma`, which is not declared as a dependency and is not available
  from CRAN or Bioconductor. The default `noob` method is unaffected.

* `qc.probe.plot()`, `plotNAs()` and `plotProbeNAs()` still use ggplot2's
  deprecated `qplot()`.
