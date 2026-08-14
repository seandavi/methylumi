# methylumi 2.59.3

## Bug fixes

* `methylumIDAT()` reads IDAT files from per-slide subdirectories, which is the
  layout Illumina's own software produces (`<idatPath>/<Slide>/<Slide>_<Array>_Grn.idat`).
  Both the file-presence check and the code that actually opened the files
  assumed a flat directory, so a nested run reported every sample as missing.
  Flat layouts are unaffected. If the same basename appears under more than one
  subdirectory the read now errors and lists the candidates rather than picking
  one. `getBarcodes()` had the same blind spot and is fixed too. (#19)

* `diagnostics()` works on objects with no color-channel annotation, such as
  the bundled GoldenGate `mldat`. Three `if (annotation(x) == "...")` tests
  errored outright on a zero-length annotation, and past those, the plotting
  loop split probes on a `COLOR_CHANNEL` column that such objects do not have,
  so every panel would have been drawn from zero probes. The long-reported
  call to the non-existent `plot.density()` was real, on the same path, and is
  also fixed. (#23)

* `methylumi.bgcorr()` explains itself when asked for `noob` correction on an
  object that has no out-of-band intensities. Objects built from GEO series
  matrices or GenomeStudio output via `methylumiR()` carry only
  methylated/unmethylated/p-value data; OOB probes come only from reading raw
  IDATs. This used to fail deep inside method dispatch with "unable to find an
  inherited method for 'intensities.OOB'". (#24)

## Deprecated and defunct

* The gamma-family background corrections — `method = "goob"`, `"gamma"` and
  `"mode"` — have been removed. They called `gamma.mle()`, `gamma.mode()` and
  `gamma.integral()` from `rGammaGamma`, a GitHub-only package that is not
  available from CRAN or Bioconductor and was never declared as a dependency,
  so these methods have errored for years. Requesting one now produces an error
  saying so. `method = "noob"` (the default) is unaffected. (#18)

## Internal changes

* `qc.probe.plot()`, `plotNAs()` and `plotProbeNAs()` no longer use ggplot2's
  deprecated `qplot()`. Plot output is unchanged except in `plotNAs()`, where
  the segment layer now maps `linewidth` rather than the also-deprecated
  `size`; `scale_linewidth` uses a linear palette where `scale_size` used an
  area one, so line widths differ slightly. (#40)

* A pkgdown site is published to <https://seandavi.github.io/methylumi/> from a
  GitHub Actions workflow. (#43)

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
