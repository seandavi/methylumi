#' Variation-based filtering of features (CpG sites) in a MethyLumiSet or
#' MethyLumiM object
#'
#' The function `varFilter` removes features exhibiting little variation across
#' samples. Such non-specific filtering can be advantageous for downstream data
#' analysis.
#'
#' This function is a counterpart of the functions `nsFilter` and `varFilter`
#' available from the `genefilter` package. See R. Bourgon et al. (2010) and
#' [genefilter::nsFilter()] for detail.
#'
#' It is proven that non-specific filtering, for which the criterion does not
#' depend on sample class, can increase the number of discoveries.
#' Inappropriate choice of test statistics, however, might have an adverse
#' effect. `limma`'s moderated \eqn{t}-statistic, for example, is based on an
#' empirical Bayes approach which models the conjugate prior of gene-level
#' variance with an inverse of \eqn{\chi^2} distribution scaled by observed
#' global variance. As the variance-based filtering removes the set of genes
#' with low variance, the scaled inverse \eqn{\chi^2} no longer provides a good
#' fit to the data passing the filter, causing the `limma` algorithm to produce
#' a posterior degree-of-freedom of infinity (Bourgon 2010). This leads to two
#' consequences: (i) the gene-level variance estimate will be ignored, and (ii)
#' the \eqn{p}-value will be overly optimistic (Bourgon 2010).
#'
#' @param eset A `MethyLumiSet` or `MethyLumiM` object.
#' @param var.func The function used as the per-feature filtering statistic.
#' @param var.cutoff A numeric value indicating the cutoff value for variation.
#'   If `filterByQuantile` is `TRUE`, features whose value of `var.func` is less
#'   than the `var.cutoff`-quantile of all `var.func` values will be removed. If
#'   `FALSE`, features whose values are less than `var.cutoff` will be removed.
#' @param filterByQuantile A logical indicating whether `var.cutoff` is to be
#'   interpreted as a quantile of all `var.func` values (the default), or as an
#'   absolute value.
#' @param ... Unused, but available for specializing methods.
#' @return A list consisting of:
#'   \describe{
#'     \item{eset}{The filtered `MethyLumiSet` or `MethyLumiM` object.}
#'     \item{filter.log}{Shows how many low-variance features were removed.}
#'   }
#' @references R. Bourgon, R. Gentleman, W. Huber, \emph{Independent filtering
#'   increases power for detecting differentially expressed genes}, PNAS,
#'   vol. 107, no. 21, pp. 9546-9551, 2010.
#' @author Chao-Jen Wong \email{cwon2@fhcrc.org}
#' @seealso [genefilter::nsFilter()]
#' @examples
#'   data(mldat)
#'   ## keep top 75 percent
#'   filt <- varFilter(mldat, var.cutoff=0.25)
#'   filt$filter.log
#'   dim(filt$eset)
#' @name varFilter
#' @aliases varFilter varFilter,MethyLumiSet-method varFilter,MethyLumiM-method
#' @usage varFilter(eset, var.func=IQR, var.cutoff=0.5, filterByQuantile=TRUE, ...)
NULL

## copy from nsFilter.R from the genefilter package
rowIQRs <- function(eSet) {
  numSamp <- ncol(eSet)
  lowQ <- rowQ(eSet, floor(0.25 * numSamp))
  upQ <- rowQ(eSet, ceiling(0.75 * numSamp))
  upQ - lowQ
}

## warpping the varFilter function from the genefilter package
.varFilter <- function(eset,
                   var.func=IQR, var.cutoff=0.5,
                   filterByQuantile=TRUE, ...)
{
    if (!is.function(var.func))
        stop("'var.func' must be a function")
            
    filter.log <- new.env(parent=emptyenv())
            
    if (deparse(substitute(var.func)) == "IQR") {
        esetIqr <- rowIQRs(exprs(eset))
    } else {
        esetIqr <- apply(exprs(eset), 1, var.func)
    }
 
    if (filterByQuantile) {
        if ( 0 < var.cutoff && var.cutoff < 1 ) {
            var.cutoff = quantile(esetIqr, var.cutoff, na.rm=TRUE)
        } else stop("Cutoff Quantile has to be between 0 and 1.")
    }
    selected <- esetIqr > var.cutoff
    eset <- eset[selected, ]
    logvar <- "numLowVar"
    assign(logvar, sum(!selected), filter.log)
    list(eset=eset, filter.log=as.list(filter.log))
}


          
setMethod("varFilter", "MethyLumiSet",
          function(eset,
                   var.func=IQR, var.cutoff=0.5,
                   filterByQuantile=TRUE, ...)
          {
              .varFilter(eset, var.func, var.cutoff, filterByQuantile)
          }
)

setMethod("varFilter", "MethyLumiM",
          function(eset,
                   var.func=IQR, var.cutoff=0.5,
                   filterByQuantile=TRUE, ...)
          {
              .varFilter(eset, var.func, var.cutoff, filterByQuantile)
          }
)
