\name{varFilter}

\alias{varFilter}
\alias{varFilter,MethyLumiSet-method}
\alias{varFilter,MethyLumiM-method}

\title{Variation-based Filtering of Features (CpG sites) in a
  MethyLumiSet or MethyLumiM object}

\description{The function \code{varFilter} removes features exhibiting
  little variation across samples. Such non-specific filtering can be
  advantageous for downstream data analysis. }

\usage{
varFilter(eset, var.func=IQR, var.cutoff=0.5, filterByQuantile=TRUE, ...)
}

\arguments{
  \item{eset}{An \code{MethyLumiSet} or \code{MethyLumiM} object.}
  \item{var.func}{The function used as the per-feature filtering
    statistics.}
  \item{var.cutoff}{A numeric value indicating the cutoff value for
    variation. If \code{filterByQuantile} is \code{TRUE}, features whose
    value of \code{var.func} is less than \code{var.cutoff}-quantile of
    all \code{var.func} value will be removed. It \code{FALSE}, features
    whose values are less than \code{var.cutoff} will be removed.}
  \item{filterByQuantile}{A logical indicating whether \code{var.cutoff}
    is to be interprested as a quantile of all \code{var.func} (the
    default), or as an absolute value.}
  \item{...}{Unused, but available for specializing methods.}
}

\value{
  The function \code{featureFilter} returns a list consisting of:
  \item{eset}{The filtered \code{MethyLumiSet} or \code{MethyLumiM} object.}
  \item{filter.log}{Shows many low-variant features are removed.}  
}

\details{
  This function is a counterpart of functions \code{nsFilter} and
  \code{varFilter} available from the \code{genefilter} package. See
  R. Bourgon et. al. (2010) and \code{\link[genefilter]{nsFilter}} for detail.

  It is proven that non-specific filtering, for which the criteria does
  not depend on sample class, can increase the number of discoverie.
  Inappropriate choice of test statistics, however, might have adverse
  effect.  \code{limma}'s moderated \eqn{t}-statistics, for example, is based on
  empirical Bayes approach which models the conjugate prior of
  gene-level variance with an inverse of \eqn{\chi^2} distribution scaled
  by observed global variance. As the variance-based filtering removes
  the set of genes with low variance, the scaled inverse \eqn{\chi^2}
  no longer provides a good fit to the data passing the filter,
  causing the \code{limma} algorithm to produce a posterior
  degree-of-freedom of infinity (Bourgon 2010). This leads to two
  consequences: (i) gene-level variance estimate will be ignore, and (ii)
  the \eqn{p}-value will be overly optimistic (Bourgon 2010). 
 
}

\references{
  R. Bourgon, R. Gentleman, W. Huber,
  \emph{Independent filtering increases power for detecting differentially
  expressed genes}, PNAS, vol. 107, no. 21, pp:9546-9551, 2010.}

\author{Chao-Jen Wong \email{cwon2@fhcrc.org}}

\seealso{\code{\link[genefilter]{nsFilter}}}

\examples{
  data(mldat)
  ## keep top 75 percent
  filt <- varFilter(mldat, var.cutoff=0.25)
  filt$filter.log
  dim(filt$eset)
}
