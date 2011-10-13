\name{featureFilter}

\alias{featureFilter}
\alias{featureFilter,MethyLumiSet-method}
\alias{featureFilter,MethyLumiM-method}

\title{Annotation-based Filtering of Features (CpG sites) in a
  MethyLumiSet or MethyLumiM object}

\description{Features with insufficient annotation
  carry little value for the subsequent data analysis. The function
  \code{featureFilter} provides options of filtering features (CpG
  sites) from a \code{MethyLumiSet} (or \code{MethyLumiM}) object based
  on available annotation data.
}

\usage{
featureFilter(eset, require.entrez=FALSE,
    require.GOBP=FALSE, require.GOCC=FALSE,
    require.GOMF=FALSE, exclude.ChrX=FALSE,
    require.closeToTSS=FALSE, range.DistToTSS=c(-500, 300),
    require.CpGisland=FALSE, ...)
}

\arguments{
  \item{eset}{A \code{MethyLumiSet} or \code{MethyLumiM} object.}
  \item{require.entrez}{If \code{TRUE}, filter out features without an Entrez
    Gene ID annotation.}
  
  \item{require.GOBP, require.GOCC, require.GOMF }{If \code{TRUE},
    filter out features whose target genes are not annotated to at least
    one GO term in BP, CC and MF ontology, respectively.}
  
  \item{exclude.ChrX}{If \code{TRUE}, filter out features in chromosome
    X to avoid gender effect.}
  
  \item{require.closeToTSS}{If \code{TRUE}, filter out features that
    are not close to transcription start site (TSS). Features without
    annotation of distance to TSS will also be removed. Can only used
    for GoldenGate platform.}

  \item{range.DistToTSS}{Ignored if \code{require.colseToTSS} is
    \code{FALSE}. A vector of numeric values of length 2, indicating the
    range of tolerable distance from transcription start site (TSS) in
    basepair (bp). If
    \code{require.clostToTSS} is \code{TRUE}, features whose distance
    to TSS falls outside this designated range will be removed. The
    default value is \eqn{c(-500, 300)}, where \eqn{-500} represents the
    distance to TSS from the left and 300 the distance from the right. }
  
  \item{require.CpGisland}{If \code{TRUE}, filter out features that are not in
    CpG islands.}
  
  \item{...}{Unused, but available for specializing methods.}
}

\value{
  The function \code{featureFilter} returns a list consisting of:
  \item{eset}{The filtered \code{MethyLumiSet} or \code{MethyLumiM} object.}
  \item{filter.log}{A list giving details of how many probe sets where
    removed for each annotation-based filtering step performed.}

}

\references{
  R. Bourgon, R. Gentleman, W. Huber,
  \emph{Independent filtering increases power for detecting differentially
  expressed genes}, PNAS, vol. 107, no. 21, pp:9546-9551.}

\author{Chao-Jen Wong \email{cwon2@fhcrc.org}}

\seealso{\code{\link[genefilter]{nsFilter}}}

\examples{
  data(mldat) ## GoldenGate Cancer Panel
  annotation(mldat) <- "GGHumanMethCancerPanelv1"
  filt <- featureFilter(mldat, require.entrez=TRUE,
            exclude.ChrX=TRUE,
            require.closeToTSS=TRUE, range.DistToTSS=c(-500, 300),
            require.CpGisland=TRUE)
  filt$filter.log
  dim(filt$eset)
}
