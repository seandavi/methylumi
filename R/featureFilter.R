#' Annotation-based filtering of features (CpG sites) in a MethyLumiSet or
#' MethyLumiM object
#'
#' Features with insufficient annotation carry little value for subsequent data
#' analysis. The function `featureFilter` provides options for filtering
#' features (CpG sites) from a `MethyLumiSet` (or `MethyLumiM`) object based on
#' available annotation data.
#'
#' @param eset A `MethyLumiSet` or `MethyLumiM` object.
#' @param require.entrez If `TRUE`, filter out features without an Entrez Gene
#'   ID annotation.
#' @param require.GOBP,require.GOCC,require.GOMF If `TRUE`, filter out features
#'   whose target genes are not annotated to at least one GO term in the BP, CC
#'   and MF ontology, respectively.
#' @param exclude.ChrX If `TRUE`, filter out features on chromosome X to avoid
#'   gender effects.
#' @param require.closeToTSS If `TRUE`, filter out features that are not close
#'   to a transcription start site (TSS). Features without an annotation of
#'   distance to TSS will also be removed. Can only be used for the GoldenGate
#'   platform.
#' @param range.DistToTSS Ignored if `require.closeToTSS` is `FALSE`. A numeric
#'   vector of length 2 indicating the range of tolerable distance from the
#'   transcription start site (TSS) in basepairs. Features whose distance to TSS
#'   falls outside this range will be removed. The default is
#'   \eqn{c(-500, 300)}, where \eqn{-500} is the distance to TSS from the left
#'   and 300 the distance from the right.
#' @param require.CpGisland If `TRUE`, filter out features that are not in CpG
#'   islands.
#' @param ... Unused, but available for specializing methods.
#' @return A list consisting of:
#'   \describe{
#'     \item{eset}{The filtered `MethyLumiSet` or `MethyLumiM` object.}
#'     \item{filter.log}{A list giving details of how many probe sets were
#'       removed for each annotation-based filtering step performed.}
#'   }
#' @references R. Bourgon, R. Gentleman, W. Huber, \emph{Independent filtering
#'   increases power for detecting differentially expressed genes}, PNAS,
#'   vol. 107, no. 21, pp. 9546-9551.
#' @author Chao-Jen Wong \email{cwon2@fhcrc.org}
#' @seealso [genefilter::nsFilter()]
#' @name featureFilter
#' @aliases featureFilter featureFilter,MethyLumiSet-method featureFilter,MethyLumiM-method
#' @usage featureFilter(eset, require.entrez=FALSE,
#'     require.GOBP=FALSE, require.GOCC=FALSE,
#'     require.GOMF=FALSE, exclude.ChrX=FALSE,
#'     require.closeToTSS=FALSE, range.DistToTSS=c(-500, 300),
#'     require.CpGisland=FALSE, ...)
NULL

.featureFilter <- function(eset,
           require.entrez=FALSE,
           require.GOBP=FALSE,
           require.GOCC=FALSE,
           require.GOMF=FALSE,
           exclude.ChrX=FALSE,
           require.closeToTSS=FALSE,
           range.DistToTSS=c(-500, 300),
           require.CpGisland=FALSE, ...)
{
     annChip <- annotation(eset)
     if (annChip == 'IlluminaHumanMethylation450k') {
       warning("HumanMethylation450k probes annotate to multiple accessions(!)")
     } else if (nchar(annChip) == 0) {
       stop("'eset' must have a valid annotation slot")
     }

     nfeat <- function(eset) length(featureNames(eset))
     filter.log <- new.env(parent=emptyenv())
     
     if (require.entrez) {
         map <- genefilter:::.findCentralMap(annChip)
         IDs <-  mget(featureNames(eset),
                      envir=getAnnMap(map, annChip),
                      ifnotfound=NA)
         haveID <-  names(IDs)[sapply(IDs, function(x) !is.na(x))]
         logvar <- paste("numRemoved", map, sep=".")
         assign(logvar, nfeat(eset) - length(haveID), envir=filter.log)
         eset <- eset[haveID, ]
     }

     if (require.closeToTSS) {
         if (is.null(range.DistToTSS)) range.DistToTSS <- c(-500, 300)
         if (!is.numeric(range.DistToTSS))
             stop("'range.DistToTSS' must be a vector of numeric values.")
         if (length(range.DistToTSS) != 2)
             stop("The length of 'range.DistToTSS' must be 2.")
   
               
         if (annChip == 'IlluminaHumanMethylation450k') {
           stop('DISTTOTSS can have multiple values for 450k probes')
         }
         map <- "DISTTOTSS"
         distotss <- mget(featureNames(eset),
                     envir=annotate::getAnnMap(map, annChip),
                     ifnotfound=NA)
         closetotss <- sapply(distotss, function(x) {
                              if (length(x) == 1 && is.na(x))
                                 FALSE ## no distance annotation available
                              else
                                 ifelse(x > range.DistToTSS[1] &
                                        x <= range.DistToTSS[2], TRUE, FALSE)
                            })
         logvar <- paste("numNotClose", map, sep=".")
         assign(logvar, sum(!closetotss), envir=filter.log)
         eset <- eset[closetotss, ]
     }

     if (require.CpGisland) {

         map <- "ISCPGISLAND"
         CpGisland <- mget(featureNames(eset),
                     envir=annotate::getAnnMap(map, annChip),
                     ifnotfound=NA)
         isCpGisland <- sapply(CpGisland, function(x) {
                               if (is.na(x))
                                   FALSE
                               else as.logical(x)
                        })
         logvar <- paste("numNot", map, sep=".")
         assign(logvar, sum(!isCpGisland), envir=filter.log)
         eset <- eset[isCpGisland, ]
     }

     if (exclude.ChrX) {
         chr <- mget(featureNames(eset),
                     envir=annotate::getAnnMap("CHR", annChip),
                     ifnotfound=NA)
         notX <- sapply(chr, function(x) x[1]!="X")
         notX[is.na(notX)] <- TRUE
         logvar <- "numChromX"
         assign(logvar, sum(!notX), envir=filter.log)
         eset <- eset[notX, ]
     }

     ## same as what's in genefiler::nsFilter
     filterGO <- function(eset, ontology) {
                  haveGo <- sapply(mget(featureNames(eset),
                                     getAnnMap("GO", annChip), ifnotfound=NA),
                                   function(x) {
                                       if (length(x) == 1 && is.na(x))
                                         FALSE
                                       else {
                                           onts <- subListExtract(x,
                                                   "Ontology", simplify=TRUE)
                                           ontology %in% onts
                                       }
                                   })
                  logvar <- paste("numNoGO", ontology, sep=".")
                  assign(logvar, sum(!haveGo), envir=filter.log)
                  eset[haveGo, ]
              }

              if (require.GOBP) {
                  eset <- filterGO(eset, "BP")
              }

              if (require.GOCC) {
                  eset <- filterGO(eset, "CC")
              }

              if (require.GOMF) {
                  eset <- filterGO(eset, "MF")
              }
     
    return(list(eset=eset, filter.log=as.list(filter.log)))
}


setMethod("featureFilter", signature(eset="MethyLumiSet"),
  function(eset,
           require.entrez=FALSE,
           require.GOBP=FALSE,
           require.GOCC=FALSE,
           require.GOMF=FALSE,
           exclude.ChrX=FALSE,
           require.closeToTSS=FALSE,
           range.DistToTSS=c(-500, 300),
           require.CpGisland=FALSE, ...)
    {
    .featureFilter(eset,
                   require.entrez=require.entrez,
                   require.GOBP=require.GOBP,
                   require.GOCC=require.GOCC,
                   require.GOMF=require.GOMF,
                   exclude.ChrX=exclude.ChrX,
                   require.closeToTSS=require.closeToTSS,
                   range.DistToTSS=range.DistToTSS,
                   require.CpGisland=require.CpGisland)
    }
)

setMethod("featureFilter", signature(eset="MethyLumiM"),
  function(eset,
           require.entrez=FALSE,
           require.GOBP=FALSE,
           require.GOCC=FALSE,
           require.GOMF=FALSE,
           exclude.ChrX=FALSE,
           require.closeToTSS=FALSE,
           range.DistToTSS=c(-500, 300),
           require.CpGisland=FALSE, ...)
    {
    .featureFilter(eset,
                   require.entrez=require.entrez,
                   require.GOBP=require.GOBP,
                   require.GOCC=require.GOCC,
                   require.GOMF=require.GOMF,
                   exclude.ChrX=exclude.ChrX,
                   require.closeToTSS=require.closeToTSS,
                   range.DistToTSS=range.DistToTSS,
                   require.CpGisland=require.CpGisland)
    }
)
