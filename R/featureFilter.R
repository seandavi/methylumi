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
