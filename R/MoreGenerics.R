# revamped for 1.9.10 and 2.0.0; checks to ensure that changes propagate
setReplaceMethod('sampleNames', signature(object="methylData",value='character'),#{{{
  function(object, value) {
    pd <- phenoData(object)
    sampleNames(pd) <- value
    ad <- assayData(object)
    sampleNames(ad) <- value
    prd <- protocolData(object)
    if (nrow(prd) == 0) prd <- pd[, integer(0)]
    else sampleNames(prd) <- value
    if( class(object) == 'MethyLumiSet' && !is.null(object@QC) ) { 
      # {{{
      qc = object@QC
      sampleNames(qc) = value
      object@QC = qc 
      # }}}
    } else if( class(object) == 'MethyLumiM' && !is.null(object@controlData) ) {
      # {{{
      cDat = object@controlData
      sampleNames(cDat) = value
      object@controlData = cDat
      # }}}
    }
    if('OOB' %in% slotNames(object)) { # {{{
      if(!is.null(object@OOB)) { 
        oob = object@OOB
        sampleNames(oob) <- value
        object@OOB = oob 
      } 
    } # }}}
    object@phenoData <- pd
    object@protocolData <- prd
    Biobase:::unsafeSetSlot(object, "assayData", ad)
  }) #}}}

# mostly QC and annotation functions
if(!isGeneric("diagnostics")) setGeneric("diagnostics",  # {{{
           function(x) standardGeneric('diagnostics')) # }}}
setMethod("diagnostics", signature(x="methylData"), function(x) { # {{{
  methylumi.diagnostics(x)
}) # }}}

# useful for comparing detection and bgcorrection results
setGeneric('sampleNAs', # {{{ should propagate through each data type
           function(object) standardGeneric('sampleNAs')) # }}}
setMethod("sampleNAs", signature(object="MethyLumiSet"), function(object){ # {{{
  colSums(is.na(betas(object)))
}) # }}}
setMethod("sampleNAs", signature(object="MethyLumiM"), function(object){ # {{{
  colSums(is.na(exprs(object)))
}) # }}}
setGeneric('probeNAs', # {{{ should propagate through each data type
           function(object) standardGeneric('probeNAs')) # }}}
setMethod("probeNAs", signature(object="MethyLumiSet"), function(object){ # {{{
  rowSums(is.na(betas(object)))
}) # }}}
setMethod("probeNAs", signature(object="MethyLumiM"), function(object){ # {{{
  rowSums(is.na(exprs(object)))
}) # }}}
setGeneric('plotNAs', # {{{ 
           function(object) standardGeneric('plotNAs')
           ) # }}}
setMethod("plotNAs", signature(object="methylData"), function(object){ # {{{
  pval <- pval.detect(object)
  sortorder <- order(sampleNames(object))
  sortedNames <- sampleNames(object)[sortorder]
  NAs <- data.frame(sample=sortedNames, index=1:length(sortedNames), 
                    dropouts=sampleNAs(object)[sortorder], 
                    slot=as.factor(sapply(sortedNames, function(x){
                      pop(strsplit(x, '_')[[1]])
                    })))
  NAs <- NAs[order(NAs$dropouts),]
  require('ggplot2')
  ggplot2::qplot(data=NAs, x=index, y=dropouts, size=dropouts, colour=slot,
                 geom=c('segment','point'), yend=0, xend=index, xlab='Sample #',
                 main=paste('Probe dropouts, colored by position, p >', pval))
}) # }}}
setGeneric('plotProbeNAs', # {{{ 
           function(object) standardGeneric('plotProbeNAs')
) # }}}
setMethod("plotProbeNAs",signature(object="methylData"),function(object){ # {{{
  require('ggplot2')
  pval <- pval.detect(object)
  x <- data.frame(drops=probeNAs(object)/dim(object)[2], 
                  mu=rowMeans(betas(object),na.rm=T))
  ggplot2::qplot(geom='jitter', x=mu, y=drops, ylab='failed probes',xlab='mean',
                 main=paste('Probe dropouts, colored by mean beta, p >', pval),
                 data=x, colour=mu)
}) # }}}

## FIXME: update these methods to work with older methylumi objects as well
setMethod('controlTypes', signature(object="MethyLumiSet"), #{{{
  function(object) controlTypes(object@QC)
) # }}}
setMethod('controlTypes', signature(object="MethyLumiM"), #{{{
  function(object) controlTypes(object@controlData)
) # }}}
setMethod('controlTypes', signature(object="MethyLumiQC"), #{{{
  function(object) levels( as.factor(fData(object)$Type) )
) # }}}

setMethod('betas', signature(object="MethyLumiM"), function(object) { # {{{
  (2**exprs(object))/(1+(2**exprs(object)))
}) # }}}
setMethod('pvals', signature(object="MethyLumiM"), function(object) { # {{{
  detection(object)
}) # }}}
setMethod('QCdata', signature(object="MethyLumiM"), # {{{
  function(object) controlData(object)) # }}}
setMethod('getHistory', signature(object="MethyLumiM"), # {{{
  function(object) object@history ) # }}}

if(!isGeneric('produceMethylationGEOSubmissionFile')) setGeneric('produceMethylationGEOSubmissionFile', # {{{
  function(object) standardGeneric('produceMethylationGEOSubmissionFile')) # }}}
setMethod('produceMethylationGEOSubmissionFile', signature(object="MethyLumiM"), # {{{
  function(object) {
    require(lumi)
    lumi:::produceMethylationGEOSubmissionFile(object)
  }) # }}}
setMethod('produceMethylationGEOSubmissionFile', signature(object="MethyLumiSet"), # {{{
  function(object) {
    require(lumi)
    lumi:::produceMethylationGEOSubmissionFile(as(object,'MethyLumiM'))
  }) # }}}
