require('Biobase')
setClass('MethyLumi', contains='eSet')
setClass('MethyLumiQC', contains='MethyLumi')
setClassUnion("QCDataOrNULL",c('NULL',"MethyLumiQC"))
setClassUnion("methylData",c('MethyLumi','ExpressionSet'))

## Generic methods
if (is.null(getGeneric("getHistory"))) setGeneric("getHistory", function(object) standardGeneric("getHistory"))
if (is.null(getGeneric("betas"))) setGeneric("betas", function(object) standardGeneric("betas"))
if (is.null(getGeneric("betas<-"))) setGeneric("betas<-", function(object, value) standardGeneric("betas<-"))
if (is.null(getGeneric("mvals"))) setGeneric("mvals", function(object) standardGeneric("mvals"))
if (is.null(getGeneric("mvals<-"))) setGeneric("mvals<-", function(object, value) standardGeneric("mvals<-"))
if (is.null(getGeneric("pvals"))) setGeneric("pvals", function(object) standardGeneric("pvals"))
if (is.null(getGeneric("pvals<-"))) setGeneric("pvals<-", function(object, value) standardGeneric("pvals<-"))
if (is.null(getGeneric("unmethylated"))) setGeneric("unmethylated", function(object) standardGeneric("unmethylated"))
if (is.null(getGeneric("unmethylated<-"))) setGeneric("unmethylated<-", function(object, value) standardGeneric("unmethylated<-"))
if (is.null(getGeneric("unmethylated.OOB"))) setGeneric("unmethylated.OOB", function(object) standardGeneric("unmethylated.OOB"))
if (is.null(getGeneric("unmethylated.OOB<-"))) setGeneric("unmethylated.OOB<-", function(object,value) standardGeneric("unmethylated.OOB<-"))
if (is.null(getGeneric("methylated"))) setGeneric("methylated", function(object) standardGeneric("methylated"))
if (is.null(getGeneric("methylated<-"))) setGeneric("methylated<-", function(object, value) standardGeneric("methylated<-"))
if (is.null(getGeneric("methylated.OOB"))) setGeneric("methylated.OOB", function(object) standardGeneric("methylated.OOB"))
if (is.null(getGeneric("methylated.OOB<-"))) setGeneric("methylated.OOB<-", function(object,value) standardGeneric("methylated.OOB<-"))
if (is.null(getGeneric("hist"))) setGeneric("hist", function(x,...) standardGeneric("hist"))
if (is.null(getGeneric("corplot"))) setGeneric("corplot", function(x,...) standardGeneric("corplot"))
if (is.null(getGeneric("plotDensity"))) setGeneric("plotDensity", function(x,...) standardGeneric("plotDensity"))
if (is.null(getGeneric("QCdata"))) setGeneric("QCdata", function(object) standardGeneric("QCdata"))
if (is.null(getGeneric("QCdata<-"))) setGeneric("QCdata<-", function(object, value) standardGeneric("QCdata<-"))

# generic methods specific for MethyLumiM class
if (is.null(getGeneric("methylated.N"))) setGeneric("methylated.N", function(object) standardGeneric("methylated.N"))
if (is.null(getGeneric("methylated.N<-"))) setGeneric("methylated.N<-", function(object, value) standardGeneric("methylated.N<-"))
if (is.null(getGeneric("unmethylated.N"))) setGeneric("unmethylated.N", function(object) standardGeneric("unmethylated.N"))
if (is.null(getGeneric("unmethylated.N<-"))) setGeneric("unmethylated.N<-", function(object, value) standardGeneric("unmethylated.N<-"))
if (is.null(getGeneric("controlData"))) setGeneric("controlData", function(object) standardGeneric("controlData"))
if (is.null(getGeneric("controlData<-"))) setGeneric("controlData<-", function(object, value) standardGeneric("controlData<-"))
if (is.null(getGeneric("detection"))) setGeneric("detection", function(object) standardGeneric("detection"))
if (is.null(getGeneric("detection<-"))) setGeneric("detection<-", function(object, value) standardGeneric("detection<-"))


if (is.null(getGeneric("featureFilter"))) # {{{
    setGeneric("featureFilter", 
               function(eset,
                        require.entrez=FALSE,
                        require.GOBP=FALSE,
                        require.GOCC=FALSE,
                        require.GOMF=FALSE,
                        exclude.ChrX=FALSE,
                        require.closeToTSS=FALSE,
                        range.DistToTSS=c(-500, 300),
                        require.CpGisland=FALSE, ...)
               standardGeneric("featureFilter")) # }}}

if (is.null(getGeneric("varFilter"))) # {{{
    setGeneric("varFilter", 
               function(eset,
                        var.func=IQR, var.cutoff=0.5,
                        filterByQuantile=TRUE, ...)
               standardGeneric("varFilter")) # }}}

setGeneric('total.intensity', # {{{
             function(object) standardGeneric('total.intensity')) # }}}
setMethod('total.intensity', signature(object="MethyLumiSet"), # {{{
          function(object) return(methylated(object)+unmethylated(object))) #}}}
setMethod('total.intensity', signature(object="MethyLumiM"), # {{{
          function(object) return(methylated(object)+unmethylated(object))) #}}}
setMethod('methylated', signature(object="MethyLumiQC"), # {{{
          function(object) return(assayDataElement(object, 'methylated'))) # }}}
setMethod('unmethylated', signature(object="MethyLumiQC"), # {{{
          function(object) return(assayDataElement(object, 'unmethylated')))#}}}

setGeneric('unmethylated.N', # {{{
              function(object) standardGeneric('unmethylated.N')) # }}}
setMethod("unmethylated.N", signature(object="methylData"),function(object){#{{{
            return(assayDataElement(object,"unmethylated.N"))
}) # }}}
setGeneric('unmethylated.N<-', # {{{
              function(object,value) standardGeneric('unmethylated.N<-')) # }}}
setReplaceMethod("unmethylated.N", signature(object="methylData",value="matrix"), function(object, value) { # {{{
            assayDataElementReplace(object,"unmethylated.N",value)
          }) # }}}
setGeneric('methylated.N', # {{{
              function(object) standardGeneric('methylated.N')) # }}}
setMethod("methylated.N", signature(object="methylData"), function(object){ #{{{
            return(assayDataElement(object,"methylated.N"))
}) # }}}
setGeneric('methylated.N<-', # {{{
              function(object, value) standardGeneric('methylated.N<-')) # }}}
setReplaceMethod("methylated.N", signature(object="methylData",value="matrix"), function(object,value) { # {{{
            assayDataElementReplace(object,"methylated.N",value)
          }) # }}}

## FIXME: use COLOR_CHANNEL for this, and add 'both'/'d2', or else retire it
setGeneric('getProbesByChannel', # {{{
  function(object, ...) standardGeneric('getProbesByChannel')) # }}}
setMethod('getProbesByChannel', signature(object="methylData"), # {{{
  function(object, channel=NULL, channels=c('Cy3','Cy5')) {
    if(is.null(channel)) { 
      perchannel <- lapply(channels, function(channel)
        intensitiesByChannel(object, channel=channel, allele=allele))
      names(perchannel) <- channels
      return(perchannel)
    } 
    if( tolower(channel) %in% c('cy5','red') ) {
      if(annotation(object)=='IlluminaHumanMethylation27k') return(cy5(object))
      else stop('450k Design I probes are red or green; Design II are both')
    }
    if( tolower(channel) %in% c('cy3','green') ) {
      if(annotation(object)=='IlluminaHumanMethylation27k') return(cy3(object))
      else stop('450k Design I probes are red or green; Design II are both')
    }
  }) # }}}

setGeneric('intensitiesByChannel', # {{{
  function(object, ...) standardGeneric('intensitiesByChannel')) # }}}
setMethod('intensitiesByChannel', signature(object="MethyLumiSet"), # {{{
  function(object, channel=NULL, allele=NULL, 
           channels=c('Cy3','Cy5'), alleles=c('methylated','unmethylated')) {
  if(is.null(channel)) { # {{{
    perchannel <- lapply(channels, function(channel)
      intensitiesByChannel(object, channel=channel, allele=allele))
    names(perchannel) <- channels
    return(perchannel)
  } # }}} 
  if(is.null(allele)) { # {{{
    perallele <- lapply(alleles, function(allele) 
      intensitiesByChannel(object, channel=channel, allele=allele))
    names(perallele) <- alleles
    return(perallele)
  } # }}}
    return(assayDataElement(object,allele)[getProbesByChannel(object,channel),])
  }) # }}}
setMethod('intensitiesByChannel', signature(object="MethyLumiM"), # {{{
  function(object, channel=NULL, allele=NULL, 
           channels=c('Cy3','Cy5'), alleles=c('methylated','unmethylated')) {
  if(is.null(channel)) { # {{{
    perchannel <- lapply(channels, function(channel)
      intensitiesByChannel(object, channel=channel, allele=allele))
    names(perchannel) <- channels
    return(perchannel)
  } # }}} 
  if(is.null(allele)) { # {{{
    perallele <- lapply(alleles, function(allele) 
      intensitiesByChannel(object, channel=channel, allele=allele))
    names(perallele) <- alleles
    return(perallele)
  } # }}}
    return(assayDataElement(object,allele)[getProbesByChannel(object,channel),])
  }) # }}}
setMethod('intensitiesByChannel', signature(object="MethyLumiQC"), # {{{
  function(object, channel=NULL, channels=c('Cy3','Cy5')) {
  if(is.null(channel)) { # {{{
    perchannel <- lapply(channels, function(channel)
      intensitiesByChannel(object, channel=channel))
    names(perchannel) <- channels
    return(perchannel)
  } # }}}
  if( channel=='Cy3' ) return(Cy3(object))
  if( channel=='Cy5' ) return(Cy5(object))
  else stop(paste("Don't know where to find channel", channel, "!"))
}) # }}}

setGeneric('Cy3', # {{{
             function(object) standardGeneric('Cy3')) # }}}
setMethod('Cy3', signature(object="MethyLumiQC"), # {{{
          function(object) return(assayDataElement(object,'methylated'))) #}}}
setGeneric("Cy3<-", # {{{
              function(object, value) standardGeneric('Cy3<-')) # }}}
setReplaceMethod('Cy3', signature(object="MethyLumiQC", value="matrix"), # {{{
  function(object, value) {
    assayDataElementReplace(object, 'methylated', value)
  }) # }}}

setGeneric('Cy3.SD', # {{{
           function(object) standardGeneric('Cy3.SD')) # }}}
setMethod('Cy3.SD', signature(object="MethyLumiSet"), # {{{
          function(object) Cy3.SD(object@QC) ) #}}}
setMethod('Cy3.SD', signature(object="MethyLumiQC"), # {{{
          function(object) {
            if( is.element('methylated.SD', assayDataElementNames(object)) ) {
              return(assayDataElement(object,'methylated.SD'))
            } else {
              return(NULL)
            }
          }) #}}}
setGeneric('Cy3.N', # {{{
           function(object) standardGeneric('Cy3.N')) # }}}
setMethod('Cy3.N', signature(object="MethyLumiSet"), # {{{
          function(object) Cy3.N(object@QC) ) #}}}
setMethod('Cy3.N', signature(object="MethyLumiQC"), # {{{
          function(object) {
            if( is.element('NBeads', assayDataElementNames(object)) ) {
              return(assayDataElement(object,'NBeads'))
            } else {
              return(NULL)
            }
          }) #}}}

setGeneric('Cy5', # {{{
             function(object) standardGeneric('Cy5')) # }}}
setMethod('Cy5', signature(object="MethyLumiQC"), # {{{
          function(object) return(assayDataElement(object,'unmethylated'))) #}}}
setGeneric("Cy5<-", # {{{
              function(object, value) standardGeneric('Cy5<-')) # }}}
setReplaceMethod('Cy5', signature(object="MethyLumiQC", value="matrix"), # {{{
  function(object, value) {
    assayDataElementReplace(object, 'unmethylated', value)
  }) # }}}

setGeneric('Cy5.SD', # {{{
          function(object) standardGeneric('Cy5.SD')) # }}}
setMethod('Cy5.SD', signature(object="MethyLumiSet"), # {{{
          function(object) Cy5.SD(object@QC)) #}}}
setMethod('Cy5.SD', signature(object="MethyLumiQC"), # {{{
          function(object) {
            if( is.element('unmethylated.SD', assayDataElementNames(object)) ) {
              return(assayDataElement(object,'unmethylated.SD'))
            } else {
              return(NULL)
            }
          }) #}}}
setGeneric('Cy5.N', # {{{
           function(object) standardGeneric('Cy5.N')) # }}}
setMethod('Cy5.N', signature(object="MethyLumiSet"), # {{{
          function(object) Cy5.N(object@QC)) #}}}
setMethod('Cy5.N', signature(object="MethyLumiQC"), # {{{
          function(object) {
            if( is.element('NBeads', assayDataElementNames(object)) ) {
              return(assayDataElement(object,'NBeads'))
            } else {
              return(NULL)
            }
          }) #}}}

setGeneric('negctls', # {{{
           function(object, channel) standardGeneric('negctls')) # }}}
setMethod('negctls', signature(object="MethyLumiSet",channel='character'), # {{{
          function(object,channel) return(negctls(controlData(object),channel)))# }}}
setMethod('negctls', signature(object="MethyLumiSet",channel='missing'), # {{{
          function(object) {
            channels = list(Cy3='Cy3',Cy5='Cy5')
            lapply(channels, function(channel) negctls(controlData(object), channel))
          }) # }}}
setMethod('negctls', signature(object="MethyLumiM",channel='character'), # {{{
          function(object,channel) return(negctls(controlData(object),channel)))# }}}
setMethod('negctls', signature(object="MethyLumiM",channel='missing'), # {{{
          function(object) {
            channels = list(Cy3='Cy3',Cy5='Cy5')
            lapply(channels, function(channel) negctls(controlData(object), channel))
          }) # }}}
setMethod('negctls', signature(object="MethyLumiQC",channel='character'), # {{{
          function(object, channel) {
            negs <- grep('negative', tolower(featureNames(object)))
            if(channel %in% c('Cy3','Grn','G')) {
              return(methylated(object)[negs,])
            } else if(channel %in% c('Cy5','Red','R')) {
              return(unmethylated(object)[negs,])
            }
          }) #}}}
setMethod('negctls', signature(object="MethyLumiQC",channel='missing'), # {{{
          function(object) {
            channels = list(Cy3='Cy3',Cy5='Cy5')
            lapply(channels, function(channel) negctls(object, channel))
          }) #}}}

setGeneric('negctls.SD', # {{{
  function(object, channel) standardGeneric('negctls.SD')) # }}}
setMethod('negctls.SD',signature(object="MethyLumiSet",channel='character'),#{{{
  function(object,channel) return(negctls.SD(QCdata(object),channel)))#}}}
setMethod('negctls.SD',signature(object="MethyLumiM",channel='character'),#{{{
  function(object,channel) return(negctls.SD(controlData(object),channel)))#}}}
setMethod('negctls.SD', signature(object="MethyLumiQC",channel='character'),#{{{
          function(object, channel) {
            negs <- grep('negative', tolower(featureNames(object)))
            if(channel %in% c('Cy3','Grn','G')) {
              return(Cy3.SD(object)[negs,])
            } else if(channel %in% c('Cy5','Red','R')) {
              return(Cy5.SD(object)[negs,])
            }
          }) #}}}
setMethod('negctls.SD', signature(object="MethyLumiQC",channel='missing'),#{{{
          function(object) {
            channels = list(Cy3='Cy3',Cy5='Cy5')
            lapply(channels, function(channel) negctls.SD(object, channel))
          }) #}}}

setGeneric('negctls.stderr', # {{{
  function(object, channel) standardGeneric('negctls.stderr')) # }}}
setMethod('negctls.stderr', signature(object="MethyLumiSet",channel='character'),#{{{
          function(object, channel) negctls.stderr(QCdata(object), channel))#}}}
setMethod('negctls.stderr', signature(object="MethyLumiQC",channel='character'),#{{{
          function(object, channel) {
            negs <- grep('negative', tolower(featureNames(object)))
            if(channel %in% c('Cy3','Grn','G')) {
              return( (Cy3.SD(object)/sqrt(Cy3.N(object)))[negs,] )
            } else if(channel %in% c('Cy5','Red','R')) {
              return( (Cy5.SD(object)/sqrt(Cy5.N(object)))[negs,] )
            }
          }) #}}}
setMethod('negctls.stderr', signature(object="MethyLumiSet", channel='missing'),  # {{{
          function(object) {
            channels = list(Cy3='Cy3',Cy5='Cy5')
            lapply(channels, function(channel) negctls.stderr(object, channel))
          }) #}}}
setMethod('negctls.stderr', signature(object="MethyLumiQC", channel='missing'),  # {{{
          function(object) {
            channels = list(Cy3='Cy3',Cy5='Cy5')
            lapply(channels, function(channel) negctls.stderr(object, channel))
          }) #}}}

setGeneric('negnorm', # {{{
           function(object, channel) standardGeneric('negnorm')) # }}}
setMethod('negnorm', signature(object="MethyLumiSet",channel='character'), # {{{
          function(object,channel) return(negnorm(QCdata(object),channel)))# }}}
setMethod('negnorm', signature(object="MethyLumiM",channel='character'), # {{{
  function(object,channel) return(negnorm(controlData(object),channel)))# }}}
setMethod('negnorm', signature(object="MethyLumiQC",channel='character'), # {{{
          function(object, channel) {
            if(channel %in% c('Cy3','Grn','G')) {
              negs = grep('negative',tolower(featureNames(object))) 
              norms = grep('norm.*?(a|t|red)',tolower(featureNames(object))) 
              return(Cy3(object)[c(negs,norms), ])
            } else if(channel %in% c('Cy5','Red','R')) {
              negs = grep('negative',tolower(featureNames(object))) 
              norms = grep('norm.*?(g|c|green)',tolower(featureNames(object))) 
              return(Cy5(object)[c(negs,norms),])
            }
          }) #}}}
setMethod('negnorm', signature(object="MethyLumiSet",channel='missing'), # {{{
          function(object,channel) {
            channels = list(Cy3='Cy3',Cy5='Cy5')
            lapply(channels, function(channel) negnorm(object@QC, channel))
          }) #}}}
setMethod('negnorm', signature(object="MethyLumiM",channel='missing'), # {{{
  function(object,channel) {
    channels = list(Cy3='Cy3',Cy5='Cy5')
    lapply(channels, function(channel) negnorm(controlData(object), channel))
  }) #}}}
setMethod('negnorm', signature(object="MethyLumiQC",channel='missing'), # {{{
          function(object, channel) {
            channels = list(Cy3='Cy3',Cy5='Cy5')
            lapply(channels, function(channel) negnorm(object, channel))
          }) #}}}

setGeneric('normctls', # {{{
  function(object, ...) standardGeneric('normctls')) # }}}
setMethod('normctls', signature(object="MethyLumiSet"), # {{{
          function(object, ...) return(normctls(QCdata(object), ...)))#}}}
setMethod('normctls', signature(object="MethyLumiM"), # {{{
          function(object, ...) return(normctls(controlData(object), ...)))#}}}
setMethod('normctls', signature(object="MethyLumiQC"), # {{{
          function(object, channel=NULL) {

            if(is.null(channel)) {
              return(lapply(list(Cy3='Cy3',Cy5='Cy5'), function(ch){
                normctls(object, ch)
              }))
            }

            ## This will depend on whether the chip is a 450k or 27k chip...
            if(annotation(object) == 'IlluminaHumanMethylation27k'){
              if(channel %in% c('Cy3','Grn')) searchfor <- 'Norm.G'
              if(channel %in% c('Cy5','Red')) searchfor <- 'Norm.R'
            } 
            if(annotation(object) == 'IlluminaHumanMethylation450k'){
              if(channel %in% c('Cy3','Grn')) searchfor <- 'Norm_(C|G)' # ?
              if(channel %in% c('Cy5','Red')) searchfor <- 'Norm_(A|T)' # ?
              # for sanity checking and exploration
              if(channel %in% c('A','T','G','C')) {
                searchfor <- paste('Norm',channel,sep='_')
              }
            } 

            probes <- grep(tolower(searchfor), tolower(featureNames(object)))
            if(channel %in% c('Cy3','Grn','G','C')) {
              return(methylated(object)[probes,])
            } else if(channel %in% c('Cy5','Red','A','T')) {
              return(unmethylated(object)[probes,])
            }
            
          }) #}}}
