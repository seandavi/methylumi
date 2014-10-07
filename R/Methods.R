setMethod('total.intensity', signature(object="MethyLumiSet"), # {{{
          function(object) return(methylated(object)+unmethylated(object))) #}}}
setMethod('total.intensity', signature(object="MethyLumiM"), # {{{
          function(object) return(methylated(object)+unmethylated(object))) #}}}
setMethod('methylated', signature(object="MethyLumiQC"), # {{{
          function(object) return(assayDataElement(object, 'methylated'))) # }}}
setMethod('unmethylated', signature(object="MethyLumiQC"), # {{{
          function(object) return(assayDataElement(object, 'unmethylated')))#}}}

setMethod("unmethylated.N", signature(object="methylData"),function(object){#{{{
  return(assayDataElement(object,"unmethylated.N"))
}) # }}}

setReplaceMethod("unmethylated.N", signature(object="methylData",value="matrix"), function(object, value) { # {{{
  assayDataElementReplace(object,"unmethylated.N",value)
}) # }}}

setMethod("methylated.N", signature(object="methylData"), function(object){ #{{{
  return(assayDataElement(object,"methylated.N"))
}) # }}}

setReplaceMethod("methylated.N", signature(object="methylData",value="matrix"), function(object,value) { # {{{
  assayDataElementReplace(object,"methylated.N",value)
}) # }}}

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
            return(assayDataElement(object,allele)[getProbesByChannel(object,channel),,drop=FALSE])
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
            return(assayDataElement(object,allele)[getProbesByChannel(object,channel),,drop=FALSE])
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

setMethod('Cy3', signature(object="MethyLumiQC"), # {{{
          function(object) return(assayDataElement(object,'methylated'))) #}}}

setReplaceMethod('Cy3', signature(object="MethyLumiQC", value="matrix"), # {{{
                 function(object, value) {
                   assayDataElementReplace(object, 'methylated', value)
                 }) # }}}

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

setMethod('Cy5', signature(object="MethyLumiQC"), # {{{
          function(object) return(assayDataElement(object,'unmethylated'))) #}}}
setReplaceMethod('Cy5', signature(object="MethyLumiQC", value="matrix"), # {{{
                 function(object, value) {
                   assayDataElementReplace(object, 'unmethylated', value)
                 }) # }}}
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

setMethod('negctls', signature(object="MethyLumiSet",channel='character'), # {{{
          function(object,channel) return(negctls(QCdata(object),channel)))# }}}
setMethod('negctls', signature(object="MethyLumiSet",channel='missing'), # {{{
          function(object) {
            channels = list(Cy3='Cy3',Cy5='Cy5')
            lapply(channels, function(channel) negctls(QCdata(object), channel))
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
            negs0 <- grep('negative', tolower(featureNames(object)))
            negs1 <- grep('-99', fData(object)$Color_Channel) ## can change
            negs <- setdiff(negs0, negs1)
            if(channel %in% c('Cy3','Grn','G')) {
              return(methylated(object)[negs, , drop=FALSE])
            } else if(channel %in% c('Cy5','Red','R')) {
              return(unmethylated(object)[negs, , drop=FALSE])
            }
          }) #}}}
setMethod('negctls', signature(object="MethyLumiQC",channel='missing'), # {{{
          function(object) {
            channels = list(Cy3='Cy3',Cy5='Cy5')
            lapply(channels, function(channel) negctls(object, channel))
          }) #}}}

setMethod('negctls.SD',signature(object="MethyLumiSet",channel='character'),#{{{
          function(object,channel) return(negctls.SD(QCdata(object),channel)))#}}}
setMethod('negctls.SD',signature(object="MethyLumiM",channel='character'),#{{{
          function(object,channel) return(negctls.SD(controlData(object),channel)))#}}}
setMethod('negctls.SD', signature(object="MethyLumiQC",channel='character'),#{{{
          function(object, channel) {
            negs <- grep('negative', tolower(featureNames(object)))
            if(channel %in% c('Cy3','Grn','G')) {
              return(Cy3.SD(object)[negs,,drop=FALSE])
            } else if(channel %in% c('Cy5','Red','R')) {
              return(Cy5.SD(object)[negs,,drop=FALSE])
            }
          }) #}}}
setMethod('negctls.SD', signature(object="MethyLumiQC",channel='missing'),#{{{
          function(object) {
            channels = list(Cy3='Cy3',Cy5='Cy5')
            lapply(channels, function(channel) negctls.SD(object, channel))
          }) #}}}

setMethod('negctls.stderr', signature(object="MethyLumiSet",channel='character'),#{{{
          function(object, channel) negctls.stderr(QCdata(object), channel))#}}}
setMethod('negctls.stderr', signature(object="MethyLumiQC",channel='character'),#{{{
          function(object, channel) {
            negs <- grep('negative', tolower(featureNames(object)))
            if(channel %in% c('Cy3','Grn','G')) {
              return( (Cy3.SD(object)/sqrt(Cy3.N(object)))[negs,,drop=FALSE] )
            } else if(channel %in% c('Cy5','Red','R')) {
              return( (Cy5.SD(object)/sqrt(Cy5.N(object)))[negs,,drop=FALSE] )
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
              return(Cy5(object)[c(negs,norms),,drop=FALSE])
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
              return(methylated(object)[probes,,drop=FALSE])
            } else if(channel %in% c('Cy5','Red','A','T')) {
              return(unmethylated(object)[probes,,drop=FALSE])
            }
            
          }) #}}}

