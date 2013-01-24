if(is.null(getGeneric('pval.detect'))) setGeneric('pval.detect', # {{{ cutoff
  function(object) standardGeneric('pval.detect')) # }}}
setMethod('pval.detect', signature(object="methylData"), function(object){ # {{{
  if( is(object, 'MethyLumiM') ) {
    round(max(pvals(object)[which(!is.na(exprs(object)))]),2)
  } else { 
    round(max(pvals(object)[which(!is.na(betas(object)))]),2)
  }
}) # }}}

if(is.null(getGeneric('pval.detect<-'))) setGeneric('pval.detect<-', # {{{ 
  function(object, ..., value) standardGeneric('pval.detect<-')
) # }}}
setReplaceMethod('pval.detect', signature(object="methylData", value="numeric"), function(object, ..., value){ # {{{

  require(matrixStats)
  if(is(object, 'MethyLumiSet')) stopifnot('QC' %in% slotNames(object))
  if(is(object, 'MethyLumiM')) stopifnot('controlData' %in% slotNames(object))
  channels <- c(Cy3='Cy3',Cy5='Cy5')

  ## determine which channel each probe is in 
  require(paste(annotation(object),'db',sep='.'), character.only=TRUE)
  colorchan <- toTable(get(paste(annotation(object), 'COLORCHANNEL', sep='')))
  probes <- list( # {{{ by color channel and design type
    Cy3=colorchan$Probe_ID[which(colorchan$Color_Channel=='Grn')],
    Cy5=colorchan$Probe_ID[which(colorchan$Color_Channel=='Red')],
    New=colorchan$Probe_ID[which(colorchan$Color_Channel=='Both')]
  ) # }}}
  rm(colorchan)

  # instead of a normal approximation, use the ECDF of the negative controls
  # interestingly, this can come in handy when dealing with FFPE samples 
  ecdfs <- lapply(sampleNames(object), function(i) { # {{{
    per.channel <- lapply(channels, 
                          function(ch) { # {{{ ECDF of true negative controls
                            ids <- rownames(negctls(object, ch))
                            color <- fData(QCdata(object))[ids,'Color_Channel']
                            keep <- which(color != '-99')
                            background <- negctls(object, ch)[keep, i, drop=F]
                            ecdf(background)
                          } # }}}
    )
    names(per.channel) <- names(channels)
    return(per.channel)
  }) # }}}
  names(ecdfs) <- sampleNames(object)

  if( is(object, 'MethyLumiM') ) pvals.scratch <- detection(object)
  if( is(object, 'MethyLumiSet') ) pvals.scratch <- pvals(object)

  dval <- function(probes, subject, type, allele) { # {{{
    probes = intersect(featureNames(object), probes)
    if( type == 'New' ) {
      cbind(M=ecdfs[[subject]][['Cy3']](methylated(object)[probes,subject]),
            U=ecdfs[[subject]][['Cy5']](unmethylated(object)[probes,subject]))
    } else {
      cbind(M=ecdfs[[subject]][[type]](methylated(object)[probes,subject]),
            U=ecdfs[[subject]][[type]](unmethylated(object)[probes,subject]))
    }
  } # }}}

  ## FIXME: trivially parallelizable, or farm out to C++?
  for( i in sampleNames(object) ) { # {{{
    for( j in names(probes) ) {
      probesj = intersect(probes[[j]], featureNames(object))
      pvals.scratch[ probesj, i ] <- rowMins(1 - dval( probesj, i, j ))
    }
  } # }}}
  if(class(object) == 'MethyLumiSet') { # {{{
    betas(object) <- pmax(methylated(object),1)/pmax(total.intensity(object),1)
    pvals(object) <- pvals.scratch
    is.na(betas(object))[which(pvals(object) > value, arr.ind=TRUE)] <- TRUE
  } # }}}
  if(class(object) == 'MethyLumiM') { # {{{
    exprs(object)<-log2(pmax(methylated(object),1)/pmax(unmethylated(object),1))
    detection(object) <- pvals.scratch
    is.na(exprs(object))[which(detection(object) > value, arr.ind=TRUE)] <- TRUE
  } # }}}
  return(object)

}) # }}}
