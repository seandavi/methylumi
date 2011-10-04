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
  function(object, value, indep=T, minBeads=3) standardGeneric('pval.detect<-')
) # }}}
setReplaceMethod('pval.detect', signature(object="methylData", value="numeric"), function(object, value, indep=T, minBeads=3){ # {{{

  require(matrixStats) 

  if(is(object, 'MethyLumiSet')) stopifnot('QC' %in% slotNames(object))
  if(is(object, 'MethyLumiM')) stopifnot('controlData' %in% slotNames(object))

  channels <- c('Cy3','Cy5')
  names(channels) <- channels

  ## determine which channel each probe is in 
  require(paste(annotation(object),'db',sep='.'), character.only=TRUE)
  colorchan <- toTable(get(paste(annotation(object), 'COLORCHANNEL', sep='')))
  probes <- list( # {{{ by color channel and design type
    Cy3=colorchan$Probe_ID[which(colorchan$Color_Channel=='Grn')],
    Cy5=colorchan$Probe_ID[which(colorchan$Color_Channel=='Red')],
    New=colorchan$Probe_ID[which(colorchan$Color_Channel=='Both')]
  ) # }}}
  rm(colorchan)

  ecdfs <- lapply(sampleNames(object), function(i) { # {{{
    per.channel <- lapply(channels, function(ch) ecdf(negctls(object, ch)[, i]))
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
  for( i in sampleNames(object) ) {
    for( j in names(probes) ) {
      probesj = intersect(probes[[j]], featureNames(object))
      pvals.scratch[ probesj, i ] <- rowMins(1 - dval( probesj, i, j ))
    }
  }

  if(class(object) == 'MethyLumiSet') {
    betas(object) <- pmax(methylated(object),1)/pmax(total.intensity(object),1)
    pvals(object) <- pvals.scratch
    is.na(betas(object))[which(pvals(object) > value, arr.ind=TRUE)] <- TRUE
    if( 'methylated.N' %in% assayDataElementNames(object) &&
        'unmethylated.N' %in% assayDataElementNames(object) ) {
      minN = pmax(assayDataElement(object, 'methylated.N'),
                  assayDataElement(object, 'unmethylated.N'))
      is.na(betas(object))[ which(minN < minBeads, arr.ind=TRUE) ] <- TRUE 
    }
  }
  if(class(object) == 'MethyLumiM') {
    exprs(object)<-log2(pmax(methylated(object),1)/pmax(unmethylated(object),1))
    detection(object) <- pvals.scratch
    is.na(exprs(object))[which(detection(object) > value, arr.ind=TRUE)] <- TRUE
    if( 'methylated.N' %in% assayDataElementNames(object) &&
        'unmethylated.N' %in% assayDataElementNames(object) ) {
      minN = pmax(assayDataElement(object, 'methylated.N'),
                  assayDataElement(object, 'unmethylated.N'))
      is.na(betas(object))[ which(minN < minBeads, arr.ind=TRUE) ] <- TRUE 
    }
  }
  return(object)

}) # }}}

setGeneric('zval.detect<-', # {{{ derp
  function(object, value) standardGeneric('zval.detect<-')
) # }}}
setReplaceMethod('zval.detect', signature(object="methylData", value="numeric"), function(object, value){ # {{{

  require(matrixStats) 
  minBeads = 3 # default
  parallel = ifelse(require('multicore'), TRUE, FALSE)

  if(is(object, 'MethyLumiSet')) stopifnot('QC' %in% slotNames(object))
  if(is(object, 'MethyLumiM')) stopifnot('controlData' %in% slotNames(object))

  channels <- c('Cy3','Cy5')
  names(channels) <- channels
  require(paste(annotation(object),'db',sep='.'), character.only=TRUE)
  colorchan <- toTable(get(paste(annotation(object), 'COLORCHANNEL', sep='')))
  probes <- list( # {{{ by color channel and design type
    Cy3=colorchan$Probe_ID[which(colorchan$Color_Channel=='Grn')],
    Cy5=colorchan$Probe_ID[which(colorchan$Color_Channel=='Red')],
    New=colorchan$Probe_ID[which(colorchan$Color_Channel=='Both')]
  ) # }}}
  rm(colorchan)

  neg.controls <- lapply(channels, function(channel) negctls(object, channel))
  negmeans <- lapply(neg.controls, function(x) colMeans(x, na.rm=T))
  negsds <- lapply(neg.controls, function(x) colSds(x, na.rm=T))
  Mstrong <- (methylated(object) > unmethylated(object))
  prb.i <- (Mstrong*methylated(object)) + ((!Mstrong)*unmethylated(object))
  prb.n <- (Mstrong*methylated.N(object))+((!Mstrong)*unmethylated.N(object))

  # could probably do a better job just keying this against channels()
  if( is(object, 'MethyLumiM')   ) pvalues <- detection(object)
  if( is(object, 'MethyLumiSet') ) pvalues <- pvals(object)
  nchips = dim(object)[2]
  for(ch in names(probes)) {
    p <- probes[[ch]]
    if( ch == 'New' ) {
      pd2 <- function(n) { # {{{
        pv.tmp = c(rep(1.0, length(p)))
        names(pv.tmp) = p 
        p.M = p[which(Mstrong[p,n]==TRUE)]
        p.U = p[which(Mstrong[p,n]==FALSE)]
        neg.M = negmeans[['Cy3']][n]
        neg.U = negmeans[['Cy5']][n]
        div.M = negsds[['Cy3']][n]/sqrt(prb.n[p.M,n])
        div.U = negsds[['Cy5']][n]/sqrt(prb.n[p.U,n])
        pv.tmp[p.M] = (1 - pnorm((prb.i[p.M,n]-neg.M)/div.M))
        pv.tmp[p.U] = (1 - pnorm((prb.i[p.U,n]-neg.U)/div.U))
        return(as.numeric(pv.tmp))
      } # }}}
      if(parallel) {
        pvalues[p,] <- data.matrix(as.data.frame(mclapply(1:nchips, pd2),opt=T))
      } else {
        pvalues[p,] <- sapply(1:nsamples, pd2)
      }
      cat("Computed Illumina-style design II p-values for", nchips, "chips\n")
    } else {
      neg = negmeans[[ch]]
      div = negsds[[ch]]/sqrt(prb.n[p,])
      pd1 <- function(n) { # {{{
        function(n) 1 - pnorm((prb.i[p,n]-neg[n])/div[n])
      } # }}}
      if(parallel) {
        pvalues[p,] <- data.matrix(as.data.frame(mclapply(1:nchips, pd1),opt=T))
      } else {
        pvalues[p,] <- sapply(1:dim(object)[2], pd1)
      }
      cat("Computed Illumina-style", ch, "p-values for", nchips, "chips\n")
    }
  }

  if(class(object) == 'MethyLumiSet') {
    betas(object) <- pmax(methylated(object),1)/pmax(total.intensity(object),1)
    pvals(object) <- pvalues
    is.na(betas(object))[which(pvals(object) > value, arr.ind=TRUE)] <- TRUE
    if( 'methylated.N' %in% assayDataElementNames(object) &&
        'unmethylated.N' %in% assayDataElementNames(object) ) {
      minN = pmin(assayDataElement(object, 'methylated.N'),
                  assayDataElement(object, 'unmethylated.N'))
      is.na(betas(object))[ which(minN < minBeads, arr.ind=TRUE) ] <- TRUE 
    }
  }
  if(class(object) == 'MethyLumiM') {
    exprs(object)<-log2(pmax(methylated(object),1)/pmax(unmethylated(object),1))
    detection(object) <- pvalues
    is.na(exprs(object))[which(detection(object) > value, arr.ind=TRUE)] <- TRUE
    if( 'methylated.N' %in% assayDataElementNames(object) &&
        'unmethylated.N' %in% assayDataElementNames(object) ) {
      minN = pmin(assayDataElement(object, 'methylated.N'),
                  assayDataElement(object, 'unmethylated.N'))
      is.na(exprs(object))[ which(minN < minBeads, arr.ind=TRUE) ] <- TRUE 
    }
  }
  return(object)

}) # }}}

## FIXME: finish Gamma version (Pr(fg > bg|g, a, d, b) = pbeta(a/(a+b), d, g))
if(!isGeneric('gval.detect<-')) setGeneric('gval.detect<-', # {{{ set cutoff
  function(object, value) standardGeneric('gval.detect<-')
) # }}}
setReplaceMethod('gval.detect', signature(object="methylData", value="numeric"), function(object, value){ # {{{

  stop('Gamma (empirical) detection calls are not yet implemented')

  #require('matrixStats') # needed for colSds
  #if(is(object, 'MethyLumiSet')) stopifnot('QC' %in% slotNames(object))
  #if(is(object, 'MethyLumiM')) stopifnot('controlData' %in% slotNames(object))
  #channels <- c('Cy3','Cy5')
  #names(channels) <- channels
  
  #if(annotation(object) == 'IlluminaHumanMethylation450k') {
    #warning("The cy3() and cy5() methods below are the cause of 450k problems")
  #}
  #probes <- list(Cy3=cy3(object), Cy5=cy5(object))

  ### FIXME: for 'design II' probes, either normalize them all,
  ###        or else have 'cy3()' and 'cy5()' correctly pull by design
  ###
  #neg.controls <- lapply(channels, function(channel) negctls(object, channel))
  #negmeans <- lapply(neg.controls, function(x) colMeans(x, na.rm=T))
  #negsds <- lapply(neg.controls, function(x) colSds(x, na.rm=T))

  ### Make this work properly with Gamma deconvolution and/or Lumi bgcorrect too
  ###
  #zvalues <- pvals(object) # will end up flipping this
  #if(indep) {
    #Mstrong <- methylated(object) > unmethylated(object)
    #prb.i <- (Mstrong*methylated(object)) + ((!Mstrong)*unmethylated(object))
    #prb.n <- (Mstrong*methylated.N(object))+((!Mstrong)*unmethylated.N(object))
  #} else {            
    #prb.n <- (methylated.N(object)+unmethylated.N(object))/2
    #prb.i <- ((methylated(object)*methylated.N(object)) +
              #(unmethylated(object)*unmethylated.N(object)))/(prb.n*2)
  #}
  #registered <- pmin(methylated.N(object), unmethylated.N(object)) > minBeads

  ## could probably do a better job just keying this against channels()
  #for(ch in names(channels)) {
    #p <- probes[[ch]]
    #zvalues[p,] <- sapply(1:length(sampleNames(object)), function(n) {
      #pnorm((prb.i[p,n]-negmeans[[ch]][n])/((negsds[[ch]][n])/sqrt(prb.n[p,n])))
    #})
  #}

  #if(class(object) == 'MethyLumiSet') {
    #pvals(object) <- 1 - zvalues 
    #is.na(exprs(object))[which(pvals(object) > value, arr.ind=TRUE)] <- TRUE
    #is.na(exprs(object))[which(registered == FALSE, arr.ind=TRUE)] <- TRUE
  #}
  #if(class(object) == 'MethyLumiM') {
    #detection(object) <- 1 - zvalues
    #is.na(exprs(object))[which(detection(object) > value, arr.ind=TRUE)] <- TRUE
    #is.na(exprs(object))[which(registered == FALSE, arr.ind=TRUE)] <- TRUE
  #}
  #return(object)
}) # }}}
