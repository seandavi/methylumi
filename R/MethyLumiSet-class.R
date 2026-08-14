# $HeadUrl$
# $Id$
##=================================================
## Define MethyLumiSet class:
setClass('MethyLumi', contains='eSet')
setClass('MethyLumiQC', contains='MethyLumi')
setClassUnion("QCDataOrNULL",c('NULL',"MethyLumiQC"))
setClassUnion("methylData",c('MethyLumi','ExpressionSet'))

setClass('MethyLumiSet', # {{{
         representation(QC="QCDataOrNULL",
                        # OOB="OOBDataOrNULL",
                        history='data.frame'),
         prototype=list(QC=NULL,
           history=data.frame(
             submitted   = I(vector()),
             finished    = I(vector()),
             command     = I(vector())
             )),
         contains='MethyLumi') # }}}
setMethod('initialize', 'MethyLumiSet', # {{{
          function(.Object,
                   assayData=assayDataNew(
                     betas       = betas,
                     ...),
                   phenoData   = annotatedDataFrameFrom(assayData,byrow=FALSE),
                   featureData = annotatedDataFrameFrom(assayData,byrow=TRUE),
                   experimentData = new("MIAME"),
                   annotation= character(),
                   betas = new("matrix"),
				   ...
                   )
                   {
            callNextMethod(.Object,
                           assayData=assayData,
                           phenoData=phenoData,
                           featureData=featureData,
                           experimentData=experimentData,
                           annotation=annotation)
          }
          ) # }}}
setMethod('initialize', 'MethyLumiQC', # {{{
          function(.Object,
                   assayData=assayDataNew(
                     ...),
                   phenoData   = annotatedDataFrameFrom(assayData,byrow=FALSE),
                   featureData = annotatedDataFrameFrom(assayData,byrow=TRUE),
                   experimentData = new("MIAME"),
                   annotation= character(),
                   betas = new("matrix")
                   )
                   {
            callNextMethod(.Object,
                           assayData=assayData,
                           phenoData=phenoData,
                           featureData=featureData,
                           experimentData=experimentData,
                           annotation=annotation)
          }
          ) # }}}
setValidity("MethyLumiSet", function(object) { # {{{
    msg <- Biobase:::validMsg(NULL, Biobase:::isValidVersion(object, "eSet"))
    msg <- Biobase:::validMsg(msg, assayDataValidMembers(assayData(object), c("betas")))
    if (is.null(msg)) TRUE else msg
}) # }}}
setValidity("MethyLumiQC", function(object) { # {{{
    msg <- Biobase:::validMsg(NULL, Biobase:::isValidVersion(object, "eSet"))
#    msg <- Biobase:::validMsg(msg, assayDataValidMembers(assayData(object), c("avgsignal")))
    if (is.null(msg)) TRUE else msg
}) # }}}

if(!isGeneric('betas')) { # {{{
  setGeneric('betas', function(object) standardGeneric('betas')) 
} # }}}
setMethod("betas", signature(object="MethyLumiSet"), # {{{
          function(object) assayDataElement(object,"betas")) # }}}
setGeneric('betas<-', # {{{
    function(object,value) standardGeneric('betas<-')) # }}}
setReplaceMethod("betas", signature(object="MethyLumiSet",value="matrix"), # {{{
function(object, value) assayDataElementReplace(object, "betas", value)) # }}}

setMethod('exprs', signature(object='MethyLumiSet'), function(object) { # {{{
  log2( pmax(pmin(betas(object), 0.99999), 0.000001) / 
        pmax(pmin(1-betas(object), 0.99999), 0.000001) )
}) # }}}
setGeneric('mvals', # {{{
    function(object) standardGeneric('mvals')) # }}}
setMethod('mvals', signature(object='MethyLumiSet'), function(object) { # {{{
  log2( pmax(pmin(betas(object), 0.99999), 0.000001) / 
        pmax(pmin(1-betas(object), 0.99999), 0.000001) )
}) # }}}

setGeneric('pvals',
           function(object) standardGeneric('pvals'))
setMethod("pvals", signature(object="MethyLumi"), # {{{
          function(object) assayDataElement(object,"pvals")) # }}}
setReplaceMethod("pvals", signature(object="MethyLumi",value="matrix"), # {{{
                 function(object, value) assayDataElementReplace(object, "pvals", value)) # }}}

setMethod("QCdata", signature(object="MethyLumiSet"), # {{{
          function(object) object@QC) # }}}
setReplaceMethod("QCdata", signature(object="MethyLumiSet",value="MethyLumiQC"), function(object, value) { # {{{
    object@QC <- value
    return(object)
  }) # }}}
setMethod("controlData", signature(object="MethyLumiSet"), # {{{
          function(object) return(object@QC)) # }}}
setReplaceMethod("controlData", signature(object="MethyLumiSet",value="MethyLumiQC"), function(object, value) { # {{{
  object@QC <- value
  return(object)
}) # }}} 

setMethod("getHistory",signature(object="MethyLumiSet"), function(object) {# {{{
    object@history
}) # }}}

setMethod("unmethylated",signature(object="MethyLumiSet"),function(object){# {{{
            return(assayDataElement(object,"unmethylated"))
          }) # }}}
setReplaceMethod("unmethylated", signature(object="MethyLumiSet",value="matrix"), function(object,value) { # {{{
            assayDataElementReplace(object,"unmethylated",value)
          }) # }}}
setMethod("methylated",signature(object="MethyLumiSet"),function(object) { # {{{
            return(assayDataElement(object,"methylated"))
            }) # }}}
setReplaceMethod("methylated", signature(object="MethyLumiSet",value="matrix"), function(object,value) { # {{{
            assayDataElementReplace(object,"methylated",value)
          }) # }}}

setMethod("unmethylated.OOB", signature(object="MethyLumiSet"),function(object){ # {{{
            return(assayDataElement(object,"unmethylated.OOB"))
          }) # }}}
setReplaceMethod("unmethylated.OOB", signature(object="MethyLumiSet",value="matrix"), function(object,value) { # {{{ 
            assayDataElementReplace(object,"unmethylated.OOB",value)
          }) # }}}
setMethod("methylated.OOB", signature(object="MethyLumiSet"), function(object) { # {{{ 
            return(assayDataElement(object,"methylated.OOB"))
            }) # }}}
setReplaceMethod("methylated.OOB", signature(object="MethyLumiSet",value="matrix"), function(object,value) { # {{{
            assayDataElementReplace(object,"methylated.OOB",value)
          }) # }}}

setMethod("show",signature(object="MethyLumiSet"), function(object) { # {{{ 
	cat('\nObject Information:\n')
	callNextMethod()
	cat('Major Operation History:\n')
	print(getHistory(object)) 
}) # }}}
setMethod("summary",signature(object="MethyLumi"), function(object, ...) { # {{{
	show(object)
}) # }}}

if(is.null(getGeneric('intensities.OOB'))) { # {{{
  setGeneric('intensities.OOB',function(x, channel) {
    standardGeneric('intensities.OOB')
  })
} # }}}
setMethod("intensities.OOB",signature(x="MethyLumiSet", channel="character"), function(x, channel) { # {{{
    if(!('COLOR_CHANNEL' %in% fvarLabels(x))) { 
      if(annotation(x) == 'IlluminaHumanMethylation27k') {
        fData(x)$COLOR_CHANNEL = mget(featureNames(x), 
                                      IlluminaHumanMethylation27kCOLORCHANNEL)
      }
      if(annotation(x) == 'IlluminaHumanMethylation450k') {
        fData(x)$COLOR_CHANNEL = mget(featureNames(x), 
                                      IlluminaHumanMethylation450kCOLORCHANNEL)
      }
    } 
    if(channel == 'Cy3') { 
      probes = which(fData(x)$COLOR_CHANNEL == 'Red') 
    } else if(channel == 'Cy5') { 
      probes = which(fData(x)$COLOR_CHANNEL == 'Grn')
    } 
    return(rbind( assayDataElement(x, 'methylated.OOB')[probes,],
                  assayDataElement(x, 'unmethylated.OOB')[probes,] ) )
}) # }}}
setMethod("intensities.OOB",signature(x="MethyLumiSet", channel="missing"),# {{{
  function(x) lapply(list(Cy3='Cy3',Cy5='Cy5'),function(y) intensities.OOB(x,y))
) # }}}

if(is.null(getGeneric('intensities.OOB.allelic'))) { # {{{
  setGeneric('intensities.OOB.allelic',function(x, channel, allele) {
    standardGeneric('intensities.OOB.allelic')
  })
} # }}}
setMethod("intensities.OOB.allelic",signature(x="MethyLumiSet", channel="character", allele="character"), function(x, channel, allele) { # {{{
    if(!('COLOR_CHANNEL' %in% fvarLabels(x))) { 
      if(annotation(x) == 'IlluminaHumanMethylation27k') { # {{{
        fData(x)$COLOR_CHANNEL = mget(featureNames(x), 
                                      IlluminaHumanMethylation27kCOLORCHANNEL)
      } # }}}
      if(annotation(x) == 'IlluminaHumanMethylation450k') { # {{{
        fData(x)$COLOR_CHANNEL = mget(featureNames(x), 
                                      IlluminaHumanMethylation450kCOLORCHANNEL)
      } # }}}
    } 
    element = paste(allele, 'OOB', sep='.')
    if(channel == 'Cy3') probes = which(fData(x)$COLOR_CHANNEL == 'Red') 
    if(channel == 'Cy5') probes = which(fData(x)$COLOR_CHANNEL == 'Grn') 
    return( assayDataElement(x, element)[probes, ])
}) # }}}
setMethod("intensities.OOB.allelic",signature(x="MethyLumiSet", channel="missing", allele="missing"),# {{{
  function(x) {
    lapply(list(Cy3='Cy3',Cy5='Cy5'), function(y) {
      lapply(list(M='methylated', U='unmethylated'), function(z) {
        intensities.OOB.allelic(x, y, z)
      })
    })
  }) # }}}

if(is.null(getGeneric('intensities.IB'))) { # {{{
  setGeneric('intensities.IB',function(x, channel) {
    standardGeneric('intensities.IB')
  })
} # }}}
setMethod("intensities.IB",signature(x="MethyLumiSet", channel="character"), function(x, channel) { # {{{
    if(!('COLOR_CHANNEL' %in% fvarLabels(x))) { 
      if(annotation(x) == 'IlluminaHumanMethylation27k') {
        fData(x)$COLOR_CHANNEL = mget(featureNames(x), 
                                      IlluminaHumanMethylation27kCOLORCHANNEL)
      }
      if(annotation(x) == 'IlluminaHumanMethylation450k') {
        fData(x)$COLOR_CHANNEL = mget(featureNames(x), 
                                      IlluminaHumanMethylation450kCOLORCHANNEL)
      }
    } 
    if(channel == 'Cy3') { 
      probes = which(fData(x)$COLOR_CHANNEL == 'Grn') 
    } else if(channel == 'Cy5') { 
      probes = which(fData(x)$COLOR_CHANNEL == 'Red')
    } 
    return(rbind( assayDataElement(x, 'methylated')[probes,],
                  assayDataElement(x, 'unmethylated')[probes,] ) )
}) # }}}
setMethod("intensities.IB",signature(x="MethyLumiSet", channel="missing"), # {{{
  function(x) lapply(list(Cy3='Cy3',Cy5='Cy5'),function(y) intensities.IB(x,y))
) # }}}

if(is.null(getGeneric('intensities.M'))) { # {{{
  setGeneric('intensities.M',function(x, channel) {
    standardGeneric('intensities.M')
  })
} # }}} 
setMethod("intensities.M",signature(x="MethyLumiSet", channel="character"),# {{{
  function(x, channel) {
    if(!('COLOR_CHANNEL' %in% fvarLabels(x))) { 
      if(annotation(x) == 'IlluminaHumanMethylation27k') {
        fData(x)$COLOR_CHANNEL = mget(featureNames(x), 
                                      IlluminaHumanMethylation27kCOLORCHANNEL)
      }
      if(annotation(x) == 'IlluminaHumanMethylation450k') {
        fData(x)$COLOR_CHANNEL = mget(featureNames(x), 
                                      IlluminaHumanMethylation450kCOLORCHANNEL)
      }
    } 
    if(channel == 'Cy3') { 
      probes = which(fData(x)$COLOR_CHANNEL == 'Grn') 
    } else if(channel == 'Cy5') { 
      probes = which(fData(x)$COLOR_CHANNEL == 'Red')
    } 
    return(assayDataElement(x, 'methylated')[probes,]) 
}) # }}}
setMethod("intensities.M",signature(x="MethyLumiSet", channel="missing"), # {{{
  function(x) lapply(list(Cy3='Cy3',Cy5='Cy5'),function(y) intensities.M(x,y))
) # }}}

if(is.null(getGeneric('intensities.U'))) { # {{{
  setGeneric('intensities.U',function(x, channel) {
    standardGeneric('intensities.U')
  })
} # }}}
setMethod("intensities.U",signature(x="MethyLumiSet", channel="character"),# {{{
  function(x, channel) {
    if(!('COLOR_CHANNEL' %in% fvarLabels(x))) { 
      if(annotation(x) == 'IlluminaHumanMethylation27k') {
        fData(x)$COLOR_CHANNEL = mget(featureNames(x), 
                                      IlluminaHumanMethylation27kCOLORCHANNEL)
      }
      if(annotation(x) == 'IlluminaHumanMethylation450k') {
        fData(x)$COLOR_CHANNEL = mget(featureNames(x), 
                                      IlluminaHumanMethylation450kCOLORCHANNEL)
      }
    } 
    if(channel == 'Cy3') { 
      probes = which(fData(x)$COLOR_CHANNEL == 'Grn') 
    } else if(channel == 'Cy5') { 
      probes = which(fData(x)$COLOR_CHANNEL == 'Red')
    } 
    return(assayDataElement(x, 'unmethylated')[probes,]) 
}) # }}}
setMethod("intensities.U",signature(x="MethyLumiSet", channel="missing"), # {{{
  function(x) lapply(list(Cy3='Cy3',Cy5='Cy5'),function(y) intensities.U(x,y))
) # }}}

#if(is.null(getGeneric('boxplot'))) {  # {{{
  #setGeneric('boxplot',function(x,...) standardGeneric('boxplot'))
#} # }}}
#setMethod("boxplot",signature(x="MethyLumiSet"), function(x, range=0, main, logMode=TRUE, ...) { # {{{ 
  #tmp <- description(x)
  #if (missing(main) && (is(tmp, "MIAME")))
    #main <- tmp@title
  #exprs <- exprs(x)
  #if (nrow(x) > 5000) {
			#index <- seq(1, nrow(x), len=5000)
              #} else {
		#index <- 1:nrow(x)
              #}
  #if (logMode & max(exprs(x), na.rm=TRUE) > 50) {
    #exprs <- log2(exprs)
  #} 
  #dataMatrix <- exprs[index,]
  #labels <- colnames(dataMatrix)
  #if (is.null(labels)) labels <- as.character(1:ncol(dataMatrix))
  ### set the margin of the plot
  #mar <- c(max(nchar(labels))/2 + 4.5, 5, 5, 3)
  #old.mar <- par('mar')
  #old.xaxt <- par('xaxt')
  #par(xaxt='n')
  #par(mar=mar)
  #boxplot(dataMatrix ~ col(dataMatrix), main=main, range=range, xlab='', ylab='amplitude', ...)
  #par(xaxt='s')
  #axis(1, at=1:ncol(dataMatrix), labels=labels, tick=TRUE, las=2)
  #par(mar=old.mar)
  #par(xaxt=old.xaxt)
#}) # }}}

if (is.null(getGeneric("pairs"))) { # {{{
  setGeneric("pairs", function(x,...) standardGeneric("pairs"))
} # }}}
setMethod("pairs", signature(x="MethyLumiSet"), function(x,...,logMode=FALSE,maxpairs=5,fold=0.1) { # {{{
    
	upperPanel <- function(x, y, fold=fold) {
		if (length(x) > 3000) {
			ind <- sample(1:length(x), 3000)
			x <- x[ind]; y <- y[ind]
		}
		points(x, y)
		abline(0, 1, col="red", lty=1)
		if (logMode) {
			abline(log2(fold), 1, col="green", lty=2)
			abline(log2(1/fold), 1, col="green", lty=2)
		} else {
			abline(fold, 1, col="green", lty=2)
			abline(-fold, 1, col="green", lty=2)
		}
	}

	lowerPanel <- function(x, y, cex=1.44, fold=fold) {
		if (logMode) {
			up <- length(which((x-y) > log2(fold)))
			down <- length(which((y-x) > log2(fold)))
		} else {
			up <- length(which((x-y) > fold))
			down <- length(which((y-x) > fold))
		}
		ex <- par("fin")[1]*0.9
		txt <- paste("Cor =", as.character(round(cor(x,y),2)),"\n")
		txt <- paste(txt, up, " (> ", fold, ", up)\n", sep="")
		txt <- paste(txt, down, " (> ", fold, ", down)\n", sep="")
		text(mean(range(x)), mean(range(x)), labels=txt, cex=ex)
	}

	## put histograms on the diagonal
	diagPanel <- function(x, ...) {
	    usr <- par("usr"); on.exit(par(usr))
	    par(usr = c(usr[1:2], 0, 1.5) )
	    h <- hist(x, plot = FALSE)
	    breaks <- h$breaks; nB <- length(breaks)
	    y <- h$counts; y <- y/max(y)
	    rect(breaks[-nB], 0, breaks[-1], y, col="cyan", ...)
	}
        par(ask=TRUE)
        for(sampstart in seq(1,ncol(x)-maxpairs,maxpairs)) {
          if(logMode & (max(exprs(x), na.rm=TRUE) > 50)) {
            pairs(log2(exprs(x)[,sampstart:(sampstart+maxpairs-1)]),upper.panel=upperPanel, diag.panel=diagPanel, 
                  lower.panel=lowerPanel, ...)
          } else {
            pairs(exprs(x)[,sampstart:(sampstart+maxpairs-1)],upper.panel=upperPanel, diag.panel=diagPanel, 
                  lower.panel=lowerPanel, ...)
          }
        }
        par(ask=FALSE)
})  # }}}

#' Plot the sample intensities
#'
#' The Illumina methylation platforms all show a significant dye bias. The
#' `plotSampleIntensities` method shows the density plots for the two channels,
#' allowing direct visualization of the effect.
#'
#' @param x An object of class MethyLumi or a subclass.
#' @param beta.cuts Cutoffs for low and high beta values.
#' @param s Sample number to plot.
#' @keywords hplot
#' @examples
#' data(mldat)
#' plotSampleIntensities(mldat, s = 1)
#' @name plotSampleIntensities
#' @aliases plotSampleIntensities
#' @usage plotSampleIntensities(x,beta.cuts,s)
NULL

if (is.null(getGeneric("plotSampleIntensities"))) { # {{{
  setGeneric("plotSampleIntensities", function(x,beta.cuts=c(0.2,0.8),s=1) {
      standardGeneric("plotSampleIntensities")
    }) 
} # }}}
myDensity <- function(x) { # {{{
  x <- x[!is.na(x)]
  return(density(x))
} # }}}
setMethod('plotSampleIntensities', signature(x='MethyLumiSet'), function(x, beta.cuts=c(0.2,0.8),s=1) { # {{{
  message("The purpose of this method is better served by diagnostics()")
  cy3h <- myDensity(unmethylated(x)[betas(x)[,s]<beta.cuts[1],s])
  cy3l <- myDensity(unmethylated(x)[betas(x)[,s]>beta.cuts[2],s])
  cy5l <- myDensity(methylated(x)[betas(x)[,s]<beta.cuts[1],s])
  cy5h <- myDensity(methylated(x)[betas(x)[,s]>beta.cuts[2],s])
  ymax <- max(c(cy5h$y,cy3h$y,cy3l$y,cy5l$y))
  xmax <- max(c(cy5h$x,cy3h$x,cy3l$x,cy5l$x))
  xmin <- min(c(cy5h$x,cy3h$x,cy3l$x,cy5l$x))
  cy5l$y <- cy5l$y/(max(cy5l$y)/ymax)
  cy5h$y <- cy5h$y/(max(cy5h$y)/ymax)
  cy3h$y <- cy3h$y/(max(cy3h$y)/ymax)
  cy3l$y <- cy3l$y/(max(cy3l$y)/ymax)
  plot(cy3h,ylim=c(0,ymax),xlim=c(xmin,xmax),col='green',axes=FALSE,xlab="Intensity",
       ylab="Relative density",
       main=sprintf("Intensities at High(%4.2f) and Low(%4.2f) betas",beta.cuts[2],beta.cuts[1]),
       sub=sprintf("Sample %d",s))
  lines(cy5l,col='red')
  lines(cy3l,col='green')
  lines(cy5h,col='red')
  box()
  axis(1)
})  # }}}

setMethod("hist",signature(x="MethyLumiSet"),function(x,...) { # {{{
  samples = dim(x)[2]
  extra = ifelse(exists('extra'), as.character(extra), '')
  per.side = ceiling(sqrt(samples))
  if(per.side > 5) warning("This plot is not going to be easy to read...")
  if(samples == 8) { # useful special case for me
    par(mfrow=c(2, 4))
  } else { 
    par(mfrow=c(per.side, per.side))
  }
  for(i in seq_along(sampleNames(x))) {
    hist(betas(x)[,i],xlab="Beta",main=paste(sampleNames(x)[i],extra),breaks=99)
  }
}) # }}}
setMethod("hist",signature(x="MethyLumiQC"),function(x,...) { # {{{
  samples = dim(x)[2]
  if(samples > 5) stop("Too many samples, choose a subset for a decent plot")
  else par(mfrow=c(samples,2))
  max.neg = max( max(negctls(x, 'Cy5')), max(negctls(x, 'Cy3')) )
  xl = c(0, max.neg) # x limits
  for(i in seq_along(sampleNames(x))) {
    hist(negctls(x,'Cy3')[,i], breaks=50, col='green', border='green',
         main=paste(sampleNames(x)[i], "negative controls"), xlim=xl,
         xlab='Nonspecific Cy3 fluorescence from negative controls', ...)
    hist(negctls(x,'Cy5')[,i], breaks=50, col='red', border='red',
         main=paste(sampleNames(x)[i], "negative controls"), xlim=xl,
         xlab='Nonspecific Cy5 fluorescence from negative controls', ...)
  }
}) # }}}

setMethod("[", "MethyLumiSet", function(x, i, j, ..., drop = FALSE) { # {{{

  history.submitted <- as.character(Sys.time())

  ## eSet's "[" method drops the QC slot, so keep a reference before calling it
  ## and reattach below. Subsetting it afterwards, as the previous code did, is
  ## a no-op: x@QC is already NULL by then, and NULL[, j, drop=FALSE] returns
  ## NULL silently rather than erroring, which is how this went unnoticed. See
  ## issue #21.
  qc <- if ('QC' %in% slotNames(x)) x@QC else NULL

  x <- callNextMethod()

  ## QC probes are control probes: they live in their own feature space and are
  ## not indexed by i. Only sample selection applies to them.
  if (!is.null(qc)) x@QC <- if (missing(j)) qc else qc[, j, drop=FALSE]

  ddim <- dim(x)
  if (!missing(i) & !missing(j)) {
    history.command <- paste('Subset of',ddim[1],'features &',ddim[2],'samples')
  } else if (!missing(i)) {
    history.command <- paste('Subset of', ddim[1], 'features.')
  } else if (!missing(j)) {
    history.command <- paste('Subset of', ddim[2], 'samples.')
  } else {
    history.command <- 'Subset with no indices.'
  }

  # history tracking
  history.finished <- as.character(Sys.time())
  x@history<- rbind(x@history, data.frame(submitted=history.submitted,finished=history.finished,command=history.command))
  
  return(x)
}) # }}}

if(is.null(getGeneric("combine"))) { # {{{
 	setGeneric("combine", function(x, y, ...)	standardGeneric("combine"))
} # }}}
.combine.methylumiQC <- function(x,y) { # {{{

  if (class(x)!=class(y)) { # {{{
    stop(paste("objects must be the same class, but are ", 
               class(x), ", ", class(y), sep=""))
  } # }}}
  sampleNames(x) <- sampleNames(assayData(x)) ## fixes an old glitch
  sampleNames(y) <- sampleNames(assayData(y)) ## fixes an old glitch
  if(!identical(sort(featureNames(x)), sort(featureNames(y)))) { # {{{
    feats <- intersect(featureNames(x), featureNames(y))
    if(length(feats) > 0) {
      message('Using the shared subset of ',length(feats),' QCdata features...')
      x <- x[feats,]
      featureData(x) <- combine(featureData(x[feats,]), featureData(y[feats,]))
      assayData(x) <- combine(assayData(x[feats,]), assayData(y[feats,]))
    } else {
      stop('No QCdata features in common were found!')
    } # }}}
  } else { # {{{
    assayData(x) <- combine(assayData(x), assayData(y))
    featureData(x) <- combine(featureData(x), featureData(y))
  } # }}}
  phenoData(x) <- combine(phenoData(x), phenoData(y))
  protocolData(x) <- combine(protocolData(x), protocolData(y))
  experimentData(x) <- combine(experimentData(x), experimentData(y))
	return(x)

}  # }}}
setMethod("combine", signature=c(x="MethyLumiQC", y="MethyLumiQC"), function(x,y) { # {{{
  .combine.methylumiQC(x,y)
}) # }}}
.combine.methylumiSets <- function(x,y) { # {{{
  if (class(x)!=class(y)) stop(paste("Classes differ:",class(x),"!=",class(y)))
  sampleNames(x) <- sampleNames(assayData(x)) ## fixes an old glitch
  sampleNames(y) <- sampleNames(assayData(y)) ## fixes an old glitch
  history.submitted <- as.character(Sys.time()) ## start logging here
  if(!identical(sort(featureNames(x)), sort(featureNames(y)))) { # {{{
    feats <- intersect(featureNames(x), featureNames(y))
    if(length(feats) > 0) {
      message('Using the shared subset of', length(feats), 'features...')
      x <- x[feats,]
      assayData(x) <- combine(assayData(x[feats,]), assayData(y[feats,]))
      featureData(x) <- combine(featureData(x[feats,]), featureData(y[feats,]))
    } # }}}
  } else { # {{{
    assayData(x) <- combine(assayData(x), assayData(y))
    featureData(x) <- combine(featureData(x), featureData(y))
  } # }}}
  pdat <- try(combine(pData(x), pData(y)), silent=TRUE)
  if(inherits(pdat, "try-error")) {
    message('Attempting to coerce phenoData columns into usable types...')
    sharedCols <- intersect(colnames(pData(x)), colnames(pData(y)))
    for(i in sharedCols) {
      class(pData(y)[,i]) <- class(c(pData(x)[,i], pData(y)[,i]))
      class(pData(x)[,i]) <- class(pData(y)[,i])
    }
    pdat <- try(combine(pData(x), pData(y)), silent=TRUE)
    if(inherits(pdat, "try-error")) browser()
  } 
  pData(x) <- pdat
  protocolData(x) <- combine(protocolData(x), protocolData(y))
  experimentData(x) <- combine(experimentData(x), experimentData(y))
	if (is.null(x@QC) | is.null(y@QC)) { # {{{
    message('Dropped control probes: any(is.null(QCdata(x),QCdata(y))) == TRUE')
  } else {
    QCdata(x) = .combine.methylumiQC(QCdata(x), QCdata(y))
  }# }}}
  x@history <- rbind(x@history, y@history)
  x@history <- rbind(x@history,  data.frame(
    submitted=history.submitted,
    finished=as.character(Sys.time()),
    command=paste('Added',dim(y)[2],'samples (',dim(x)[2],'unique samples).')
  ))
	return(x)
}  # }}}
setMethod("combine",signature=c(x="MethyLumiSet",y="MethyLumiSet"),function(x,y){ # {{{
  .combine.methylumiSets(x,y)
}) # }}}

if(is.null(getGeneric("combine27k450k"))) { # {{{
 	setGeneric("combine27k450k", function(x, y, ...) {
    if(length(list(...)) > 0) callGeneric(x, do.call(callGeneric, list(y, ...)))
    else standardGeneric("combine27k450k")
  })
} # }}}
.combine.methylumiQC.27k.450k <- function(x,y,...) { # {{{

  message("This method needs polishing -- should subset and/or impute both...")
  if (class(x)!=class(y)) { # {{{
    stop(paste("objects must be the same class, but are ", 
               class(x), ", ", class(y), sep=""))
  } # }}}
  if (any(sort(featureNames(x)) != sort(featureNames(y)))) { # {{{
    stop('The two data sets have different row names!')
  } # }}}
  assayData(x) <- combine(assayData(x), assayData(y))
  phenoData(x) <- combine(phenoData(x), phenoData(y))
  featureData(x) <- combine(featureData(x), featureData(y))
  protocolData(x) <- combine(protocolData(x), protocolData(y))
  experimentData(x) <- combine(experimentData(x), experimentData(y))
	return(x)

}  # }}}
setMethod("combine27k450k", signature=c(x="MethyLumiSet", y="MethyLumiSet"), function(x,y) { # {{{
  if (class(x)!=class(y)) { # {{{
    stop(paste("objects must be the same class, but are ", 
               class(x), ", ", class(y), sep=""))
  } # }}}

  history.submitted <- as.character(Sys.time())
  platform.x = gsub('IlluminaHumanMethylation', '', annotation(x))
  platform.y = gsub('IlluminaHumanMethylation', '', annotation(y))
  if(any(fvarLabels(x)[-1] %in% fvarLabels(y)[-1]) && platform.x != platform.y){
    fvarLabels(x) = paste(fvarLabels(x), platform.x, sep='.')
    fvarLabels(y) = paste(fvarLabels(y), platform.y, sep='.')
  }
  x = subsetCommonProbes(x)
  y = subsetCommonProbes(y)
  n.x = dim(x)[2]
  n.y = dim(y)[2]
  n.xy = length(unique(sampleNames(x),sampleNames(y)))
  if(n.xy<(n.x+n.y)) {
    sampleNames(x) = paste(sampleNames(x), platform.x, sep='.')
    sampleNames(y) = paste(sampleNames(y), platform.y, sep='.')
    n.xy = length(unique(sampleNames(x),sampleNames(y)))
  }
  assayData(x) <- combine(assayData(x), assayData(y))
  phenoData(x) <- combine(phenoData(x), phenoData(y))
  featureData(x) <- combine(featureData(x), featureData(y))
  protocolData(x) <- combine(protocolData(x), protocolData(y))
  experimentData(x) <- combine(experimentData(x), experimentData(y))
	if (!is.null(x@QC) | !is.null(y@QC)) { # {{{
    QCdata(x)
  } # }}}
#  if (!is.null(x@OOB)|!is.null(y@OOB)) { # {{{
#    OOB(x) = .combine.methylumiOOB(OOB(x), OOB(y))
#  } # }}}
  history.finished <- as.character(Sys.time())
  history.command <- paste('Combined common probes from', n.x, 'samples',
                           'with', n.y, 'additional samples', 
                           paste('(', n.xy, ' unique)', sep=''))
  x@history <- rbind(x@history, y@history)
  x@history <- rbind(x@history, 
                     data.frame(submitted=history.submitted,
                                finished=history.finished,
                                command=history.command))
	return(x)

}) # }}}

setMethod("corplot","MethyLumiSet",function(x,...) {  # {{{
  corvals <- cor(betas(x))
  ordering=hclust(as.dist(corvals))$order
  image(corvals[ordering,ordering])
})  # }}}
#' Normalize a MethyLumiSet, accounting for dye bias
#'
#' The Illumina GoldenGate methylation platform uses two colors, one to
#' represent the unmethylated state and the other to represent the methylated
#' state. This function corrects that dye bias and recalculates the betas based
#' on the corrected intensities.
#'
#' For HumanMethylation27 data, the function does nothing.
#'
#' For HumanMethylation450 data, the function delegates to
#' `normalizeViaControls()` the task of scaling red and green intensities
#' against a reference array (chip) which uses the closest-to-equal chip
#' (i.e. `which.min(abs(R.G.ratio - 1))`).
#'
#' The code to do this is based on code from the 'minfi' package and uses the
#' built-in red and green normalization control probes on the hm450 arrays to
#' scale the channels of the samples, so that a consistent degree of dye bias
#' is maintained for Infinium II probes across an experiment or set of
#' experiments.
#'
#' @details
#' For HumanMethylation450 data, the function delegates to
#' `normalizeViaControls()` the task of scaling red and green intensities
#' against a reference array (chip) which defaults to the first chip in a set.
#' The code to do this is based on code from the 'minfi' package and uses the
#' built-in normalization controls to scale the channels of the samples, so
#' that a consistent degree of dye bias is maintained for Infinium II probes
#' across an experiment or set of experiments. The remainder of the
#' documentation below is specific to GoldenGate data.
#'
#' The Illumina GoldenGate methylation platform uses two colors, one to
#' represent the unmethylated state and the other to represent the methylated
#' state. This function corrects that dye bias and recalculates the betas based
#' on the corrected intensities.
#'
#' As a first step, the medians for each of Cy3 and Cy5 are calculated at high
#' and low betas, representing the (nearly) fully methylated state and the
#' (nearly) fully unmethylated states. Values of Cy3 and Cy5 that are negative
#' are set to zero for this process. Then, the Cy5 medians are adjusted to
#' match those of the Cy3 channel, thereby correcting the dye bias.
#'
#' To map the new intensities back to betas, one of two map functions can be
#' used. The default is the `atan(Cy3/Cy5)`. The ratio maps using the function
#' `(Cy3/Cy3+Cy5)`. The differences should be very small, but we feel that the
#' atan map function is probably the mathematically appropriate way of doing
#' this.
#'
#' @param x A MethyLumiSet object.
#' @param beta.cuts Two numeric values with the first less than the second and
#'   between 0 and 1, representing the beta cutoffs that will be used when
#'   determining the median intensities to which to correct. See details below.
#' @param mapfun Either `"atan"` or `"ratio"`. See details below.
#' @return A new `MethyLumiSet` that contains the corrected betas and the
#'   adjusted intensities.
#' @author Sean Davis <seandavi@gmail.com>
#' @keywords manip
#' @examples
#' ## Read in sample information
#' samps <- read.table(system.file("extdata/samples.txt",
#'                                 package = "methylumi"), sep = "\t", header = TRUE)
#' ## Perform the actual data reading
#' ## This is an example of reading data from an
#' ## Sentrix Array format file (actually two files,
#' ## one for data and one for QC probes)
#' mldat <- methylumiR(system.file('extdata/exampledata.samples.txt', package = 'methylumi'),
#'                     qcfile = system.file('extdata/exampledata.controls.txt', package = "methylumi"),
#'                     sampleDescriptions = samps)
#' mldatnorm <- normalizeMethyLumiSet(mldat)
normalizeMethyLumiSet <- function(x,beta.cuts=c(0.2,0.8),mapfun=c('atan','ratio')) { # {{{

  if( length(annotation(x)) > 0 ) { 
    if( annotation(x) == 'IlluminaHumanMethylation450k' ) { 
      message('Normalizing via Illumina controls...')
      return(normalizeViaControls(x))
    }
    if( annotation(x) == 'IlluminaHumanMethylation27k' ) { 
      message('HumanMethylation27 data encountered, skipping...')
      return(x)
    }
  }

  mapfun=match.arg(mapfun)
  history.submitted <- as.character(Sys.time())
  good <- rep(TRUE,ncol(x))
  cy3 <- unmethylated(x)
  cy3[cy3<0] <- NA
  cy5 <- methylated(x)
  cy5[cy5<0] <- NA
  for(i in 1:ncol(cy5)) {
    cy3inc <- (!is.na(betas(x)[,i]) & !is.na(cy3[,i]))
    cy5inc <- (!is.na(betas(x)[,i]) & !is.na(cy5[,i]))
    cy3vec <- cy3[cy3inc,i]
    cy5vec <- cy5[cy5inc,i]
    cy3h <- median(cy3vec[betas(x)[cy3inc,i]<beta.cuts[1]])
    cy3l <- median(cy3vec[betas(x)[cy3inc,i]>beta.cuts[2]])
    cy5l <- median(cy5vec[betas(x)[cy5inc,i]<beta.cuts[1]])
    cy5h <- median(cy5vec[betas(x)[cy5inc,i]>beta.cuts[2]])
    corfactor <- (cy3h-cy3l)/(cy5h-cy5l)
    cy5[,i] <- cy5[,i]*(corfactor)
    cy5vec <- cy5[cy5inc,i]
    newcy5l <- median(cy5vec[betas(x)[cy5inc,i]<beta.cuts[1]])
    if(newcy5l<cy3l) {
      cy5[,i] <- cy5[,i]+(cy3l-newcy5l)
    } else {
      cy3[,i] <- cy3[,i]+(newcy5l-cy3l)
    }
    if(corfactor<0) {
      good[i] <- FALSE
      warning(sprintf("Sample %d has medians that do not make sense for a normal sample\n(cy3l=%f ,cy5l=%f ,cy3h=%f ,cy5h=%f)\nRemoving sample!  Check quality control.",
                      i,cy3l,cy5l,cy3h,cy5h))
      cy5[,i] <- NA
      cy3[,i] <- NA
    }
#    print(sprintf("cy5l %f cy3l %f cy5h %f cy3h %f",median(cy5[,i][betas(x)[,i]<beta.cuts[1]]),
#                  median(cy3[,i][betas(x)[,i]>beta.cuts[2]]),
#                  median(cy5[,i][betas(x)[,i]>beta.cuts[2]]),
#                  median(cy3[,i][betas(x)[,i]<beta.cuts[1]])))
  }
  newbeta <- 0
  if(mapfun=='atan') {
    newbeta <- atan((cy5)/(cy3))/(pi/2)
  } else {
    newbeta <- cy5/(cy5+cy3+100)
  }
  assaydata <- new.env(hash=TRUE,parent=emptyenv())
  assaydata[['unmethylated']] <- cy3
  assaydata[['methylated']] <- cy5
  assaydata[['betas']] <- newbeta
                                        #assaydata[['pvals']] <- assayData(x)$pvals
  history.finished <- as.character(Sys.time())
  history.command <- capture.output(print(match.call(normalizeMethyLumiSet)))  
  ret <- new("MethyLumiSet",phenoData=phenoData(x),featureData=featureData(x),
             assayData=assaydata,annotation=annotation(x))
  QCdata(ret) <- QCdata(x)
  ret <- ret[,good]
  ret@history <- rbind(getHistory(x), data.frame(submitted=history.submitted,finished=history.finished,command=history.command))
  return(ret)
} # }}}

setGeneric("parplot", function(object,...) { # {{{
           standardGeneric("parplot")
}) # }}}
.parallel <- function(object,quantiles,what,...) { # {{{
  parallelplot(apply(what(object),2,function(x,quantiles) {
    a <- ecdf(x)
    return(a(quantiles))
    },quantiles),...)
} # }}}
setMethod("parplot", signature(object="MethyLumi"), function(object,quantiles=seq(0,1,0.2),what=c("betas","exprs","pvals"),...) { # {{{
                     what=match.arg(what)
                     what=get(what)
                     .parallel(object,quantiles=quantiles,what=what,...)
                   })  #}}}

#' Methods for dealing with control data for Illumina methylation data
#'
#' The `qcplot` function simply generates a plot of the control probe
#' information for a given `controlType`.
#'
#' The descriptions of the various control types can be obtained from the
#' Illumina methylation user's guides.
#'
#' @param object An object of class [MethyLumiSet-class] or
#'   [MethyLumiQC-class].
#' @param controltype A single character value representing the bead type to
#'   plot from the quality control data. The available types are accessible via
#'   the `controlTypes` method.
#' @param ... Passed to the plot function.
#' @author Sean Davis <seandavi@gmail.com>
#' @seealso [MethyLumiSet-class], [MethyLumiQC-class]
#' @keywords hplot
#' @examples
#' data(mldat)
#' controlTypes(mldat)
#' qcplot(mldat, controlTypes(mldat)[3])
#' @name qcplot
#' @aliases qcplot controlTypes
#' @usage qcplot(object,controltype,...)
#' controlTypes(object,...)
NULL

if (is.null(getGeneric("qcplot"))) {  # {{{
  setGeneric("qcplot", function(object,controltype,...) {
      standardGeneric("qcplot")
  })
}  # }}}
.qcplot <- function(object,controltype,...) { # {{{
  rows <- grep(controltype,featureNames(object))
  arraytype <- "goldengate"
  if(length(grep("Signal_Red",assayDataElementNames(object),ignore.case=TRUE))>0)
    arraytype <- "infinium"
  ## Had to change the stuff below to use 'methylated' and
  ## 'unmethylated' since I now rename these columns everywhere
#  datElements <- switch(arraytype,
#                        infinium=c('Signal_Red','Signal_Grn'),
#                        goldengate=c('Signal CY3','Signal CY5')
#                        )
  datElements <- c('unmethylated','methylated')
  plot(dotplot(t(assayDataElement(object,datElements[1])
          [grep(controltype,featureNames(object)),]),
               xlab=datElements[1],main=controltype,auto.key=TRUE),
       split=c(1,1,2,1))
  plot(dotplot(t(assayDataElement(object,datElements[2])
          [grep(controltype,featureNames(object)),]),
               auto.key=TRUE,xlab=datElements[2],main=controltype),
       split=c(2,1,2,1),newpage=FALSE)
}  # }}}
setMethod("qcplot", signature(object="MethyLumiQC"),  #{{{
          function(object,controltype="NON",...) {
            return(.qcplot(object,controltype,...))}
          )  # }}}
setMethod("qcplot", signature(object="MethyLumiSet"), # {{{
          function(object,controltype="NON",...) {
            qcplot(QCdata(object),controltype)}
          ) # }}}
