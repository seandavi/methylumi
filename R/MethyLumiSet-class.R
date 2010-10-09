# $HeadUrl$
# $Id$
##=================================================
## Define MethyLumiSet object:
setClass('MethyLumi',
         contains='eSet')

setClass('MethyLumiQC',
         contains='MethyLumi')

setClassUnion("QCDataOrNULL",c('NULL',"MethyLumiQC"))

setClass('MethyLumiSet',
         representation(QC="QCDataOrNULL",
                        history='data.frame'),
         prototype=list(QC=NULL,
           history=data.frame(
             submitted   = I(vector()),
             finished    = I(vector()),
             command     = I(vector())
             )),
         contains='MethyLumi')


setMethod('initialize', 'MethyLumiSet',
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
          )

setMethod('initialize', 'MethyLumiQC',
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
          )


setValidity("MethyLumiSet", function(object) 
{
    msg <- Biobase:::validMsg(NULL, Biobase:::isValidVersion(object, "eSet"))
    msg <- Biobase:::validMsg(msg, assayDataValidMembers(assayData(object), c("betas")))
    if (is.null(msg)) TRUE else msg
})

setValidity("MethyLumiQC", function(object) 
{
    msg <- Biobase:::validMsg(NULL, Biobase:::isValidVersion(object, "eSet"))
#    msg <- Biobase:::validMsg(msg, assayDataValidMembers(assayData(object), c("avgsignal")))
    if (is.null(msg)) TRUE else msg
})


##=================================================
## methods
if (is.null(getGeneric("getHistory"))) setGeneric("getHistory", function(object) standardGeneric("getHistory"))
if (is.null(getGeneric("betas"))) setGeneric("betas", function(object) standardGeneric("betas"))
if (is.null(getGeneric("betas<-"))) setGeneric("betas<-", function(object, value) standardGeneric("betas<-"))
if (is.null(getGeneric("pvals"))) setGeneric("pvals", function(object) standardGeneric("pvals"))
if (is.null(getGeneric("pvals<-"))) setGeneric("pvals<-", function(object, value) standardGeneric("pvals<-"))
if (is.null(getGeneric("exprs"))) setGeneric("exprs", function(object) standardGeneric("exprs"))
if (is.null(getGeneric("exprs<-"))) setGeneric("exprs<-", function(object, value) standardGeneric("exprs<-"))
if (is.null(getGeneric("unmethylated"))) setGeneric("unmethylated", function(object) standardGeneric("unmethylated"))
if (is.null(getGeneric("unmethylated<-"))) setGeneric("unmethylated<-", function(object, value) standardGeneric("unmethylated<-"))
if (is.null(getGeneric("methylated"))) setGeneric("methylated", function(object) standardGeneric("methylated"))
if (is.null(getGeneric("methylated<-"))) setGeneric("methylated<-", function(object, value) standardGeneric("methylated<-"))
if (is.null(getGeneric("hist"))) setGeneric("hist", function(x,...) standardGeneric("hist"))
if (is.null(getGeneric("corplot"))) setGeneric("corplot", function(x,...) standardGeneric("corplot"))
if (is.null(getGeneric("plotDensity"))) setGeneric("plotDensity", function(x,...) standardGeneric("plotDensity"))
if (is.null(getGeneric("QCdata"))) setGeneric("QCdata", function(object) standardGeneric("QCdata"))
if (is.null(getGeneric("QCdata<-"))) setGeneric("QCdata<-", function(object, value) standardGeneric("QCdata<-"))


setMethod("betas", signature(object="MethyLumiSet"),
          function(object) assayDataElement(object,"betas"))

setReplaceMethod("betas", signature(object="MethyLumiSet",value="matrix"),
                 function(object, value) assayDataElementReplace(object, "betas", value))

setMethod("pvals", signature(object="MethyLumi"),
          function(object) assayDataElement(object,"pvals"))

setReplaceMethod("pvals", signature(object="MethyLumi",value="matrix"),
                 function(object, value) assayDataElementReplace(object, "pvals", value))

setMethod("QCdata", signature(object="MethyLumiSet"),
          function(object) object@QC)

setReplaceMethod("QCdata", signature(object="MethyLumiSet",value="MethyLumiQC"),
                 function(object, value) {object@QC <- value
                                          return(object)})

setMethod("exprs", signature(object="MethyLumiSet"),
          function(object) betas(object))

setReplaceMethod("exprs", signature(object="MethyLumiSet",value="matrix"),
                 function(object, value) assayDataElementReplace(object, "betas", value))

setMethod("getHistory",signature(object="MethyLumiSet"), function(object) object@history)

if (is.null(getGeneric("summary"))) setGeneric("summary", function(object,...) standardGeneric("summary"))

setMethod("unmethylated", signature(object="MethyLumiSet"),
          function(object) {
            return(assayDataElement(object,"unmethylated"))
          })
setReplaceMethod("unmethylated", signature(object="MethyLumiSet",value="matrix"),
          function(object,value) {
            assayDataElementReplace(object,"unmethylated",value)
          })
setMethod("methylated", signature(object="MethyLumiSet"),
          function(object) {
            return(assayDataElement(object,"methylated"))
            })
setReplaceMethod("methylated", signature(object="MethyLumiSet",value="matrix"),
          function(object,value) {
            assayDataElementReplace(object,"methylated",value)
          })


setMethod("show",signature(object="MethyLumiSet"), function(object) 
{
	cat('\nObject Information:\n')
	callNextMethod()
	cat('Major Operation History:\n')
	print(getHistory(object)) 
})

setMethod("summary",signature(object="MethyLumi"), function(object) 
{
	show(object)
})


##geneNames method
if (is.null(getGeneric("combine")))
  	setGeneric("combine", function(x, y, ...)
		standardGeneric("combine"))

setMethod("combine", signature=c(x="MethyLumiSet", y="MethyLumiSet"), function(x, y) 
{
    if (class(x) != class(y))
      stop(paste("objects must be the same class, but are ",
                 class(x), ", ", class(y), sep=""))
	
	if (any(sort(featureNames(x)) != sort(featureNames(y)))) stop('Two data sets have different row names!')

   	history.submitted <- as.character(Sys.time())

    assayData(x) <- combine(assayData(x), assayData(y))
    # phenoData(x) <- combine(phenoData(x), phenoData(y))
    # featureData(x) <- combine(featureData(x), featureData(y))
    experimentData(x) <- combine(experimentData(x),experimentData(y))
	
	## combine pheno data
	if (!is.null(phenoData(x)) | !is.null(phenoData(y))) {
		phenoData.x <- phenoData(x)
		phenoData.y <- phenoData(y)
		
		pData(phenoData.x) <- merge(pData(phenoData.x), pData(phenoData.y), all=TRUE)

		metaInfo <- rbind(varMetadata(phenoData.x), varMetadata(phenoData.y))
		varMetadata(phenoData.x) <- metaInfo[!duplicated(c(rownames(varMetadata(phenoData.x)),
		 		rownames(varMetadata(phenoData.y)))), ,drop=FALSE]
		phenoData(x) <- phenoData.x
	}	
	
	## combine feature data
	if (!is.null(featureData(x)) | !is.null(featureData(y))) {
		feature.x <- featureData(x)
		feature.y <- featureData(y)
		
		repInfo <- merge(pData(feature.x), pData(feature.y), by='targetID', all=TRUE, suffixes = c(".x",".y"))
		if ('presentCount' %in% intersect(colnames(pData(feature.x)), colnames(pData(feature.y)))) {
			colInd <- which(colnames(repInfo) %in% c('presentCount.x', 'presentCount.y'))
			presentCount <- rowSums(repInfo[, colInd])
			repInfo <- repInfo[, -colInd, drop=FALSE]
			repInfo <- data.frame(repInfo, presentCount=presentCount)
		}
		pData(feature.x) <- repInfo
		metaInfo <- rbind(varMetadata(feature.x), varMetadata(feature.y))
		varMetadata(feature.x) <- metaInfo[!duplicated(c(rownames(varMetadata(feature.x)), rownames(varMetadata(feature.y)))), ,drop=FALSE]
		featureData(x) <- feature.x
	}

    # history tracking
    history.finished <- as.character(Sys.time())
	#history.command <- match.call()
    history.command <- capture.output(print(match.call(combine)))  
	x@history<- rbind(x@history, y@history)
  x@history<- rbind(x@history, data.frame(submitted=history.submitted,finished=history.finished,command=history.command))
	return(x)
})



if(!isGeneric('boxplot')) {
  setGeneric('boxplot',function(x,...) standardGeneric('boxplot'))
}

##some special handling of main is needed
setMethod("boxplot",signature(x="MethyLumiSet"),
function(x, range=0, main, logMode=TRUE, ...)
{
  tmp <- description(x)
  if (missing(main) && (is(tmp, "MIAME")))
    main <- tmp@title
  exprs <- exprs(x)
  if (nrow(x) > 5000) {
	  	index <- seq(1, nrow(x), len=5000)
              } else {
		index <- 1:nrow(x)
              }
  if (logMode & max(exprs(x), na.rm=TRUE) > 50) {
    exprs <- log2(exprs)
  } 
  dataMatrix <- exprs[index,]
  labels <- colnames(dataMatrix)
  if (is.null(labels)) labels <- as.character(1:ncol(dataMatrix))
  ## set the margin of the plot
  mar <- c(max(nchar(labels))/2 + 4.5, 5, 5, 3)
  old.mar <- par('mar')
  old.xaxt <- par('xaxt')
  par(xaxt='n')
  par(mar=mar)
  boxplot(dataMatrix ~ col(dataMatrix), main=main, range=range, xlab='', ylab='amplitude', ...)
  par(xaxt='s')
  axis(1, at=1:ncol(dataMatrix), labels=labels, tick=TRUE, las=2)
  par(mar=old.mar)
  par(xaxt=old.xaxt)
})

if (is.null(getGeneric("pairs"))) setGeneric("pairs", function(x,...) standardGeneric("pairs"))

setMethod("pairs", signature(x="MethyLumiSet"), 
	function(x,...,logMode=FALSE,maxpairs=5,fold=0.1) 
{
    
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
})

if (is.null(getGeneric("plotSampleIntensities"))) setGeneric("plotSampleIntensities", function(x,beta.cuts=c(0.2,0.8),s=1) standardGeneric("plotSampleIntensities"))

myDensity <- function(x) {
x <- x[!is.na(x)]
return(density(x))
}

setMethod('plotSampleIntensities', signature(x='MethyLumiSet'), 
function(x, beta.cuts=c(0.2,0.8),s=1) {
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
})


setMethod("hist",signature(x="MethyLumiSet"),function(x,...) {
  hist(betas(x),xlab="Beta",...)
})

setMethod("[", "MethyLumiSet", function(x, i, j, ..., drop = FALSE) 
{
  if (missing(drop)) drop <- FALSE
  history.submitted <- as.character(Sys.time())
  
  ## do default processing of 'ExpressionSet'
  x <- callNextMethod()
  
  ddim <- dim(x)
  if (!missing(i) & !missing(j)) {
    history.command <- paste('Subsetting', ddim[1], 'features and', ddim[2], 'samples.')		
  } else if (!missing(i)) {
    history.command <- paste('Subsetting', ddim[1], 'features.')
  } else if (!missing(j)) {
    history.command <- paste('Subsetting', ddim[2], 'samples.')
    QCdata(x) <- QCdata(x)[,j]
  } else {
    return(x)
  }
  
                                        # history tracking
  history.finished <- as.character(Sys.time())
  x@history<- rbind(x@history, data.frame(submitted=history.submitted,finished=history.finished,command=history.command))
  
  return(x)
})

setMethod("corplot","MethyLumiSet",function(x,...) {
  corvals <- cor(betas(x))
  ordering=hclust(as.dist(corvals))$order
  image(corvals[ordering,ordering])
})


normalizeMethyLumiSet <- function(x,beta.cuts=c(0.2,0.8),mapfun=c('atan','ratio'))
{
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
}


#####
###
### Graphics
###
#####

### Parallel plot of ecdf
###    This plots the amount of data at each cut given in quantiles
###
if (is.null(getGeneric("parplot"))) setGeneric("parplot", function(object,...) standardGeneric("parplot"))

.parallel <- function(object,quantiles,what,...) {
  parallel(apply(what(object),2,function(x,quantiles) {
    a <- ecdf(x)
    return(a(quantiles))
    },quantiles),...)
}


setMethod("parplot", signature(object="MethyLumi"),
          function(object,
                   quantiles=seq(0,1,0.2),
                   what=c("betas","exprs","pvals"),...)
                   {
                     what=match.arg(what)
                     what=get(what)
                     .parallel(object,quantiles=quantiles,what=what,...)
                   })

if (is.null(getGeneric("qcplot"))) setGeneric("qcplot", function(object,controltype,...) standardGeneric("qcplot"))


.qcplot <- function(object,controltype,...) {
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
}

setMethod("qcplot", signature(object="MethyLumiQC"),
          function(object,controltype="NON",...) {
            return(.qcplot(object,controltype,...))}
          )

setMethod("qcplot", signature(object="MethyLumiSet"),
          function(object,controltype="NON",...) {
            qcplot(QCdata(object),controltype)}
          )

if (is.null(getGeneric("controlTypes"))) setGeneric("controlTypes", function(object,...) standardGeneric("controlTypes"))

setMethod("controlTypes",signature(object="MethyLumiQC"),
          function(object) {
            return(unique(sapply(strsplit(featureNames(object),'\\.'),function(x) x[1])))})

setMethod("controlTypes",signature(object="MethyLumiSet"),
          function(object,...) {
            return(controlTypes(QCdata(object)))})

