## ---------------------------------------------------------------
## define a new class lumiMethyM
setClass('MethyLumiM', 
	representation(controlData='QCDataOrNULL', history='data.frame'), 
	prototype=list(controlData = NULL, history=data.frame(
		submitted   = I(vector()),
		finished    = I(vector()),
		command     = I(vector()),
		lumiVersion = I(vector())
	)), 
	contains='ExpressionSet')

#
setMethod('initialize', 'MethyLumiM', function(.Object, 
	exprs = new('matrix'),
	methylated = new('matrix'),		
	unmethylated = new('matrix'),
	detection = new('matrix'),  # detection pvalues
	methylated.N = new('matrix'),  # number of M beads
	unmethylated.N = new('matrix'),  # number of U beads
	controlData = NULL,
    ...,
    assayData)
{
	if (missing(assayData)) {
		cmd <- 'assayData <- assayDataNew(exprs=exprs, methylated=methylated, unmethylated=unmethylated'
		nSample <- ncol(exprs)

		if (ncol(detection) == nSample) {
		  cmd <- paste(cmd, ', detection=detection')
		}
		if (ncol(methylated.N) == nSample) {
		  cmd <- paste(cmd, ', methylated.N=methylated.N')
		}
		if (ncol(unmethylated.N) == nSample) {
		  cmd <- paste(cmd, ', unmethylated.N=unmethylated.N')
		}

		cmd <- paste(cmd, ')')
		eval(parse(text=cmd))
	} else if (!missing(exprs)) 
		stop("only one of 'assayData' or ('exprs', 'methylated' and 'unmethylated') allowed")

	callNextMethod(.Object, assayData=assayData, ...)
})


setValidity("MethyLumiM", function(object) 
{
    msg <- Biobase:::validMsg(NULL, Biobase:::isValidVersion(object, "ExpressionSet"))
    msg <- Biobase:::validMsg(msg, assayDataValidMembers(assayData(object), c("exprs")))
    msg <- Biobase:::validMsg(msg, assayDataValidMembers(assayData(object), c("methylated")))
    msg <- Biobase:::validMsg(msg, assayDataValidMembers(assayData(object), c("unmethylated")))
    if (is.null(msg)) TRUE else msg
})


setAs("eSet", "MethyLumiM", function(from) {

  	history.submitted <- as.character(Sys.time())

	# from <- asS4(from)
	
	if (exists('methylated', assayData(from)) && exists('methylated', assayData(from))) {
		M <- estimateM(from, returnType="matrix")
	} else {
		stop("Cannot convert as MethyLumiM object because methylated and unmethylated slots do not exist!\n")
	}
	
	if (.hasSlot(from, "history")) {
		history <- from@history
		if (is.null(history$lumiVersion)) history$lumiVersion <- rep(NA, nrow(history))
	} else {
		history <- data.frame(submitted = I(vector()), finished = I(vector()), command = I(vector()), lumiVersion = I(vector()))
	}
	
	aData <- assayData(from)
	detection <- NULL
	if (exists('pvals', assayData(from))) {
		detection <- aData[['pvals']]
		storageMode(aData) <- "environment"
		aData[['detection']] <- detection
		storageMode(aData) <- "lockedEnvironment"
	} else if (exists('detection', assayData(from))) {
		detection <- aData[['detection']]
	}
	
	methy <- aData[['methylated']]
	unmethy <- aData[['unmethylated']]
	methy.N <- unmethy.N <- NULL
	if (is.element('methylated.N', assayDataElementNames(from)) && is.element('unmethylated.N', assayDataElementNames(from))) {
		methy.N <- aData[['methylated.N']]
		dimnames(methy.N) <- dimnames(aData[['methylated']])
		unmethy.N <- aData[['unmethylated.N']]
		dimnames(unmethy.N) <- dimnames(aData[['unmethylated']])
	}
	storageMode(aData) <- "environment"
	ts <- ls(envir=aData)
	rm(list=ts, envir=aData)
	aData[['exprs']] <- M
	aData[['methylated']] <- methy
	aData[['unmethylated']] <- unmethy
	if (!is.null(detection)) aData[['detection']] <- detection
	if (!is.null(methy.N)) aData[['methylated.N']] <- methy.N
	if (!is.null(unmethy.N)) aData[['unmethylated.N']] <- unmethy.N
	storageMode(aData) <- "lockedEnvironment"	
	to <- new("MethyLumiM", assayData=aData, phenoData=phenoData(from), featureData=featureData(from), annotation=annotation(from), experimentData=experimentData(from), protocolData=protocolData(from))		

	# check whether there are QC data available
	if (.hasSlot(from, "QC")) {
		if (is(QCdata(from), "MethyLumiQC")) {
			to@controlData <- from@QC
		} 
	}

	history.finished <- as.character(Sys.time())
	history.command <- capture.output(print(match.call(setAs)))  
	lumiVersion <- packageDescription('lumi')$Version
	to@history <- rbind(history, 
                      data.frame(submitted=history.submitted, 
                                 finished=history.finished, 
                                 command=history.command, 
                                 lumiVersion=lumiVersion))
	return(to)
})


setMethod("methylated", signature(object="MethyLumiM"), function(object) {
	if ('methylated' %in% assayDataElementNames(object)) {
		return(assayDataElement(object,"methylated"))
	} else {
		return(NULL)
	}
})

setReplaceMethod("methylated", signature(object="MethyLumiM"), function(object, value) {
	if (is.null(value)) {
		assay <- assayData(object)
		if (exists('methylated', envir=assay)) {
			oldMode <- storageMode(assay)
			storageMode(assay) <- 'environment'
			rm(methylated, envir=assay)
			storageMode(assay) <- oldMode
			assayData(object) <- assay
		}
		return(object)
	} else {
		assayDataElementReplace(object, "methylated", value)
	}
})	


setMethod("unmethylated", signature(object="MethyLumiM"), function(object) {
	if ('unmethylated' %in% assayDataElementNames(object)) {
		return(assayDataElement(object,"unmethylated"))
	} else {
		return(NULL)
	}
})

setReplaceMethod("unmethylated", signature(object="MethyLumiM"), function(object, value) {
	if (is.null(value)) {
		assay <- assayData(object)
		if (exists('methylated', envir=assay)) {
			oldMode <- storageMode(assay)
			storageMode(assay) <- 'environment'
			rm(unmethylated, envir=assay)
			storageMode(assay) <- oldMode
			assayData(object) <- assay
		}
		return(object)
	} else {
		assayDataElementReplace(object, "unmethylated", value)
	}
})	
	


setMethod("detection", signature(object="MethyLumiM"), function(object) {
	if ('detection' %in% assayDataElementNames(object)) {
		return(assayDataElement(object,"detection"))
	} else {
		return(NULL)
	}
})

setReplaceMethod("detection", signature(object="MethyLumiM"), function(object, value) {
	if (is.null(value)) {
		assay <- assayData(object)
		if (exists('detection', envir=assay)) {
			oldMode <- storageMode(assay)
			storageMode(assay) <- 'environment'
			rm(detection, envir=assay)
			storageMode(assay) <- oldMode
			assayData(object) <- assay
		}
		return(object)
	} else {
		assayDataElementReplace(object, "detection", value)
	}
})	


setMethod("methylated.N", signature(object="MethyLumiM"), function(object) {
	if ('methylated.N' %in% assayDataElementNames(object)) {
		return(assayDataElement(object,"methylated.N"))
	} else {
		return(NULL)
	}
})

setReplaceMethod("methylated.N", signature(object="MethyLumiM"), function(object, value) {
	if (is.null(value)) {
		assay <- assayData(object)
		if (exists('methylated.N', envir=assay)) {
			oldMode <- storageMode(assay)
			storageMode(assay) <- 'environment'
			rm(methylated.N, envir=assay)
			storageMode(assay) <- oldMode
			assayData(object) <- assay
		}
		return(object)
	} else {
		assayDataElementReplace(object, "methylated.N", value)
	}
})	


setMethod("unmethylated.N", signature(object="MethyLumiM"), function(object) {
	if ('unmethylated.N' %in% assayDataElementNames(object)) {
		return(assayDataElement(object,"unmethylated.N"))
	} else {
		return(NULL)
	}
})

setReplaceMethod("unmethylated.N", signature(object="MethyLumiM"), function(object, value) {
	if (is.null(value)) {
		assay <- assayData(object)
		if (exists('unmethylated.N', envir=assay)) {
			oldMode <- storageMode(assay)
			storageMode(assay) <- 'environment'
			rm(unmethylated.N, envir=assay)
			storageMode(assay) <- oldMode
			assayData(object) <- assay
		}
		return(object)
	} else {
		assayDataElementReplace(object, "unmethylated.N", value)
	}
})	


setMethod("controlData", signature(object="MethyLumiM"), function(object) {
	if (.hasSlot(object, "controlData")) {
		return(object@controlData)
	} else {
		return(NULL)
	}
})


setReplaceMethod("controlData", signature(object="MethyLumiM"), function(object, value) {
	if (is(value, "MethyLumiQC")) {
		object@controlData <- value
	} else {
		cat("The control data should be a MethyLumiQC object!\n")
	}
	return(object)
})	


setMethod("getHistory",signature(object="MethyLumiM"), function(object) object@history)


setMethod("combine", signature=c(x="MethyLumiM", y="MethyLumiM"), function(x, y, ...) {

   	history.submitted <- as.character(Sys.time())

	if (missing(y)) return(x)
	if (length(list(...)) > 0) 
	        return(combine(x, combine(y, ...)))

  	## do default processing of 'ExpressionSet'
  	x.comb <- callNextMethod()

	## deal with control data
	if (!is.null(controlData(x)) && !is.null(controlData(y))) {
		controlData(x.comb) <- combine(controlData(x), controlData(y))
	}
	
    # history tracking
    history.finished <- as.character(Sys.time())
	#history.command <- match.call()
    history.command <- capture.output(print(match.call(combine)))  
	x.comb@history<- rbind(x@history, y@history)
	if (is.null(x.comb@history$lumiVersion) && nrow(x@history) > 0) {
		x.comb@history <- data.frame(x.comb@history, lumiVersion=rep(NA, nrow(x.comb@history)))
	} 
	lumiVersion <- packageDescription('lumi')$Version
	x.comb@history<- rbind(x.comb@history, data.frame(submitted=history.submitted,finished=history.finished,command=history.command, lumiVersion=lumiVersion))
	return(x.comb)
})


setMethod("[", "MethyLumiM", function(x, i, j, ..., drop = FALSE)  {

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
	  	if (!is.null(controlData(x)))
	  		controlData(x) <- controlData(x)[,j]
	} else {
	  	return(x)
	}
	
	# history tracking
	history.finished <- as.character(Sys.time())
	if (is.null(x@history$lumiVersion) && nrow(x@history) > 0) {
		x@history <- data.frame(x@history, lumiVersion=rep(NA, nrow(x@history)))
	}
	lumiVersion <- packageDescription('lumi')$Version
	x@history<- rbind(x@history, data.frame(submitted=history.submitted,finished=history.finished, command=history.command, lumiVersion=lumiVersion))
	
	return(x)
})



## ------------------------------------------------------------------------------
## Other functions designed for MethyLumiM class object

# estimate the M-value based on methylated and unmethylated probe intensities
estimateM <- function(methyLumiM, returnType=c("ExpressionSet", "matrix"), offset=1) {
	
	if (!assayDataValidMembers(assayData(methyLumiM), c("unmethylated", "methylated"))) {
		stop("The input should include 'methylated' and 'unmethylated' elements in the assayData slot!\n")
	}
	returnType <- match.arg(returnType)
	
	unmethy <- assayDataElement(methyLumiM, 'unmethylated') 
	methy <- assayDataElement(methyLumiM, 'methylated') 
	mm <- min(c(unmethy, methy))
	if (mm < 0.01) {
		unmethy[unmethy < 0.01] <- 0.01 
		methy[methy < 0.01] <- 0.01 
	}
	M <- log2((methy + offset) / (unmethy + offset))
	if (returnType == "matrix") {
		return(M)
	} else {
		exprs(methyLumiM) <- M
		return(methyLumiM)
	}
}


