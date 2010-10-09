# $HeadUrl$
# $Id$



###
### Params:
###    filename:
### Returns:
###    Either "finalreport" or "goldengate"
###
### Currently, the first line of the file is checked to
### see if it contains "[Header]" in a case-insensitive
### sense.  If so, the "finalreport" format is assumed.
### Otherwise, the "goldengate" format is assumed.
###
.getFileType <- function(filename) {
  f <- toupper(readLines(filename,n=1))
  if(f[1]=="[HEADER]") {
    return("finalreport")
  } else {
    return("goldengate")
  }
}
  

###
### Params:
###    dat: data frame containing colnames and data from a
###         methylation profile
### Returns:
###    list with two members, featuredata and assaydata
###
### Makes the assumption that sample names and data type are
### separated by either "." or ":"
###

getAssayDataNameSubstitutions <- function() {
  subs <- read.table(system.file("extdata/substitutions.txt",package="methylumi")
                     ,header=TRUE,as.is=TRUE,sep="\t")
  return(subs)
}

.doAssayDataNameSubstitutions <- function(assayDataNames) {
  assayDataNamesCopy <- assayDataNames
  subs <- getAssayDataNameSubstitutions()
  for(i in 1:nrow(subs)) {
    tmp <- grep(subs$regex[i],assayDataNames)
    if(length(tmp)>1) {
      warning(sprintf("Found greater than 1 match in assayDataNames for regex %s\nso will use the last one",
                      subs$regex[i]))
    }
    if(length(tmp)>0) {
      assayDataNamesCopy[max(tmp)] <- subs$result[i]
    }
  }
  return(assayDataNamesCopy)
}

.extractMethylationFdataAndAssayData <- function(dat) {
  datcolsep <- "[.:]"
  cn <- colnames(dat)
  datcnidx <- grep(datcolsep,cn)
  datcn <- cn[datcnidx]
  cnSplit <- do.call(rbind,strsplit(cn,datcolsep))
  datcnSplit <- do.call(rbind,strsplit(datcn,datcolsep))

  dattypes <- data.frame(original=unique(datcnSplit[,2]),
                        newnames=.doAssayDataNameSubstitutions(unique(datcnSplit[,2])), stringsAsFactors = F)  ## turn off stringAsFactors
 #                        newnames=.doAssayDataNameSubstitutions(unique(datcnSplit[,2]))) # changed by Pan Du, July 1, 2010
  
	## added by Pan Du, July 1, 2010
	## Check the number of columns of each data types, remove those do not match with the substituted ones, which should have the same dimensions.
	colnum <- table(datcnSplit[,2])
	## If the column numbers for each datatype are not consistent, we need to futher check it and remove the inconsistent ones.
	if (length(unique(colnum)) > 1) {
		substitutedNames <- dattypes$original[dattypes$original != dattypes$newnames]
		## check whether all substituted data types have the same number of columns
		if (length(unique(colnum[substitutedNames])) != 1) stop("Parsed column dimension doesnot match!\n")
		## check whether other data types also have the same number of columns as the substituted ones	
		rmInd <- which(colnum[dattypes$original] != colnum[substitutedNames[1]])
		dattypes <- dattypes[-rmInd,]
	}

  assaydata <- new.env(hash=TRUE,parent=emptyenv())
  featurenames <- make.unique(as.character(dat$TargetID))
  for (i in 1:nrow(dattypes)) {
    colsOfInterest <- grep(dattypes$original[i],cn)
    tmpmat <- as.matrix(dat[,colsOfInterest])
    colnames(tmpmat) <- cnSplit[colsOfInterest,1]
    rownames(tmpmat) <- featurenames
    assaydata[[as.character(dattypes$newnames[i])]] <- tmpmat
  }

  ## added by Pan Du, July 1, 2010
  storageMode(assaydata) <- "lockedEnvironment"

  fd <- dat[,(!(cn %in% datcn)), drop=FALSE]  # add drop = FALSE by Pan Du on Oct 6, 2010
  rownames(fd) <- featurenames
  featuredata <- as(fd,"AnnotatedDataFrame")
  return(list(featuredata=featuredata,assaydata=assaydata))
}


#########################################################
###
### The following few functions:
###    .getBlock
###    .getHeader
###    .getInfiniumMethylationProfile
###    .readInfiniumMethylationFile
###
### are for Final Report format.
###
#########################################################

###
### Params:
###    filename: filename of the Final Report file
###    blocks:  data.frame with block name, start, and end in line
###       numbers
###    blockname: The specific name of the block to be retrieved
###    required: boolean for if the block is required or not
### Returns:
###    A data frame for the data between the block start and end
###
.getFinalReportBlock <- function(filename,blockname,blocks,required,header=FALSE,...) {
  if (!(blockname %in% blocks$block) & required) {
    stop(sprintf('The block called "%s" in the file "%s" is not present and needs to be.  Please re-export from beadstudio',blockname,filename))
  }
  if (!(blockname %in% blocks$block)) {
    return(NULL)
  }
  skip=blocks[blocks$block==blockname,'start']
  nrows=blocks[blocks$block==blockname,'nrows']
  if(header) nrows=nrows-1
  dat <- read.table(filename,skip=skip,nrows=nrows,sep="\t",
                    header=header,comment.char="",quote="",fill=TRUE,...)
  return(dat)
}


###
### Params:
###    filename: filename of the Final Report file
###    blocks:  data.frame with block name, start, and end in line
###       numbers
### Returns:
###    A data frame with two columns, Tag and Value
###
.getFinalReportHeader <- function(filename,blocks,...) {
  tmp <- .getFinalReportBlock(filename,blockname="[HEADER]",blocks,required=TRUE,header=FALSE)
  colnames(tmp) <- c("Tag","Value")
  return(tmp)
}

.getFinalReportMethylationProfile <- function(filename,blocks,blockname,required,...) {
  dat <- .getFinalReportBlock(filename,blockname=blockname,blocks=blocks,
                              required=required,header=TRUE,check.names=FALSE)
  if(is.null(dat)) {
    return(NULL)
  }
  return(.extractMethylationFdataAndAssayData(dat))
}


.readFinalReportMethylationFile <- function(filename,...) {
  tmpdat <- readLines(filename)
  tmp <- grep('\\[.*\\]$',tmpdat)
  blocks <- data.frame(start=tmp,nrows=c(tmp[2:length(tmp)],length(tmpdat)+1)-tmp-1,block=toupper(tmpdat[tmp]))
  header <- .getFinalReportHeader(filename,blocks)
  datblock <- .getFinalReportMethylationProfile(filename,blocks=blocks,blockname="[SAMPLE METHYLATION PROFILE]",required=TRUE,header=TRUE,fill=FALSE)
  qcblock <- .getFinalReportMethylationProfile(filename,blocks=blocks,blockname="[CONTROL PROBE PROFILE]",required=FALSE,header=TRUE)
  return(list(header=header,datblock=datblock,qcblock=qcblock))
}


###
### Returns:
###    list of datblock, qcblock (set to NULL if qcfile
###    was NULL, and header (always NULL for the old
###    format.
###
### This function reads the old format methylation file
### that does not have data blocks but is simply tab-
### delimited text.
###
.readOldMethylationFile <- function(filename,qcfile=NULL,...) {
  dat <- read.delim(filename,check.names=FALSE,...)
  datblock <- .extractMethylationFdataAndAssayData(dat,...)
  qcblock <- NULL
  if(!(is.null(qcfile))) {
    dat <- read.delim(qcfile,check.names=FALSE,...)
    qcblock <- .extractMethylationFdataAndAssayData(dat,...)
  }
  return(list(datblock=datblock,qcblock=qcblock,header=NULL))
}

###
### User function
###
### Autodetects file format, reads data, makes substititutions
### of assayData element names for the various formats so that
### methods accessors will work, maps sampleDescriptions
### appropriately, and returns a MethyLumiSet object
###
methylumiR <-
  function(filename,qcfile=NULL,sampleDescriptions=NULL,...) {

    history.submitted <- as.character(Sys.time())
    ## Determine the file type from the file itself
    ftype <- .getFileType(filename)
    fulldat <- NULL
    ## Parse the files
    ## and return a list with three members:
    ##   datblock: contains actual data
    ##   qcblock : contains QC data, if available, NULL otherwise
    ##   header: present only for Final Report type data
    if(ftype=="finalreport") {
      if(!(is.null(qcfile))) {
        warning("qcfile specification is not supported for Final Report type BeadStudio Export files")
      }
      fulldat <- .readFinalReportMethylationFile(filename)
    } else {
      fulldat <- .readOldMethylationFile(filename,qcfile)
    }
    
    sampleName <-  sampleNames(fulldat$datblock$assaydata)
	
	## added by Pan Du, July 1, 2010
	if (any(duplicated(sampleName))) {
		warning("Duplicated column names found!\n Suffix indexes were appended!\n")
		sampleName <- make.names(sampleName, unique=T)
	}
	
    sampleID <- sampleName
    label <- NULL
    pData <- NULL

    if(!is.null(sampleDescriptions)) {
      sampleIDcol <- grep('SampleID',colnames(sampleDescriptions),ignore.case=TRUE)
      if(!sampleIDcol) {
        stop('If sampleDescriptions is supplied, one column must be named SampleID')
      }
      sampleMatch <- match(sampleName,as.character(sampleDescriptions[,sampleIDcol]))
      if(any(is.na(sampleMatch))) {
        stop('Nonmatching sampleID(s) in Beadstudio file and sampleDescriptions')
      }
      sampleDescriptions <- sampleDescriptions[sampleMatch,]
      pData <- sampleDescriptions
      colnames(pData)[sampleIDcol] <- 'sampleID'
      labelCol <- grep('SampleLabel',colnames(sampleDescriptions),ignore.case=TRUE)
      if(labelCol) {
        label <- sampleDescriptions[,labelCol]
      }
    }
    if(is.null(label)) label <- sampleName
    if (length(unique(label)) == length(label) & length(label) > 0) {
      colName <- label
    } else {
      colName <- sampleName
    }
    
    sampleNames(fulldat$datblock$assaydata) <- colName
    if(!(is.null(fulldat$qcblock))) {
      sampleNames(fulldat$qcblock$assaydata) <- colName
    }


    if(is.null(pData)) {
      pData <- data.frame(sampleID=sampleID, label=label)
    }
    rownames(pData) <- colName
    ##pdata <- new("phenoData", pData=pData, varLabels=list('sampleID', 'label'))
    varMetadata <- data.frame(labelDescription=colnames(pData))
    rownames(varMetadata) <- c(colnames(pData))
    pdata <- new("AnnotatedDataFrame", data=pData, varMetadata=varMetadata)
    
    x.lumi=new('MethyLumiSet',assayData=fulldat$datblock$assaydata,
      featureData=fulldat$datblock$featuredata,
      phenoData=pdata)
    if(!is.null(fulldat$qcblock)) {
      control.lumi <- new('MethyLumiQC',assayData=fulldat$qcblock$assaydata,
                          featureData=fulldat$qcblock$featuredata)
      x.lumi@QC <- control.lumi
    }
    history.finished <- as.character(Sys.time())

    history.command <- paste(capture.output(print(match.call(methylumiR))),collapse="")
    
    x.lumi@history<- rbind(x.lumi@history,
                           data.frame(submitted=history.submitted,
                                      finished=history.finished,
                                      command=history.command))
    return(x.lumi)
  }

extractBarcodeAndPosition <- function(sentrixids) {
  l <- do.call(rbind,strsplit(sentrixids,'_'))
  l <- data.frame(l)
  if(ncol(l)!=3) stop("sentrix ids should contain two '_' characters")
  colnames(l) <- c('sentrix','row','column')
  l$rowNumber <- as.numeric(gsub('R','',l$row))
  l$columnNumber <- as.numeric(gsub('C','',l$column))
  return(data.frame(l))
}


