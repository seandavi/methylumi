## FIXME: the real answer is to farm the data out to ff or hdf5, like Benilton
stripBeadNs <- function(from) { # {{{ make 450k datasets smaller
  if( 'methylated.N' %in% assayDataElementNames(from) ||
      'unmethylated.N' %in% assayDataElementNames(from) ) { 
    history.submitted <- as.character(Sys.time())
    aData <- assayData(from)
    storageMode(aData) <- "environment"
    rm('unmethylated.N', envir=aData)
    rm('methylated.N', envir=aData)
    storageMode(aData) <- "lockedEnvironment"
    assayData(from) <- aData
    history.finished <- as.character(Sys.time())
    history.command <- 'Removed probe-level bead numbers from assayData.'
    oldHistory <- from@history
    newHistory <- data.frame(submitted=history.submitted, 
                             finished=history.finished,
                             command=history.command)
    if(is(from, 'MethyLumiM')) newHistory$lumiVersion = NA
    from@history<- rbind(oldHistory, newHistory)
  }
  return(from)
} # }}}
stripBeadSDs <- function(from) { # {{{
  if( 'methylated.SD' %in% assayDataElementNames(from) ||
      'unmethylated.SD' %in% assayDataElementNames(from) ) { 
    history.submitted <- as.character(Sys.time())
    aData <- assayData(from)
    storageMode(aData) <- "environment"
    rm('unmethylated.SD', envir=aData)
    rm('methylated.SD', envir=aData)
    storageMode(aData) <- "lockedEnvironment"
    assayData(from) <- aData
    history.finished <- as.character(Sys.time())
    history.command <- 'Removed probe-level stderrs from assayData.'
    oldHistory <- from@history
    newHistory <- data.frame(submitted=history.submitted, 
                             finished=history.finished,
                             command=history.command)
    if(is(from, 'MethyLumiM')) newHistory$lumiVersion = NA
    from@history<- rbind(oldHistory, newHistory)
  }
  return(from)
} # }}} 
stripOOB <- function(from) { # {{{
  if( 'methylated.OOB' %in% assayDataElementNames(from) ||
      'unmethylated.OOB' %in% assayDataElementNames(from) ) { 
    history.submitted <- as.character(Sys.time())
    aData <- assayData(from)
    storageMode(aData) <- "environment"
    rm('unmethylated.OOB', envir=aData)
    rm('methylated.OOB', envir=aData)
    storageMode(aData) <- "lockedEnvironment"
    assayData(from) <- aData
    history.finished <- as.character(Sys.time())
    history.command <- 'Removed out-of-band intensities from assayData.'
    oldHistory <- from@history
    newHistory <- data.frame(submitted=history.submitted, 
                             finished=history.finished,
                             command=history.command)
    if(is(from, 'MethyLumiM')) newHistory$lumiVersion = NA
    from@history<- rbind(oldHistory, newHistory)
  }
  return(from)
} # }}}
stripMethyLumiSet <- function(from) { # {{{
  return(stripOOB(stripBeadSDs(stripBeadNs(from))))
} # }}}
