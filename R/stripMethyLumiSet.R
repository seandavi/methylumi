## FIXME: the real answer is to farm the data out to ff or hdf5, like Benilton
#' Strip excessive probe-level data from MethyLumiSets
#'
#' 450k datasets with probe-level standard errors, out-of-band intensities and
#' bead numbers can become huge. These functions help to manage their growth in
#' memory, at least until preprocessing and QC is completed, whereupon the
#' summary data can be exported to a RangedData-based object of some sort for
#' integration.
#'
#' @param object an object of class MethyLumi or a subclass.
#' @return The object, with the corresponding assay data elements removed.
#' @author Tim Triche, Jr. <tim.triche@gmail.com>
#' @name MethyLumi-strippers
#' @aliases stripOOB stripBeadNs stripBeadSDs stripMethyLumiSet
NULL

#' @rdname MethyLumi-strippers
stripBeadNs <- function(object) { # {{{ make 450k datasets smaller
  if( 'methylated.N' %in% assayDataElementNames(object) ||
      'unmethylated.N' %in% assayDataElementNames(object) ) { 
    history.submitted <- as.character(Sys.time())
    aData <- assayData(object)
    storageMode(aData) <- "environment"
    rm('unmethylated.N', envir=aData)
    rm('methylated.N', envir=aData)
    storageMode(aData) <- "lockedEnvironment"
    assayData(object) <- aData
    history.finished <- as.character(Sys.time())
    history.command <- deparse(match.call())
    oldHistory <- object@history
    newHistory <- data.frame(submitted=history.submitted, 
                             finished=history.finished,
                             command=history.command)
    if(is(object, 'MethyLumiM')) newHistory$lumiVersion = NA
    object@history<- rbind(oldHistory, newHistory)
  }
  return(object)
} # }}}
#' @rdname MethyLumi-strippers
stripBeadSDs <- function(object) { # {{{
  if( 'methylated.SD' %in% assayDataElementNames(object) ||
      'unmethylated.SD' %in% assayDataElementNames(object) ) { 
    history.submitted <- as.character(Sys.time())
    aData <- assayData(object)
    storageMode(aData) <- "environment"
    rm('unmethylated.SD', envir=aData)
    rm('methylated.SD', envir=aData)
    storageMode(aData) <- "lockedEnvironment"
    assayData(object) <- aData
    history.finished <- as.character(Sys.time())
    history.command <- deparse(match.call())
    oldHistory <- object@history
    newHistory <- data.frame(submitted=history.submitted, 
                             finished=history.finished,
                             command=history.command)
    if(is(object, 'MethyLumiM')) newHistory$lumiVersion = NA
    object@history<- rbind(oldHistory, newHistory)
  }
  return(object)
} # }}} 
#' @rdname MethyLumi-strippers
stripOOB <- function(object) { # {{{
  if( 'methylated.OOB' %in% assayDataElementNames(object) ||
      'unmethylated.OOB' %in% assayDataElementNames(object) ) { 
    history.submitted <- as.character(Sys.time())
    aData <- assayData(object)
    storageMode(aData) <- "environment"
    rm('unmethylated.OOB', envir=aData)
    rm('methylated.OOB', envir=aData)
    storageMode(aData) <- "lockedEnvironment"
    assayData(object) <- aData
    history.finished <- as.character(Sys.time())
    history.command <- deparse(match.call())
    oldHistory <- object@history
    newHistory <- data.frame(submitted=history.submitted, 
                             finished=history.finished,
                             command=history.command)
    if(is(object, 'MethyLumiM')) newHistory$lumiVersion = NA
    object@history<- rbind(oldHistory, newHistory)
  }
  return(object)
} # }}}
#' @rdname MethyLumi-strippers
stripMethyLumiSet <- function(object) { # {{{
  return(stripOOB(stripBeadSDs(stripBeadNs(object))))
} # }}}
