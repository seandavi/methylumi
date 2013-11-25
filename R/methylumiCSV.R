## backwards compatibility with methylumiCSV 
##
CSVtoDF <- function(x, parallel=F, chans=c(Cy3='GRN',Cy5='RED')) { # {{{ tidy up
  if(length(x) > 1) {
    names(x) <- x
    DFs <- lapply(x, CSVtoDF)
    names(DFs) <- x
    return(DFs)
  } else {
    df <- read.csv(paste(x, 'csv', sep='.'))
    for(chan in names(chans)) {
      names(df) <- gsub(paste('Mean',chans[[chan]],sep='.'),
                        paste(chan,'Mean',sep='.'), 
                        names(df))
      names(df) <- gsub(paste('Dev',chans[[chan]],sep='.'),
                        paste(chan,'SD',sep='.'), 
                        names(df))
      df[[ paste(chan, 'NBeads', sep='.') ]] <- df$N
    }
    rownames(df) <- df$Illumicode
    return(df)
  }
} # }}}

methylumiCSV <- function(barcodes,pdat=NULL,n=T,n.sd=F,oob=T,parallel=F, ...){ # {{{

    if(any(duplicated(barcodes))|any(grepl('csv',ignore.case=T,barcodes))) {
      message('Warning: filtering out duplicate barcodes and raw filenames')
      barcodes = unique(gsub('_(Red|Grn).idat','',barcodes, ignore.case=TRUE))
    }
    n = length(barcodes)
    nch = DFsToNChannelSet(CSVtoDF(barcodes), IDAT=FALSE, parallel=parallel)
    annotation(nch) = 'IlluminaHumanMethylation27k'
    mlumi = NChannelSetToMethyLumiSet(nch,parallel=parallel,n.sd=n.sd,oob=oob,
                                      caller=deparse(match.call()))
    sampleNames(mlumi) = as.character(barcodes)
    sampleNames(mlumi@QC) = sampleNames(mlumi)
    rm(nch)
    gc()
    if(!is.null(pdat) & class(pdat)=='data.frame') {
      if(identical(rownames(pdat), sampleNames(mlumi))) pData(mlumi) = pdat
      else warning('Your pData does not have barcodes as row names. Skipping.')
    }
    return(mlumi)

} # }}}

lumiCSV <- function(barcodes, pdat=NULL, parallel=F, ...){ # {{{ one-liner
    as(methylumiCSV(barcodes=barcodes,pdat=pdat,parallel=parallel),'MethyLumiM')
} # }}}

## read "raw" M/U files like from GEO 
## 
MUtoDF <- function(filename) { # {{{ for raw MU dumps when no IDATs are present
  stopifnot(file.exists(filename))
  data.matrix(read.table(filename, header=TRUE, sep="\t", row.names=1))
} # }}}

methylumiMU <- function(filename,annotation='IlluminaHumanMethylation450k',...){# {{{
  x <- MUtoDF(filename)
  pats <- c(m='\\.Methylated\\.signal',
            u='\\.Unmethylated\\.Signal',
            p='\\.Detection\\.Pval')
  xx <- lapply(pats, function(y) data.matrix(x[, grep(y, ignore=T, names(x))]))
  for(i in names(pats)) {
    colnames(xx[[i]]) <- gsub(pats[i],'',colnames(xx[[i]]), ignore.case=T)
  }
  xx$betas <- xx$m / (xx$m + xx$u)
  stopifnot(identical(colnames(xx[[1]]), colnames(xx[[2]])))
  stopifnot(identical(rownames(xx[[1]]), rownames(xx[[2]])))
  stopifnot(identical(colnames(xx[[2]]), colnames(xx[[3]])))
  stopifnot(identical(rownames(xx[[2]]), rownames(xx[[3]])))
  history.submitted = as.character(Sys.time())
  aDat <- with(xx, 
               assayDataNew(betas=betas, methylated=m, unmethylated=u, pvals=p))
  x.lumi = new("MethyLumiSet", assayData=aDat)
  x.lumi@annotation <- annotation
  pdat <- data.frame(ID=colnames(xx[[i]]))
  rownames(pdat) <- pdat$ID
  pData(x.lumi) <- pdat 
  fdat <- data.frame(Probe_ID=rownames(xx[[1]]))
  fdat$type <- substr(fdat$Probe_ID, 1, 2)
  rownames(fdat) <- fdat$Probe_ID
  fData(x.lumi) <- fdat
  history.finished <- as.character(Sys.time())
  history.command <- 'methylumiMU()'
  x.lumi@history <- rbind(x.lumi@history, 
                          data.frame(submitted = history.submitted, 
                                     finished = history.finished, 
                                     command = history.command))
  return(x.lumi)
} # }}}

lumiMU <- function(filename,annotation='IlluminaHumanMethylation450k',...){# {{{
    as(methylumiMU(filename=filename, annotation=annotation), 'MethyLumiM')
} # }}}
