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
