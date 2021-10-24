## for HM27k, all probes are of "design I": single-channel, two-address pairings
## for HM450k, probes are either "design I" or "design II" as noted in manifest!
## IDATs from GoldenGate methylation arrays are not supported at this time.
cy3 <- function(object) { # {{{
  if(is.element('Color_Channel', fvarLabels(object)) && 
     !is.element('COLOR_CHANNEL', fvarLabels(object))) {
    fvarLabels(object)<-gsub('Color_Channel','COLOR_CHANNEL',fvarLabels(object))
  }
  if(!is.element('COLOR_CHANNEL', fvarLabels(object))) {
    annotChip <- paste(annotation(object),'db',sep='.')
    object <- addColorChannelInfo(object, annotChip)
  }
  return(which(fData(object)$COLOR_CHANNEL=='Grn'))
} # }}}
cy5 <- function(object) { # {{{
  if(is.element('Color_Channel', fvarLabels(object)) && 
     !is.element('COLOR_CHANNEL', fvarLabels(object))) {
    fvarLabels(object)<-gsub('Color_Channel','COLOR_CHANNEL',fvarLabels(object))
  }
  if(!is.element('COLOR_CHANNEL', fvarLabels(object))) {
    annotChip <- paste(annotation(object),'db',sep='.')
    object <- addColorChannelInfo(object, annotChip)
  }
  return(which(fData(object)$COLOR_CHANNEL=='Red'))
} # }}}

## Utility function for dealing with single samples (still not 100% perfect...)
columnMatrix <- function(x, row.names=NULL) { # {{{
  if(is.null(dim(x)[2])) dim(x) = c( length(x), 1 )
  if(!is.null(row.names)) rownames(x) = row.names
  return(x)
} # }}}

## require()s the appropriate package for annotating a chip & sets up mappings
getMethylationBeadMappers <- function(chipType=c('450k','27k'), genome=c('hg19','hg18')) { # {{{
  
  genome <- match.arg(genome) ## default to FDb.InfiniumMethylation.hg19
  pkg <- paste0('FDb.InfiniumMethylation.', genome)
  require(pkg, character.only=TRUE) ## and

  chipType <- sub('^IlluminaHumanMethylation', '', chipType)
  if(class(chipType) %in% c('NChannelSet','MethyLumiSet','MethyLumiM')) {
    chipType <- sub('^IlluminaHumanMethylation', '', annotation(chipType))
    chipType <- sub('.db$', '', annotation(chipType))
  }
  chipType <- match.arg(chipType) # backwards compatibility purposes

  getControls <- switch(chipType,
                        '27k'=function() {
                          data(hm27.controls)
                          return(hm27.controls)
                        },
                        '450k'=function() {
                          data(hm450.controls)
                          return(hm450.controls)
                        })
 
  ## addressA=U, addressB=M
  getProbes <- switch(chipType,
                      '27k'=function(color=NULL, ...) {
                              data(hm27.ordering)
                              what <- c('Probe_ID','M','U')
                              r <- split(hm27.ordering[,what],
                                         hm27.ordering$col)
                              if (is.null(color)) return(r)
                              else return(r[[substr(color,1,1)]]) 
                            },
                      '450k'=function(design=NULL, color=NULL, ...) {
                               data(hm450.ordering)
                               what <- c('Probe_ID','M','U')
                               r <- split(hm450.ordering[,c(what,'col')],
                                          hm450.ordering$DESIGN)
                               r$I <- split(r$I[, what], r$I$col)
                               r$II <- r$II[, what]
                               if (is.null(design)) {
                                 return(r)
                               } else if (design == 'I' && !is.null(color)) {
                                 return(r[[design]][[substr(color,1,1)]]) 
                               } else { 
                                 return(r[[design]])
                               }
                             })

  ord <- c('Probe_ID','DESIGN','COLOR_CHANNEL')  
  getOrdering <- switch(chipType,
                        '27k'=function() return(hm27.ordering[, ord]),
                        '450k'=function() return(hm450.ordering[, ord]))

  mapper <- list(probes=getProbes,
                 controls=getControls,
                 ordering=getOrdering)
  return(mapper)

} # }}}

#' process a single IDAT (just the mean intensities)
#' 
#' @param x character(1)
#' @param fileExts named list
#' @param idatPath character(1)
#'  
IDATtoMatrix <- function(x,fileExts=list(Cy3="Grn",Cy5="Red"),idatPath='.'){#{{{
  chs = names(fileExts)
  names(chs) = fileExts
  processed = lapply(fileExts, function(ch) {
    ext = paste(ch, 'idat', sep='.')
    dat = readIDAT(file.path(idatPath, paste(x, ext, sep='_')))
    Quants = data.matrix(dat$Quants)
    colnames(Quants) = paste(chs[ch], colnames(Quants), sep='.')
    return(list(Quants=Quants, 
                RunInfo=dat$RunInfo,
                ChipType=dat$ChipType))
  })
  probe.data = do.call(cbind, lapply(processed, function(x) x[['Quants']]))
  probe.data <- probe.data[ , c('Cy3.Mean','Cy5.Mean')]
  attr(probe.data, 'RunInfo') = processed[[1]][['RunInfo']]
  attr(probe.data, 'ChipType') = processed[[1]][['ChipType']]
  return(probe.data)
} # }}}

#' convert multiple idats to matrices
#' 
#' @param barcodes character()
#' @param fileExts character()
#' @param parallel logical(1)
#' @param idatPath character(1)
#' 
IDATsToMatrices <- function(barcodes, fileExts=list(Cy3="Grn", Cy5="Red"), parallel=F, idatPath='.') { # {{{
  names(barcodes) = as.character(barcodes)
  if(parallel) {
    mats = .mclapply(barcodes,IDATtoMatrix,fileExts=fileExts,idatPath=idatPath)
  } else {
    mats = lapply(barcodes, IDATtoMatrix, fileExts=fileExts, idatPath=idatPath)
  }
  names(mats) = as.character(barcodes)
  return(mats)
} # }}}

## might consider doing the following in C++ via Rcpp to avoid copying data
extractAssayDataFromList <- function(assay, mats, fnames) { # {{{
  d <- do.call('cbind', lapply(mats, function(x) x[ , assay ]))
  if (!is.matrix(d)) d <- data.matrix(d)
  colnames(d) <- names(mats)
  rownames(d) <- fnames
  return(d)
} # }}}

## a faster rewrite of DFsToNChannelSet() so that I can decommission it...
DataToNChannelSet <- function(mats, chans=c(Cy3='GRN',Cy5='RED'), parallel=F, protocol.data=F, IDAT=TRUE){ # {{{

  stopifnot(is(mats, 'list'))
  assayNames = paste0(names(chans), '.Mean')
  names(assayNames) = assayNames
  fnames <- rownames(mats[[1]]) 
  extract <- function(assay) extractAssayDataFromList(assay, mats, fnames)  
  assays = lapply( assayNames, extract)
  obj = new("NChannelSet",
             assayData=assayDataNew(R=assays[['Cy5.Mean']],
                                    G=assays[['Cy3.Mean']]))
  featureNames(obj) = rownames(mats[[1]])
  if(IDAT) { # {{{
    message('Attempting to extract protocolData() from list...')
    ChipType = attr(mats[[1]], 'ChipType')
    RunInfo = lapply(mats, function(d) attr(d, 'RunInfo'))
    if(protocol.data) { # {{{
      scanDates = data.frame(DecodeDate=rep(NA, length(mats)),
                             ScanDate=rep(NA, length(mats)))
      rownames(scanDates) = names(mats)
      for(i in seq_along(mats)) {
        cat("decoding protocolData for", names(mats)[i], "...\n")
        if(nrow(RunInfo[[i]]) >= 2) {
          scanDates$DecodeDate[i] = RunInfo[[i]][1,1]
          scanDates$ScanDate[i]  =  RunInfo[[i]][2,1]
        }
      }
      protocoldata = new("AnnotatedDataFrame",
                          data=scanDates,
                          varMetadata=data.frame(
                            labelDescription=colnames(scanDates),
                            row.names=colnames(scanDates)
                           )
                          )
      protocolData(obj) = protocoldata
    } # }}}
  
    message('Determining chip type from IDAT protocolData...')
    if(ChipType == "BeadChip 12x1") {
      annotation(obj) = 'IlluminaHumanMethylation27k'
    } else if(ChipType == "BeadChip 12x8") {
      annotation(obj) = 'IlluminaHumanMethylation450k'
    }

  } # }}}
  if(is.null(annotation(obj))) { # {{{
    if(dim(obj)[1] == 55300) annotation(obj) = 'IlluminaHumanMethylation27k'
    else annotation(obj) = 'IlluminaHumanMethylation450k'
  } # }}}
  return(obj)

} # }}}

getControlProbes <- function(NChannelSet) { # {{{

  fD <- getMethylationBeadMappers(annotation(NChannelSet))$controls()
  ctls <- match(fD[['Address']], featureNames(NChannelSet))

  ## FIXME: make this happen in the annotations, to avoid redundancy in names!
  rownames(fD) <- ctlnames <- make.names(fD[,'Name'], unique=T)
  fvD <- data.frame(labelDescription=c(
        'Address of this control bead',
        'Purpose of this control bead',
        'Color channel for this bead',
        'Reporter group ID for this bead'
      ))
  fDat <- new("AnnotatedDataFrame", data=fD, varMetadata=fvD)
  methylated <- assayDataElement(NChannelSet,'G')[ctls,,drop=FALSE] # Cy3
  unmethylated <- assayDataElement(NChannelSet,'R')[ctls,,drop=FALSE] # Cy5
  rownames(methylated) <- rownames(unmethylated) <- ctlnames

  aDat <- assayDataNew(methylated=methylated, 
                       unmethylated=unmethylated)
  new("MethyLumiQC", assayData=aDat, featureData=fDat, 
                     annotation=annotation(NChannelSet))

} # }}}

## 27k design, both probes same channel; ~100,000 of the 450k probes as well
designItoMandU <- function(NChannelSet, parallel=F, n=F, n.sd=F, oob=T) { # {{{

  mapper <- getMethylationBeadMappers(annotation(NChannelSet))
  probes <- mapper$probes(design='I') # as list(G=..., R=...)
  channels <- c('G','R')
  names(channels) <- channels

  getIntCh <- function(NChannelSet, ch, al) { # {{{
    a = assayDataElement(NChannelSet,ch)[as.character(probes[[ch]][[al]]),,drop=FALSE]
    rownames(a) = as.character(probes[[ch]][['Probe_ID']])
    return(a)
  } # }}}

  getOOBCh <- function(NChannelSet, ch, al) { # {{{
    ch.oob <- ifelse(ch == 'R', 'G', 'R')
    a = assayDataElement(NChannelSet,ch.oob)[as.character(probes[[ch]][[al]]),,drop=FALSE]
    rownames(a) = as.character(probes[[ch]][['Probe_ID']])
    return(a)
  } # }}}

  getAllele <- function(NChannelSet, al, parallel=F, n=n, n.sd=T, oob=T) { # {{{
    fluor = lapply(channels, function(ch) getIntCh(NChannelSet, ch, al))
    fluor.oob = lapply(channels, function(ch) getOOBCh(NChannelSet, ch, al))
    res = list()
    res[[ 'I' ]] = fluor
    if(oob)  res[[ 'OOB' ]] = fluor.oob
    lapply(res, function(r) {
      names(r) = channels
      return(r)
    })
  } # }}}

  signal <- lapply(c(M='M',U='U'), function(al) {
    getAllele(NChannelSet, al, parallel=F, n=n, n.sd=n.sd, oob=oob)
  })

  retval = list(
    methylated=rbind(signal$M$I$R, signal$M$I$G),
    unmethylated=rbind(signal$U$I$R, signal$U$I$G)
  )
  if(oob) {
    retval[['methylated.OOB']] = rbind(signal$M$OOB$R, signal$M$OOB$G)
    retval[['unmethylated.OOB']] = rbind(signal$U$OOB$R, signal$U$OOB$G)
  }

  return(retval)

} # }}}

## 450k/GoldenGate design (green=methylated, red=unmethylated, single address)
designIItoMandU <- function(NChannelSet, parallel=F, n=F, n.sd=F, oob=T) { # {{{

  ## loads the annotation DB so we can run SQL queries
  mapper <- getMethylationBeadMappers(annotation(NChannelSet))
  probes2 <- mapper$probes(design='II')
  probes2$M <- probes2$U ## horrid kludge

  getIntCh <- function(NChannelSet, ch=NULL, al) { # {{{
    ch <- ifelse(al=='M', 'G', 'R')
    a <- assayDataElement(NChannelSet,ch)[as.character(probes2[[al]]), , drop=F]
    rownames(a) <- as.character(probes2[['Probe_ID']])
    return(a)
  } # }}}

  getAllele <- function(NChannelSet, al, n=F, n.sd=F, oob=F) { # {{{

    ch <- ifelse(al=='M', 'G', 'R')
    res <- list()
    res[['I']] <- getIntCh(NChannelSet,ch,al)
    if(oob) {
      res[['OOB']] <- res[['I']]
      is.na(res[['OOB']]) <- TRUE 
    }  
    return(res)

  } # }}}
  
  ## M == Grn/Cy3 and U == Red/Cy5, same address
  ##
  alleles = c(M='M',U='U')
  signal = lapply(alleles, function(a) getAllele(NChannelSet,a,n,n.sd,oob))
  
  retval = list( methylated=signal$M$I, unmethylated=signal$U$I )
  if(oob) {
    retval[['methylated.OOB']] = signal$M$OOB
    retval[['unmethylated.OOB']] = signal$U$OOB
  }
  return(retval)

} # }}}

## 12/12/14: this code is hideous, what sort of clown wrote it?  Oh yeah, I did
mergeProbeDesigns <- function(NChannelSet, parallel=F, n=F, n.sd=F, oob=T){ #{{{
  
  if(annotation(NChannelSet) == 'IlluminaHumanMethylation450k') {
    design1=designItoMandU(NChannelSet,parallel=parallel,n=n,n.sd=n.sd,oob=oob)
    ## this is the source of the problem currently:
    design2=designIItoMandU(NChannelSet,parallel=parallel,n=n,n.sd=n.sd,oob=oob)
    res <- list()
    for(i in names(design1)) {
      res[[i]] <- rbind(design1[[i]], design2[[i]])
      rownames(res[[i]]) <- c( rownames(design1[[i]]), rownames(design2[[i]]) )
    }
  } else if(annotation(NChannelSet) == 'IlluminaHumanMethylation27k') {
    res <- designItoMandU(NChannelSet, parallel=parallel,n=n,n.sd=n.sd,oob=oob)
  } else {
    stop("don't know how to process chips of type", annotation(NChannelSet))
  }
  return(res) # reorder on the way out...

} # }}}

NChannelSetToMethyLumiSet <- function(NChannelSet, parallel=F, normalize=F, pval=0.05, n=F, n.sd=F, oob=T, caller=NULL){ # {{{

  history.submitted = as.character(Sys.time())

  results = mergeProbeDesigns(NChannelSet,parallel=parallel,n.sd=n.sd,oob=oob)
  if(oob) {
    aDat <- with(results,
              assayDataNew(methylated=methylated, 
                           unmethylated=unmethylated,
                           methylated.OOB=methylated.OOB,
                           unmethylated.OOB=unmethylated.OOB,
                           betas=methylated/(methylated+unmethylated),
                           pvals=methylated/(methylated+unmethylated)))
                           # pvals are a cheat to force pval.detect()
  } else {
    aDat <- with(results,
              assayDataNew(methylated=methylated, 
                           unmethylated=unmethylated,
                           betas=methylated/(methylated+unmethylated),
                           pvals=methylated/(methylated+unmethylated)))
  }
  rm(results)
  gc()

  ## now return the MethyLumiSet (which can be directly coerced to MethyLumiM)
  x.lumi = new("MethyLumiSet", assayData=aDat)
  x.lumi@QC <- getControlProbes(NChannelSet)
  x.lumi@protocolData <- protocolData(NChannelSet)
  x.lumi@annotation <- annotation(NChannelSet)
  x.lumi@QC@annotation <- annotation(NChannelSet)
  pdat <- data.frame(barcode=sampleNames(NChannelSet))
  rownames(pdat) <- sampleNames(NChannelSet)
  pData(x.lumi) <- pdat 
  varLabels(x.lumi) <- c('barcode')
  varMetadata(x.lumi)[,1] <- c('Illumina BeadChip barcode')
  mapper <- getMethylationBeadMappers(annotation(NChannelSet))
  fdat <- mapper$ordering()
  rownames(fdat) <- fdat$Probe_ID
  x.fnames <- rownames(betas(x.lumi))
  fdat <- fdat[ x.fnames, ]

  ## Regression tests: fail noisily if there is an ordering issue
  stopifnot(identical(rownames(methylated(x.lumi)),
                      rownames(unmethylated(x.lumi))))
  stopifnot(identical(rownames(betas(x.lumi)), 
                      rownames(unmethylated(x.lumi))))
  stopifnot(identical(rownames(fdat), 
                      rownames(betas(x.lumi))))

  fData(x.lumi) <- fdat
  possibleLabels <- c('Probe_ID',
                      'DESIGN',
                      'COLOR_CHANNEL',
                      'PROBE_TYPE',
                      'SNP10',
                      'SYMBOL',
                      'CHR36',
                      'CPG36',
                      'CPGS')
  fvarLabels(x.lumi) <- possibleLabels[ 1:ncol(fdat) ]
  possibleMetadata <- c('Illumina probe ID from manifest',
                        'Infinium design type (I or II)',
                        'Color channel (for type I probes)',
                        'Probe locus type (CpG, CpH, or SNP)',
                        'SNP (dbSNP build 128) within 10bp of target?',
                        'Gene symbol (if probe is annotated to a gene)',
                        'Chromosome mapping for probe in hg18 assembly',
                        'Coordinates of interrogated cytosine in hg18',
                        'Number of CpG dinucleotides in probe sequence')
  fvarMetadata(x.lumi)[,1] <- possibleMetadata[ 1:ncol(fdat) ]
  pval.detect(x.lumi) <- pval # default value
  history.finished <- as.character(Sys.time())
  history.command <- ifelse(is.null(caller),'NChannelSet(x)',caller)
  x.lumi@history <- rbind(x.lumi@history, 
                          data.frame(submitted = history.submitted, 
                                     finished = history.finished, 
                                     command = history.command))
  if(normalize) x.lumi = normalizeMethyLumiSet(x.lumi)
  return(x.lumi)

} # }}}

methylumIDAT <- function(barcodes=NULL,pdat=NULL,parallel=F,n=F,n.sd=F,oob=T,idatPath=getwd(), ...) { # {{{
  if(is(barcodes, 'data.frame')) pdat = barcodes
  if((is.null(barcodes))&(is.null(pdat) | (!('barcode' %in% names(pdat))))){#{{{
    stop('"barcodes" or "pdat" (with pdat$barcode defined) must be supplied.')
  } # }}}
  if(!is.null(pdat) && 'barcode' %in% tolower(names(pdat))) { # {{{
    names(pdat)[ which(tolower(names(pdat))=='barcode') ] = 'barcode'
    barcodes = pdat$barcode
    if(any(grepl('idat',ignore.case=T,barcodes))) { 
      message('Warning: filtering out raw filenames') 
      barcodes = gsub('_(Red|Grn)','', barcodes, ignore.case=TRUE)
      barcodes = gsub('.idat', '', barcodes, ignore.case=TRUE)
    }
    if(any(duplicated(barcodes))) {  
      message('Warning: filtering out duplicates') 
      pdat = pdat[ -which(duplicated(barcodes)), ] 
      barcodes = pdat$barcode
    } # }}}
  } else { # {{{
    if(any(grepl('idat',ignore.case=T,barcodes))) { 
      message('Warning: filtering out raw filenames') 
      barcodes = unique(gsub('_(Red|Grn)','', barcodes, ignore.case=TRUE))
      barcodes = unique(gsub('.idat','', barcodes, ignore.case=TRUE))
    }
    if(any(duplicated(barcodes))) { 
      message('Warning: filtering out duplicate barcodes')
      barcodes = barcodes[ which(!duplicated(barcodes)) ] 
    } 
  } # }}}
  files.present = rep(TRUE, length(barcodes)) # {{{
  idats = sapply(barcodes, function(b) paste(b,c('_Red','_Grn'),'.idat',sep=''))
  for(i in colnames(idats)) for(j in idats[,i]) if(!j %in% list.files(idatPath))  {
    message(paste('Error: file', j, 'is missing for sample', i))
    files.present = FALSE
  }
  stopifnot(all(files.present)) # }}}
  hm27 = hm450 = 0 # {{{
  hm27 = sum(!grepl('_R0[123456]C0[12]', barcodes)) 
  message(paste(hm27, 'HumanMethylation27 samples found'))
  hm450 = sum(grepl('_R0[123456]C0[12]', barcodes))
  message(paste(hm450, 'HumanMethylation450 samples found'))
  if( hm27 > 0 && hm450 > 0 ) {
    stop('Cannot process both platforms simultaneously; please run separately.')
  } # }}}

  mats <- IDATsToMatrices(barcodes, parallel=parallel, idatPath=idatPath) 
  dats <- DataToNChannelSet(mats, IDAT=T, parallel=parallel)
  mlumi <- NChannelSetToMethyLumiSet(dats, parallel=parallel, oob=oob, 
                                     caller=deparse(match.call()))

  if(is.null(pdat)) { # {{{
    pdat = data.frame(barcode=as.character(barcodes))
    rownames(pdat) = pdat$barcode
    pData(mlumi) = pdat # }}}
  } else { # {{{
    pData(mlumi) = pdat
  } # }}}
  if(!is.null(mlumi@QC)) { #{{{ should be gratuitous now
    sampleNames(mlumi@QC) = sampleNames(mlumi)
  } # }}}

  # finally
  return(mlumi[ sort(featureNames(mlumi)), ])

} # }}}

lumIDAT <- function(barcodes, pdat=NULL, parallel=F, n=T, idatPath=getwd(), ...){ # {{{ 
  as(methylumIDAT(barcodes=barcodes,pdat=pdat,parallel=parallel,n=n,oob=F,idatPath=idatPath),
     'MethyLumiM')
} # }}}
