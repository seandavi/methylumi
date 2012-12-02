if(require('GenomicRanges')) {

  setAs("MIAME", "SimpleList", function(from) { # {{{
    to = list()
    for(i in slotNames(from)) if(i != '.__classVersion__') to[[i]]=slot(from, i)
    return(SimpleList(to))
  }) # }}}

  ## convenience functions 
  mOOB <- function(mset) assayDataElement(mset, 'methylated.OOB')
  uOOB <- function(mset) assayDataElement(mset, 'unmethylated.OOB')

  msetToSE <- function(from) { # {{{
    require(FDb.InfiniumMethylation.hg19) 
    chip=gsub('^IlluminaHumanMethylation','HM',gsub('k$','',annotation(from)))
    row.dat <- getPlatform(chip)
    asy.dat <- SimpleList()
    features <- intersect(featureNames(from), names(row.dat))
    if(is(from, 'MethyLumiM')) {
      asy.dat$mvals = assayDataElement(from, 'exprs')[features, ]
    } else if(is(from, 'MethyLumiSet')) {
      asy.dat$betas = assayDataElement(from, 'betas')[features, ]
    }
    if(all(c('methylated','unmethylated') %in% assayDataElementNames(from))){
      asy.dat$total = assayDataElement(from, 'methylated')[features, ] +
                      assayDataElement(from, 'unmethylated')[features, ]
    }
    SummarizedExperiment(assays=asy.dat,
                         rowData=row.dat[features],
                         colData=as(pData(from), 'DataFrame'),
                         exptData=as(experimentData(from), 'SimpleList'))
  } # }}}

  setAs("MethyLumiSet", "SummarizedExperiment", function(from) msetToSE(from))
  setAs("MethyLumiM", "SummarizedExperiment", function(from) msetToSE(from))

  if(require(minfi)) {

    require(FDb.InfiniumMethylation.hg19) ## and/or IlluminaManifest.foo.bar?

    byChannel.Both <- function(mset, annot) { # {{{
      probes <- names(split(annot, values(annot)$channel)[['Both']])
      probes <- intersect(probes, featureNames(mset))
      addrs <- values(annot[probes])$addressA
      Green = methylated(mset)[probes, ]
      Red = unmethylated(mset)[probes, ]
      rownames(Green) = rownames(Red) = addrs
      return(list(Green=Green, Red=Red))
    } # }}}
    byChannel.Grn <- function(mset, annot) { # {{{
      probes <- names(split(annot, values(annot)$channel)[['Grn']])
      probes <- intersect(probes, featureNames(mset))
      addrsA <- values(annot[probes])$addressA
      addrsB <- values(annot[probes])$addressB

      Grn.M = methylated(mset)[probes, ]
      rownames(Grn.M) = addrsB
      Grn.U = unmethylated(mset)[probes, ]
      rownames(Grn.U) = addrsA
      Green = rbind(Grn.M, Grn.U)

      Red.M = mOOB(mset)[probes, ]
      rownames(Red.M) = addrsB
      Red.U = uOOB(mset)[probes, ]
      rownames(Red.U) = addrsA
      Red = rbind(Red.M, Red.U)

      return(list(Green=Green, Red=Red))
    } # }}}
    byChannel.Red <- function(mset, annot) { # {{{
      probes <- names(split(annot, values(annot)$channel)[['Red']])
      probes <- intersect(probes, featureNames(mset))
      addrsA <- values(annot[probes])$addressA
      addrsB <- values(annot[probes])$addressB

      Grn.M = mOOB(mset)[probes, ]
      rownames(Grn.M) = addrsB
      Grn.U = uOOB(mset)[probes, ]
      rownames(Grn.U) = addrsA
      Green = rbind(Grn.M, Grn.U)

      Red.M = methylated(mset)[probes, ]
      rownames(Red.M) = addrsB
      Red.U = unmethylated(mset)[probes, ]
      rownames(Red.U) = addrsA
      Red = rbind(Red.M, Red.U)

      return(list(Green=Green, Red=Red))
    } # }}}
    byChannel <- function(mset, channel, annot) { # {{{
      stopifnot(channel %in% levels(as.factor(values(annot)$channel)))
      if(channel=='Both') byChannel.Both(mset, annot)
      else if(channel=='Grn') byChannel.Grn(mset, annot)
      else if(channel=='Red') byChannel.Red(mset, annot)
    } # }}}

    methylumiToMinfi <- function(from, annot=NULL) { # {{{
      if(!all(c('methylated','unmethylated','methylated.OOB','unmethylated.OOB')
              %in% assayDataElementNames(from))){
        stop('Cannot construct an RGChannelSet without full (OOB) intensities')
      }
      chip=gsub('^IlluminaHumanMethylation','HM',gsub('k$','',annotation(from)))
      if(is.null(annot)) annot <- getPlatform(chip)
    
      chs <- levels(as.factor(values(annot)$channel))
      names(chs) <- chs
      chans <- lapply(chs, byChannel, mset=from, annot=annot)
      Grn.ctls <- methylated(QCdata(from))
      Red.ctls <- methylated(QCdata(from))
      ctl.addrs <- as.character(fData(QCdata(from))$Address)
      rownames(Grn.ctls) <- rownames(Red.ctls) <- ctl.addrs
      Green = do.call( rbind, lapply(chans, function(x) x[['Green']]) )
      Green = rbind(Green, Grn.ctls)
      Red = do.call( rbind, lapply(chans, function(x) x[['Red']]) )
      Red = rbind(Red, Red.ctls)

      stopifnot(identical(rownames(Red), rownames(Green)))
      rg <- RGChannelSet(Green=Green, Red=Red, phenoData=phenoData(from))
      annotation(rg) <- annotation(from)
      return(rg)

    } # }}}

    setMethod("methylated",signature(object="MethylSet"),function(object)# {{{
      return(assayDataElement(object,"Meth"))) # }}}
    setMethod("unmethylated",signature(object="MethylSet"),function(object)# {{{
      return(assayDataElement(object,"Unmeth"))) # }}}
    setMethod("betas",signature(object="MethylSet"), function(object) # {{{
      return(minfi::getBeta(object))) # }}}
    setMethod("betas", signature(object="SummarizedExperiment"), # {{{
           function(object) assays(object)$betas ) # }}}
    setMethod("betas", signature(object="GenomicMethylSet"), # {{{
           function(object) minfi::getBeta(object)) # }}} 

    setAs("SummarizedExperiment", "GenomicMethylSet", function(from) { # {{{
      message('This function is almost solely for TCGA use... beware...')
      assaynames <- names(assays(from, withDimnames=F))
      stopifnot(all(c('betas','total') %in% assaynames) ||
                all(c('methylated','unmethylated') %in% assaynames))
      if(nrow(from) > 27578) { # {{{ 450k
        annot <- c(array="IlluminaHumanMethylation450k", annotation="ilmn.v1.2")
        prepro <- c(rg.norm='methylumi.bgcorr() + normalizeMethyLumiSet()')# }}}
      } else { # {{{
        annot <- c(array="IlluminaHumanMethylation27k", annotation="ilmn.v1.2")
        prepro <- c(rg.norm='methylumi.bgcorr()')
      } # }}}
      if('betas' %in% assaynames) { # {{{
        gm <- GenomicMethylSet(gr=rowData(from),
                               Meth=(assays(from, withDim=F)$betas* 
                                     assays(from, withDim=F)$total ),
                               Unmeth=((1-assays(from, withDim=F)$betas)*
                                       assays(from, withDim=F)$total),
                               pData=colData(from),
                               annotation=annot,
                               preprocessMethod=prepro) # }}}
      } else { # {{{
        gm <- GenomicMethylSet(gr=rowData(from),
                               Meth=assays(from, withDimnames=F)$methylated,
                               Unmeth=assays(from, withDimnames=F)$unmethylated,
                               pData=colData(from),
                               annotation=annot,
                               preprocessMethod=prepro)
      } # }}}
      exptData(gm) <- exptData(from)
      return(gm)
    }) # }}}

    setAs("MethyLumiSet", "RGChannelSet", function(from) methylumiToMinfi(from))
    setAs("MethyLumiM", "RGChannelSet", function(from) methylumiToMinfi(from))
    setAs("MethyLumiSet", "MethylSet", function(from) { # {{{
      pre <- c('methylumi')
      if(any(grepl('bgcorr', as.character(getHistory(from)$command)))) 
        pre <- c(pre, 'background corrected')
      if(any(grepl('ormalize', as.character(getHistory(from)$command))))
        pre <- c(pre, 'dye bias equalized')
      to <- MethylSet(Meth=methylated(from), Unmeth=unmethylated(from))
      to@annotation <- c(array=annotation(from), annotation='ilmn.v1.2')
      to@preprocessMethod <- c(rg.norm=paste(pre, collapse=', '),
                               minfi=paste(packageVersion('minfi'),
                                           collapse='.'), 
                               manifest='0.4')
      pData(to) <- pData(from)
      fData(to) <- fData(from)
      return(to)
    }) # }}}
    setAs("MethyLumiM", "MethylSet", function(from) { # {{{
      pre <- c('methylumi')
      if(any(grepl('bgcorr', as.character(getHistory(from)$command)))) 
        pre <- c(pre, 'background corrected')
      if(any(grepl('ormalize', as.character(getHistory(from)$command))))
        pre <- c(pre, 'dye bias equalized')
      to <- MethylSet(Meth=methylated(from), Unmeth=unmethylated(from))
      to@annotation <- c(array=annotation(from), annotation='ilmn.v1.2')
      to@preprocessMethod <- c(rg.norm=paste(pre, collapse=', '),
                               minfi=paste(packageVersion('minfi'),
                                           collapse='.'), 
                               manifest='0.4')
      pData(to) <- pData(from)
      fData(to) <- fData(from)
      return(to)
    }) # }}}
    setMethod("mapToGenome", signature(object="MethyLumiSet"), # {{{
          function(object, ...) {
            mapToGenome(as(object, 'MethylSet'))
          }) # }}}
    setMethod("mapToGenome", signature(object="MethyLumiM"), # {{{
          function(object, ...) {
            mapToGenome(as(object, 'MethylSet'))
          }) # }}}

  } # }}}

}
