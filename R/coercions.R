
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


    getPlatform <- function(platform="HM450", genome="hg19") { ## {{{ 
      require(Biostrings)
      require(rtracklayer)
      if (genome == "hg19") {
        message("Fetching coordinates for hg19...")
        if(require(FDb.InfiniumMethylation.hg19)) {
          GR <- features(FDb.InfiniumMethylation.hg19)
        } else {
          message("The FDb.InfiniumMethylation.hg19 package appears to be unavailable.  Please install it to use hg19 coordinates")
          stop()
        }
      } else if (genome == "hg18") {
        # FDb.InfiniumMethyulation.hg18 is gone as of bioc 2.12
        stop('hg18 is no longer supported')
      } else {
        stop("Only hg19 is currently supported.")
      }
      GR <- keepSeqlevels(GR, paste0("chr", c(1:22, "X", "Y")))
      if ("name" %in% names(mcols(GR))) names(GR) <- mcols(GR)$name
      if (is.na(unique(genome(GR)))) genome(GR) <- genome
      seqinfo(GR) <- SeqinfoForBSGenome(unique(genome(GR)))[seqlevels(GR)]
      platform <- toupper(platform)
      if (platform %in% c("HM450", "ILLUMINAHUMANMETHYLATION450")) {
        hm450.controls <- NULL
        data(hm450.controls)
        GR <- GR[which(mcols(GR)$platform %in% c("HM450", "BOTH"))]
        mcols(GR)$channel <- Rle(as.factor(mcols(GR)$channel450))
        mcols(GR)$addressA <- as.character(mcols(GR)$addressA_450)
        mcols(GR)$addressB <- as.character(mcols(GR)$addressB_450)
        attr(GR, "controls") <- hm450.controls
      } else if (platform == "HM27") {
        hm27.controls <- NULL
        data(hm27.controls)
        GR <- GR[which(mcols(GR)$platform %in% c("HM27", "BOTH"))]
        mcols(GR)$channel <- Rle(as.factor(mcols(GR)$channel27))
        mcols(GR)$addressA <- as.character(mcols(GR)$addressA_27)
        mcols(GR)$addressB <- as.character(mcols(GR)$addressB_27)
        attr(GR, "controls") <- hm27.controls
      } else {
        stop("You need to specify either HM27 or HM450 as platform to run")
      }
      mcols(GR)$percentGC <- as.numeric(mcols(GR)$percentGC)
      mcols(GR)$probeType <- Rle(as.factor(mcols(GR)$probeType))
      mcols(GR)$platform <- Rle(as.factor(mcols(GR)$platform))
      mcols(GR)$sourceSeq <- DNAStringSet(mcols(GR)$sourceSeq)
      kept = c("addressA", "addressB", "channel", "platform", "percentGC", 
               "sourceSeq","probeType","probeStart","probeEnd","probeTarget")
      val <- mcols(GR)[, intersect(names(mcols(GR)), kept)]
      mcols(GR) <- val
      return(GR)
    } # }}}

    byChannel.Both <- function(mset, annot) { # {{{
      probes <- names(split(annot, values(annot)$channel)[['Both']])
      probes <- intersect(probes, featureNames(mset))
      addrs <- values(annot[probes])$addressA
      Green = methylated(mset)[probes, ]
      Red = unmethylated(mset)[probes, ]
      rownames(Green) = rownames(Red) = as.character(addrs)
      return(list(Green=Green, Red=Red))
    } # }}}
    byChannel.Grn <- function(mset, annot) { # {{{
      probes <- names(split(annot, values(annot)$channel)[['Grn']])
      probes <- intersect(probes, featureNames(mset))
      addrsA <- values(annot[probes])$addressA
      addrsB <- values(annot[probes])$addressB

      Grn.M = methylated(mset)[probes, ]
      rownames(Grn.M) = as.character(addrsB)
      Grn.U = unmethylated(mset)[probes, ]
      rownames(Grn.U) = as.character(addrsA)
      Green = rbind(Grn.M, Grn.U)

      Red.M = mOOB(mset)[probes, ]
      rownames(Red.M) = as.character(addrsB)
      Red.U = uOOB(mset)[probes, ]
      rownames(Red.U) = as.character(addrsA)
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
            minfi::mapToGenome(as(object, 'MethylSet'))
          }) # }}}
    setMethod("mapToGenome", signature(object="MethyLumiM"), # {{{
          function(object, ...) {
            minfi::mapToGenome(as(object, 'MethylSet'))
          }) # }}}

