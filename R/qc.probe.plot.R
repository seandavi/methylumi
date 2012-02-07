setClassUnion("methylData", c('MethyLumiSet','MethyLumiM'))
setClassUnion("ND", c('character','missing'))

# fix for QC plots in methylumi (hard to read on Infinium arrays) using ggplot2
qc.probe.plot <- function(obj,controltype="negnorm",log2=T,by.type=F,...){ # {{{
  require("ggplot2")
  if( class(obj) %in% c('MethyLumiSet','MethyLumiM') ) {
    qc <- controlData(obj)
    if(!identical(sampleNames(qc), sampleNames(obj))) {
      sampleNames(qc) <- sampleNames(obj)
    }
  } else if(class(obj) == 'MethyLumiOOB' & tolower(controltype) == 'oob') {
    qc <- obj
  } else if(class(obj) == 'MethyLumiQC') {
    qc <- obj
  } else {
    stop("Don't know how to QC this data you've given me...")
  }
  if( tolower(controltype) == 'negnorm' || missing(controltype) ) {
    rows <- grep('(Negative|Norm)', fData(qc)$Type, ignore.case=T)
  } else { 
    rows <- grep(paste('^',controltype,sep=''),fData(qc)$Type,ignore.case=T)
  }
  if( tolower(controltype) == 'oob' ) {
    # {{{ out-of-band intensities
    dat <- intensities.OOB.allelic(obj)
    rownames(dat$Cy5$M) = paste(rownames(dat$Cy5$M), 'M', sep='_')
    rownames(dat$Cy5$U) = paste(rownames(dat$Cy5$U), 'U', sep='_')
    rownames(dat$Cy3$M) = paste(rownames(dat$Cy3$M), 'M', sep='_')
    rownames(dat$Cy3$U) = paste(rownames(dat$Cy3$U), 'U', sep='_')
    dat$Cy5 = rbind(dat$Cy5$M, dat$Cy5$U)
    dat$Cy3 = rbind(dat$Cy3$M, dat$Cy3$U)
    names(dat) = gsub('Cy3','green', gsub('Cy5','red', names(dat)))
    dat <- lapply(dat, function(d) { 
      datum = data.frame(d)
      colnames(datum) = gsub('^X', '', colnames(datum))
      m.probes = grepl('_M$', rownames(d))
      datum$probe = as.factor(gsub('_(M|U)$','',rownames(datum)))
      datum$type = 'unmethylated'
      datum$type[which(m.probes)] = 'methylated'
      datum$type = as.factor(datum$type)
      return(datum)
    })
    dat$green$channel = 'Cy5'
    dat$red$channel = 'Cy3'
    # }}}
    dat.frame <- rbind(dat$red, dat$green)
    dat.frame$channel <- as.factor(dat.frame$channel)
    a.title <- "Out-of-band probe intensity plot"
  } else { 
    # {{{ actual control probes...
    probes <- featureNames(qc)[ rows ]
    type <- as.factor(fData(qc)$Type[ rows ])
    if( tolower(controltype) == 'negnorm' ) { # {{{
      if( 'NORM_C' %in% unique(fData(qc)$Type) ) { # 450k array
        colour.settings <- c(NEGATIVE='darkgray',
                             NORM_A='red',
                             NORM_T='darkred',
                             NORM_C='green',
                             NORM_G='darkgreen')
        type = factor(type, levels=names(colour.settings))
        shape.settings <- c(20, 20, 20, 20, 20)
      } else { # figure it's a 27k array
        colour.settings <- c('darkgray','darkgreen','red')
        shape.settings <- c(20, 20, 20)
      }
    } # }}}
    dat <- c(red = "unmethylated", green = "methylated")
    Cy5 <- as.data.frame(assayDataElement(qc, dat[1])[rows, ])
    Cy5$channel <- "Cy5 (Red)"
    Cy5$probe <- as.factor(probes)
    Cy3 <- as.data.frame(assayDataElement(qc, dat[2])[rows, ])
    Cy3$channel <- "Cy3 (Green)"
    Cy3$probe <- as.factor(probes)
    Cy3$type <- Cy5$type <- type
    # }}}
    dat.frame <- rbind(Cy5, Cy3)
    dat.frame$channel <- as.factor(dat.frame$channel)
    a.title <- paste(controltype, "control probe plot")
    if (tolower(controltype) == 'negnorm') {
      a.title <- 'Negative & normalization control probes'
    }
  }
  more.args = list(...)
  if('extra' %in% names(more.args)) a.title=paste(a.title,more.args[['extra']])
  qc <- melt(dat.frame, id = c("probe", "channel", "type"))
  geometry <- ifelse(tolower(controltype) == 'negnorm',
                     ifelse(tolower(controltype)=='oob', 
                            'boxplot', 
                            'jitter'),
                     'point')
  if( tolower(controltype) == 'oob' ) {
    qc$grouping = paste(qc$variable, qc$type, sep='.')
    p <- ggplot2::qplot(data = qc, colour = type, fill = type, group = grouping,
                        x = variable, y = value, geom = 'boxplot', main=a.title,
                        xlab="Sample", ylab="Intensities") + 
                        coord_flip() 
    if (log2) p <- p + scale_y_log2()
  } else {
    p <- ggplot2::qplot(data = qc, colour = type, shape = type, 
                        x = value, y = variable, geom = geometry,
                        main=a.title, ylab="Sample", xlab="Intensities")
    p <- p + scale_y_discrete( limits=rev(sampleNames(obj)) ) # more readable
    if (log2) p <- p + scale_x_log2()
  }
  if( by.type ) {
    p <- p + facet_grid( type ~ channel)
  } else { 
    p <- p + facet_grid( . ~ channel)
  }
  if( tolower(controltype) == 'negnorm' ) {
    p <- p + scale_colour_manual( value=colour.settings )
    p <- p + scale_shape_manual( value=shape.settings )
  }
  p <- p + theme_bw()
  return( p )
} # }}}

methylumi.diagnostics <- function (x, onlybg=F) { # {{{
  x.qc <- controlData(x)
  if (!is.null(x.qc)) { # {{{ realistically, use OOB
      bg <- list(red = log2(negctls(x.qc, "Cy5")), 
                 green = log2(negctls(x.qc, "Cy3")))
      message("(Using negative controls for dashed vertical background line)")
  } # }}}
  if(!('COLOR_CHANNEL' %in% fvarLabels(x))) { # {{{
    if(annotation(x) == 'IlluminaHumanMethylation27k') { 
      fData(x)$COLOR_CHANNEL = mget(featureNames(x),
                                    IlluminaHumanMethylation27kCOLORCHANNEL)
    } else if(annotation(x) == 'IlluminaHumanMethylation450k') { 
      fData(x)$COLOR_CHANNEL = mget(featureNames(x),
                                    IlluminaHumanMethylation450kCOLORCHANNEL)
    }
  } # }}}
  assays <- c("exprs", "methylated", "unmethylated")
  if(onlybg) assays <- c("methylated", "unmethylated")
  is.450k = (annotation(x) == 'IlluminaHumanMethylation450k')
  par(mfrow = c(length(assays), 2 + is.450k))
  dye <- list('red'='Cy5','green'='Cy3')
  if(is.450k) dye[['both']] = 'Design II'
  for(assay in assays) {
    for(channel in names(dye)) {
      if(channel=="red") { # {{{
        probes = which(fData(x)[['COLOR_CHANNEL']]=='Grn')
        chcolor = channel # }}}
      } else if(channel=="green") { # {{{
        probes = which(fData(x)[['COLOR_CHANNEL']]=='Red')
        chcolor = channel # }}}
      } else if(channel=="both") { # {{{
        probes=which(fData(x)[['COLOR_CHANNEL']]=='Both')
        chcolor = ifelse(assay == 'exprs', 'blue', 
                         ifelse(assay == 'methylated', 'green', 'red'))
      } # }}}
      if (assay == "exprs") { # {{{
        if (is(x, "MethyLumiM")) { # {{{
          dat <- mvals(x)[probes, ]
          cutpoint <- 0
          xlab <- "M-value"
          xlim <- c(min(dat), max(dat)) # }}}
        } else { # {{{
          dat <- betas(x)[probes, ]
          cutpoint <- 0.5
          xlab <- "Beta value"
          xlim <- c(0, 1)
        } # }}}
        Nx <- 255
        scheme <- c("lightblue", "gray", "yellow")
        colrs <- colorRampPalette(scheme, space = "Lab")(Nx)[1:Nx]
        dens <- apply(dat, 2, function(x) density(na.omit(x)))
        densmax <- max(unlist(lapply(dens, function(x) x[["y"]])))
        dx <- density(na.omit(dat), n = Nx)
        plot(dx$x, dx$y, col = colrs, xlim = xlim, ylim = c(0, 
             densmax), xlab = xlab, ylab='density', type = "h")
        for (i in 1:length(dens)) lines(dens[[i]], lty=2, col=chcolor)
        title(paste(xlab, ":", dye[[channel]], "probes")) # }}}
      } else { # {{{ methylated/unmethylated
        dat <- log2(assayDataElement(x, assay)[probes,])
        plot.density(density(na.omit(dat)), col = "white", xlim = c(0, 16), 
                     ylim = c(0, 0.8), lwd = 1, xlab = "log2(Intensity)", 
                     ylab = "Proportion", main = "", lty = 1)
        for(i in 1:dim(dat)[2]) lines(density(dat[, i]), col=chcolor, lty=1)
        if (!is.null(x.qc)) {
          if(channel=='both' && assay=='methylated') bgch = bg[['green']]
          else if(channel=='both' && assay=='unmethylated') bgch = bg[['red']]
          else bgch = bg[[channel]]
          for(bgmean in colMeans(bgch, na.rm=TRUE)) {
            abline(v=bgmean, col=paste("dark",chcolor,sep=""), lty=3, lwd=1)
          }
          title(paste(assay, "intensities:", dye[[channel]]))
        }
      } # }}}
    }
  }
} # }}}

compare.chips <- function(x, by.design=F, by.channel=F, what='betas', how=NULL, plot.type='histogram', title=NULL) { # {{{
 
  require(ggplot2) 
  chips = dim(x)[2]
  probes = dim(x)[1]
  dat = data.frame(probe=c(rep(featureNames(x), chips)), # {{{
                   assay=as.vector(switch(what, 
                                          'mvals' = exprs(x),
                                          assayDataElement(x, what))),
                   chip=as.vector(sapply(sampleNames(x),rep,times=probes)))# }}}

  if(!is.null(how)) { # {{{
    dat$sample = as.factor(sapply(pData(x)[,how],rep,times=probes)) # }}}
  } else { # {{{
    dat$sample = dat$chip
  } # }}}

  if(by.channel) { # {{{
    if('COLOR_CHANNEL' %in% fvarLabels(x)) {
      dat$channel = rep(fData(x)$COLOR_CHANNEL, chips)
    } else { 
      if( annotation(x) == 'IlluminaHumanMethylation450k' ) {
        require(IlluminaHumanMethylation450k.db)
        dat$channel = unlist(mget(dat$probe, 
                                  IlluminaHumanMethylation450kCOLOR_CHANNEL))
      } else if( annotation(x) == 'IlluminaHumanMethylation27k' ) { 
        require(IlluminaHumanMethylation27k.db)
        dat$channel = unlist(mget(dat$probe, 
                                  IlluminaHumanMethylation27kCOLOR_CHANNEL))
      }
    }
    dat$channel = as.character(dat$channel)
    dat$channel[ which(dat$channel == 'Red') ] = 'Cy5 (Infinium I)'
    dat$channel[ which(dat$channel == 'Grn') ] = 'Cy3 (Infinium I)'
    dat$channel[ which(dat$channel == 'Both') ] = 'Cy3/Cy5 (Infinium II)'
    dat$channel = as.factor(dat$channel) # }}}
  } else if( by.design ) { # {{{
    if('DESIGN' %in% fvarLabels(x)) {
      dat$design = c(rep(fData(x)$DESIGN, chips))
    } else { 
      if( annotation(x) == 'IlluminaHumanMethylation450k' ) {
        require(IlluminaHumanMethylation450k.db)
        dat$design = unlist(mget(dat$probe, IlluminaHumanMethylation450kDESIGN))
      } else if( annotation(x) == 'IlluminaHumanMethylation27k' ) { 
        dat$design = 'I' 
      }
    }
    dat$design = as.factor(paste('Infinium', as.character(dat$design)))
  } # }}}

  if(plot.type == 'histogram') { # {{{
    dat.hist = ggplot(dat, aes(x=assay, fill=sample, group=chip)) +
               geom_histogram(binwidth=.01, alpha=I(0.5),
                              position=position_identity()) # }}}
  } else if( plot.type == 'density' ) { # {{{
    dat.hist = ggplot(dat, aes(x=assay, colour=sample, group=chip)) +
               geom_density() 
  } # }}}

  if(by.channel) { # {{{
    dat.hist = dat.hist + facet_grid( . ~ channel ) # }}}
  } else if(by.design) { # {{{
    dat.hist = dat.hist + facet_grid( . ~ design )
  } # }}}
  if(!is.null(title)) { # {{{
    dat.hist = dat.hist + opts(title=title)
  } # }}}

  if(what=='betas') { # {{{
    return( dat.hist + 
            xlab('estimated percent methylation') +
            scale_x_continuous(formatter='percent') )  # }}}
  } else { # {{{
    return( dat.hist + xlab(what) + theme_bw() )
  } # }}}

} # }}}

## FIXME: export these qc.probe.plot() cleanly and replace qcplot()?
#if(require('ggplot2')) { # {{{
#  setMethod('qcplot', signature(object='MethyLumiQC'), # {{{
#    function(object,controlType,...) qc.probe.plot(object,controlType,...))#}}}
  # setMethod('qcplot', signature(object='MethyLumiSet'), # {{{
  #  function(object, controlType, ...) {
  #    if( length(annotation(object) == 0) ) qcplot(object, controlType)
  #    if( sampleNames(object@QC) != sampleNames(object) ) {
  #      sampleNames(object@QC) <- sampleNames(object) 
  #    }
  #    qc.probe.plot(QCdata(object), controlType, ...) 
  #  }) # }}}
  # setMethod('qcplot', signature(object='MethyLumiM'), # {{{
  #  function(object,controlType,...) {
  #    qc.probe.plot(QCdata(object), controlType, ...) 
  #}) # }}}
#} # }}}
