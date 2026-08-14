setClassUnion("methylData", c('MethyLumiSet','MethyLumiM'))
setClassUnion("ND", c('character','missing'))

# fix for QC plots in methylumi (hard to read on Infinium arrays) using ggplot2
qc.probe.plot <- function(obj,controltype="negnorm",log2=TRUE,by.type=FALSE,...){ # {{{
  log2_trans = log_trans(base=2)
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
    rows <- grep('(Negative|Norm)', fData(qc)$Type, ignore.case=TRUE)
  } else { 
    rows <- grep(paste('^',controltype,sep=''),fData(qc)$Type,ignore.case=TRUE)
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
    p <- ggplot2::ggplot(qc, ggplot2::aes(x = variable, y = value, colour = type,
                                          fill = type, group = grouping)) +
                        ggplot2::geom_boxplot() +
                        ggplot2::labs(title = a.title,
                                      x = "Sample", y = "Intensities") +
                        coord_flip()
    if (log2) p <- p + scale_x_continuous(trans='log2', limits=c(2**4, 2**16))
    else p <- p + scale_x_continuous(limits=c(2**4, 2**16))
  } else {
    p <- ggplot2::ggplot(qc, ggplot2::aes(x = value, y = variable,
                                          colour = type, shape = type)) +
                        ## geometry is 'jitter' for negnorm controls, else 'point'
                        (if (geometry == 'jitter') ggplot2::geom_jitter()
                         else ggplot2::geom_point()) +
                        ggplot2::labs(title = a.title,
                                      x = "Intensities", y = "Sample")
    p <- p + scale_y_discrete( limits=rev(sampleNames(obj)) ) # more readable
    if (log2) p <- p + scale_x_continuous(trans='log2', limits=c(2**4, 2**16))
    else p <- p + scale_x_continuous(limits=c(2**4, 2**16))
  }
  if( by.type ) {
    p <- p + facet_grid( type ~ channel)
  } else { 
    p <- p + facet_grid( . ~ channel)
  }
  if( tolower(controltype) == 'negnorm' ) {
    p <- p + scale_colour_manual( values=colour.settings )
    p <- p + scale_shape_manual( values=shape.settings )
  }
  p <- p + theme_bw()
  return( p )
} # }}}

methylumi.diagnostics <- function (x, onlybg=FALSE) { # {{{
  x.qc <- controlData(x)
  if (!is.null(x.qc)) { # {{{ realistically, use OOB
      bg <- list(red = log2(negctls(x.qc, "Cy5")), 
                 green = log2(negctls(x.qc, "Cy3")))
      message("(Using negative controls for dashed vertical background line)")
  } # }}}
  ## annotation() is character(0) on objects that never had a platform set
  ## (e.g. mldat, read from a Sentrix CSV), so compare with identical(), not ==.
  if(!('COLOR_CHANNEL' %in% fvarLabels(x))) { # {{{
    if(identical(annotation(x), 'IlluminaHumanMethylation27k')) {
      fData(x)$COLOR_CHANNEL = mget(featureNames(x),
                                    IlluminaHumanMethylation27kCOLORCHANNEL)
    } else if(identical(annotation(x), 'IlluminaHumanMethylation450k')) {
      fData(x)$COLOR_CHANNEL = mget(featureNames(x),
                                    IlluminaHumanMethylation450kCOLORCHANNEL)
    }
  } # }}}
  assays <- c("exprs", "methylated", "unmethylated")
  if(onlybg) assays <- c("methylated", "unmethylated")
  is.450k = identical(annotation(x), 'IlluminaHumanMethylation450k')
  # {{{ one panel column per probe set; the red/Grn swap below is deliberate
  colorchannel <- fData(x)[['COLOR_CHANNEL']]
  if(is.null(colorchannel)) {
    # ponytail: no per-probe colour channel (GoldenGate data, unannotated 27k
    # CSV); plot every probe in a single column rather than erroring out.
    dye <- list(all='all')
    probesets <- list(all=seq_len(nrow(x)))
  } else {
    dye <- list('red'='Cy5','green'='Cy3')
    probesets <- list(red=which(colorchannel=='Grn'),
                      green=which(colorchannel=='Red'))
    if(is.450k) {
      dye[['both']] = 'Design II'
      probesets[['both']] = which(colorchannel=='Both')
    }
  } # }}}
  par(mfrow = c(length(assays), length(dye)))
  for(assay in assays) {
    for(channel in names(dye)) {
      probes <- probesets[[channel]]
      chcolor <- switch(channel,
                        both = ifelse(assay == 'exprs', 'blue',
                                      ifelse(assay == 'methylated', 'green',
                                             'red')),
                        all = 'blue',
                        channel)
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
        for (i in seq_along(dens)) lines(dens[[i]], lty=2, col=chcolor)
        title(paste(xlab, ":", dye[[channel]], "probes")) # }}}
      } else { # {{{ methylated/unmethylated
        dat <- log2(assayDataElement(x, assay)[probes,])
        ## stats::plot.density is a method, not an exported function: call plot()
        plot(density(na.omit(dat)), col = "white", xlim = c(0, 16),
             ylim = c(0, 0.8), lwd = 1, xlab = "log2(Intensity)",
             ylab = "Proportion", main = "", lty = 1)
        for(i in 1:dim(dat)[2]) {
          lines(density(na.omit(dat[, i])), col=chcolor, lty=1)
        }
        if (!is.null(x.qc)) {
          if(channel=='both') bgch = bg[[ifelse(assay=='methylated','green','red')]]
          else if(channel=='all') bgch = do.call(rbind, bg)
          else bgch = bg[[channel]]
          for(bgmean in colMeans(bgch, na.rm=TRUE)) {
            abline(v=bgmean, col=paste("dark",chcolor,sep=""), lty=3, lwd=1)
          }
        }
        title(paste(assay, "intensities:", dye[[channel]]))
      } # }}}
    }
  }
} # }}}
