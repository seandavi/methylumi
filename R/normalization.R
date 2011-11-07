# utility functions
t.submit <- function() as.character(Sys.time())
t.finish <- function() as.character(format(Sys.time(), "%H:%M:%S"))

normalizeDesignIIViaSQN <- function(x, N.mix=3, weight=0.5, lg2=F){ # {{{ Wu's

  require(SQN)
  if( annotation(x) == 'IlluminaHumanMethylation27k' ) 
    stop('Subset quantile normalization should not be used on 27k arrays')
  if( is.null(QCdata(x)) ) stop('You do not seem to have any control probes.')

  history.submitted <- as.character(Sys.time())
  if('DESIGN' %in% fvarLabels(x)) {
    dII.probes = featureNames(x)[which(fData(x)$DESIGN == 'II')]
  } else {
    require(IlluminaHumanMethylation450k.db)
    fData(x)$DESIGN = mget(featureNames(x), IlluminaHumanMethylation450kDESIGN)
    dII.probes = featureNames(x)[which(fData(x)$DESIGN == 'II')]
  }
  normprobes.Cy3 = rownames(normctls(x,'Cy3'))
  normprobes.Cy5 = rownames(normctls(x,'Cy5'))
  if(length(normprobes.Cy3) != length(normprobes.Cy5)){
    stop("The numbers of normalization probes in each channel should be equal")
  } else {
    normnumber = length(normprobes.Cy5)
  }

  # could parallelize this easily with foreach()
  for( subject in sampleNames(x)) {
    Cy3.col = c(methylated(x)[dII.probes,subject],normctls(x,'Cy3')[,subject])
    Cy5.col = c(unmethylated(x)[dII.probes,subject],normctls(x,'Cy5')[,subject])
    ctrl.id = (length(dII.probes)+1):(length(dII.probes)+normnumber)
    bound = cbind(Cy3.col, Cy5.col)
    if(lg2) bound = log2(bound)
    normed = SQN(bound, N.mix, ctrl.id=ctrl.id, weight)
    if(lg2) normed = 2**normed
    methylated(x)[dII.probes, subject] = normed[1:length(dII.probes),1]
    Cy3(x@QC)[normprobes.Cy3, subject] = normed[ctrl.id,1]
    unmethylated(x)[dII.probes, subject] = normed[1:length(dII.probes),2]
    Cy5(x@QC)[normprobes.Cy5, subject] = normed[ctrl.id,2]
  }

  history.command <- "Applied subset quantile normalization to design II probes"
  history.finished <- t.finish()
  x@history<- rbind(x@history,
                    data.frame(submitted=history.submitted,
                               finished=history.finished,
                               command=history.command))
  return(x)

} # }}}

normalizeSamplesViaSQN <- function(x, N.mix=15, weight=0.9, lg2=T) { # {{{ Wu's

  require(SQN)
  if( annotation(x) == 'IlluminaHumanMethylation27k' ) 
    stop('Subset quantile normalization should not be used on 27k arrays')
  if( is.null(QCdata(x)) ) stop('You do not seem to have any control probes.')

  if(!('COLOR_CHANNEL' %in% fvarLabels(x))) { # {{{
    fData(x)$COLOR_CHANNEL = mget(featureNames(x),
                                  IlluminaHumanMethylation450kCOLORCHANNEL)
  } # }}}
  d2.probes <- which(fData(x)[['COLOR_CHANNEL']]=='Both') 
  cy3.probes <- which(fData(x)[['COLOR_CHANNEL']]=='Grn')
  cy5.probes <- which(fData(x)[['COLOR_CHANNEL']]=='Red')
  probes = list(Cy3=cy3.probes, Cy5=cy5.probes)
  normprobes.Cy3 = rownames(normctls(x,'Cy3'))
  normprobes.Cy5 = rownames(normctls(x,'Cy5'))
  if(length(normprobes.Cy3) != length(normprobes.Cy5)){
    stop("The numbers of normalization probes in each channel should be equal")
  } else normnumber = length(normprobes.Cy5)
  history.submitted <- as.character(Sys.time())

  channels = c('Cy3','Cy5')
  names(channels) = channels
  proberows = lapply(channels, function(channel) {
    c(D2=length(d2.probes), 
      M=length(probes[[channel]]), 
      U=length(probes[[channel]]),
      C=normnumber)
  })
  last = lapply(proberows, cumsum)
  first = lapply(channels, function(ch) (last[[ch]]-proberows[[ch]]) + 1)
  names(first) = names(last) = channels
  for(channel in channels) {
    if( channel == 'Cy3' ) D2.intensities = methylated(x)[d2.probes,]
    if( channel == 'Cy5' ) D2.intensities = unmethylated(x)[d2.probes,]
    bound = rbind(D2.intensities, 
                  methylated(x)[probes[[channel]],],
                  unmethylated(x)[probes[[channel]],],
                  normctls(x,channel))
    ctrl.id = first[[channel]][['C']]:last[[channel]][['C']]
    if(lg2) bound = log2(bound)
    normed = SQN(bound, N.mix, ctrl.id=ctrl.id, weight)
    if(lg2) normed = 2**normed
    d2.rows = first[[channel]][['D2']]:last[[channel]][['D2']]
    if( channel == 'Cy3' ) methylated(x)[d2.probes,] = normed[d2.rows,]
    if( channel == 'Cy5' ) unmethylated(x)[d2.probes,] = normed[d2.rows,]
    d1.M.rows = first[[channel]][['M']]:last[[channel]][['M']]
    methylated(x)[probes[[channel]], ] = normed[d1.M.rows, ]
    d1.U.rows = first[[channel]][['U']]:last[[channel]][['U']]
    unmethylated(x)[probes[[channel]], ] = normed[d1.U.rows, ]
    if(channel == 'Cy3') Cy3(x@QC)[normprobes.Cy3, ] = normed[ctrl.id,]
    if(channel == 'Cy5') Cy5(x@QC)[normprobes.Cy5, ] = normed[ctrl.id,]
  }
  history.command <- "Applied subset quantile normalization across samples."
  history.finished <- t.finish()
  x@history<- rbind(x@history,
                    data.frame(submitted=history.submitted,
                               finished=history.finished,
                               command=history.command))
  return(x)

} # }}}

normalizeViaRUV2 <- function(x) { # {{{ Terry Speed's factor extractor
  stop('Terry Speed\'s method has not yet been patched into methylumi')
} # }}}

normalizeViaControls <- function(x, reference=1) { # {{{ from minfi/KDH

  if(is.null(x@QC)) stop('Cannot normalize against controls without controls!')
  else history.submitted <- as.character(Sys.time())

  # this is easier in methylumi
  controls <- normctls(x)
  Grn.avg <- colMeans(controls$Cy3)
  Red.avg <- colMeans(controls$Cy5)

  # this is about the same 
  ref <- (Grn.avg + Red.avg)[reference]/2
  if(is.na(ref)) stop("'reference' refers to an array that is not present")
  Grn.factor <- ref/Grn.avg
  Red.factor <- ref/Red.avg

  # this is harder in methylumi
  Grn <- list( M1=methylated(x)[ which(fData(x)$COLOR_CHANNEL=='Grn'), ],
               U1=unmethylated(x)[ which(fData(x)$COLOR_CHANNEL=='Grn'), ],
               M2=methylated(x)[ which(fData(x)$COLOR_CHANNEL=='Both'), ] )
  Grn <- lapply(Grn, function(y) sweep(y, 2, FUN="*", Grn.factor))
  methylated(x)[ which(fData(x)$COLOR_CHANNEL=='Grn'), ] <- Grn$M1
  unmethylated(x)[ which(fData(x)$COLOR_CHANNEL=='Grn'), ] <- Grn$U1
  methylated(x)[ which(fData(x)$COLOR_CHANNEL=='Both'), ] <- Grn$M2

  Red <- list( M1=methylated(x)[ which(fData(x)$COLOR_CHANNEL=='Red'), ],
               U1=unmethylated(x)[ which(fData(x)$COLOR_CHANNEL=='Red'), ],
               U2=unmethylated(x)[ which(fData(x)$COLOR_CHANNEL=='Both'), ] )
  Red <- lapply(Red, function(y) sweep(y, 2, FUN="*", Red.factor))
  methylated(x)[ which(fData(x)$COLOR_CHANNEL=='Red'), ] <- Red$M1
  unmethylated(x)[ which(fData(x)$COLOR_CHANNEL=='Red'), ] <- Red$U1
  unmethylated(x)[ which(fData(x)$COLOR_CHANNEL=='Both'), ] <- Red$U2

  # now do the same to the controls (for plotting purposes)
  ctls = list(Cy3=methylated(x@QC), Cy5=unmethylated(x@QC))
  assayDataElement(x@QC,'methylated') <- sweep(ctls$Cy3, 2, FUN='*', Grn.factor)
  assayDataElement(x@QC,'unmethylated') <- sweep(ctls$Cy5,2,FUN='*', Red.factor)

  # and add an entry to the transaction log for this preprocessing step.
  history.command <- deparse(match.call())
  history.finished <- t.finish()
  x@history<- rbind(x@history,
                    data.frame(submitted=history.submitted,
                               finished=history.finished,
                               command=history.command))
  return(x)

} # }}}
