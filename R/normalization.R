# utility functions
t.submit <- function() as.character(Sys.time())
t.finish <- function() as.character(format(Sys.time(), "%H:%M:%S"))

normalizeViaSQN <- function(x, N.mix=3, weight=0.5){ # {{{ Wu & Aryee 

  require(SQN)
  if( annotation(x) == 'IlluminaHumanMethylation27k' ) 
    stop('Subset quantile normalization should not be used on 27k arrays')
  stopifnot( 'methylated.OOB' %in% assayDataNames(x) & 
             'unmethylated.OOB' %in% assayDataNames(x) )
  history.submitted <- as.character(Sys.time())
  if(!('DESIGN' %in% fvarLabels(x))) {
    require(IlluminaHumanMethylation450k.db)
    fData(x)$DESIGN = mget(featureNames(x), IlluminaHumanMethylation450kDESIGN)
  }
  dII.probes = featureNames(x)[which(fData(x)$DESIGN == 'II')]
  dII = list(Cy3=methylated(x)[dII.probes,], Cy5=unmethylated(x)[dII.probes,])
  dII.rows = lapply(dII.probes, seq_along)
  IB = intensities.IB(x)
  OB = intensities.OOB(x)
  Cy3.intensities = rbind( dII$Cy3, 
                           IB$Cy3, 
                           OB$Cy3 )
  Cy5.intensities = rbind( dII$Cy5, 
                           IB$Cy5[seq_along(rownames(IB$Cy3)),],
                           OB$Cy5[seq_along(rownames(OB$Cy3)),] )
  if(!is.null(QCdata(x))) {
    Cy3.intensities = rbind( Cy3.intensities, normctls(x)$Cy3 )
    Cy5.intensities = rbind( Cy5.intensities, normctls(x)$Cy5 )
  }
  ctrl.id = list(Cy3=setdiff(seq_along(rownames(Cy3.intensities)),dII.rows$Cy3),
                 Cy5=setdiff(seq_along(rownames(Cy5.intensities)),dII.rows$Cy5))

  # subset quantile normalize each array with itself (green vs. red)
  for( subject in sampleNames(x)) {
    
    stop('need to finish within-array SQN')
    ctrl.id = (length(dII.probes)+1):(length(dII.probes)+normnumber)
    bound = cbind(Cy3.col, Cy5.col)
    withins = SQN(bound, N.mix, ctrl.id=ctrl.id, weight)
    methylated(x)[dII.probes, subject] = withins[1:length(dII.probes),1]
    Cy3(x@QC)[normprobes.Cy3, subject] = withins[ctrl.id,1]
    unmethylated(x)[dII.probes, subject] = withins[1:length(dII.probes),2]
    Cy5(x@QC)[normprobes.Cy5, subject] = withins[ctrl.id,2]
  }

  # now normalize all the arrays against each other in each channel
  for( channel in names(dII) ) {

    stop('need to finish between-array SQN')
    betweens = SQN(bound, N.mix, ctrl.id=ctrl.id, weight)
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

normalizeViaControls <- function(x, reference=1) { # {{{ from Kasper  

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

normalizeViaRUV2 <- function(x) { # {{{ Terry Speed's factor extractor
  stop('Terry Speed\'s method has not yet been patched into methylumi')
} # }}}
