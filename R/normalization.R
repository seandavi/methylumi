# utility functions
t.submit <- function() as.character(Sys.time())
t.finish <- function() as.character(format(Sys.time(), "%H:%M:%S"))

normalizeViaControls <- function(x, reference=NULL) { # {{{ originally by Kasper

  if(is.null(x@QC)) stop('Cannot normalize against controls without controls!')
  else history.submitted <- as.character(Sys.time())

  # this is easier in methylumi
  controls <- normctls(x)
  Grn.avg <- colMeans(controls$Cy3)
  Red.avg <- colMeans(controls$Cy5)
  R.G.ratio = Red.avg/Grn.avg
  if(is.null(reference)) reference = which.min( abs(R.G.ratio-1) )
  message(paste('Using sample number', reference, 'as reference level...'))

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

  # now fix the beta values, and set the appropriate ones to NA
  betas(x) <- methylated(x) / total.intensity(x)
  is.na(betas(x)) <- (pvals(x) > 0.05)

  # and add an entry to the transaction log for this preprocessing step.
  history.command <- deparse(match.call())
  history.finished <- t.finish()
  x@history<- rbind(x@history,
                    data.frame(submitted=history.submitted,
                               finished=history.finished,
                               command=history.command))
  return(x)

} # }}}
