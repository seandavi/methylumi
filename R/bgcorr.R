# utility functions
t.submit <- function() as.character(Sys.time())
t.finish <- function() as.character(format(Sys.time(), "%H:%M:%S"))

# generic dispatcher for background correction of Infinium methylation arrays
methylumi.bgcorr<-function(x, method='noob', offset=15, controls=NULL, correct=T, parallel=F, ...) { # {{{

  allelic = FALSE

  # GoldenGate: not supported (it could be, perhaps, but why?) 
  stopifnot(grepl('IlluminaHumanMethylation', annotation(x)))

  if(any(methylated(x)<=0)) { # {{{
    methylated(x)[which(methylated(x)==0)] <- 1
  } # }}}
  if(any(unmethylated(x)<=0)) { # {{{
    unmethylated(x)[which(unmethylated(x)==0)] <- 1
  } # }}}

  # 'noob' and 'goob' are just shorthand for 'normexp, OOB', and 'gamma, OOB'
  if(tolower(method) %in% c('noob','goob','mode')) controls = intensities.OOB(x)
  if(tolower(method) == 'lumi') controls = intensities.M(x)
  if(tolower(method) == 'noob') method = 'normexp'
  if(tolower(method) == 'goob') method = 'gamma'

  # controls must be a list(Cy3=matrix, Cy5=matrix) if supplied
  if(!is.null(controls) &&  # {{{
     !( 'Cy5' %in% names(controls) & 'Cy3' %in% names(controls))) { 
    stop('If you supply controls, they must be a list(Cy3=matrix, Cy5=matrix)')
  } # }}}
  if(is.null(controls)) { # {{{
    if(is(x, 'MethyLumiSet')) { 
      stopifnot('QC' %in% slotNames(x))
      stopifnot(!is.null(QCdata(x)))
      controls = negnorm(controlData(x))
    } else if(is(x, 'MethyLumiM')) {
      stopifnot('controlData' %in% slotNames(x))
      stopifnot(!is.null(controlData(x)))
      controls = negnorm(controlData(x))
    }
  } # }}}

  history.submitted = as.character(Sys.time())
  if(!('COLOR_CHANNEL' %in% fvarLabels(x))) { # {{{
    if(annotation(x) == 'IlluminaHumanMethylation27k') { 
      fData(x)$COLOR_CHANNEL = mget(featureNames(x),
                                   IlluminaHumanMethylation27kCOLORCHANNEL)
    } else if(annotation(x) == 'IlluminaHumanMethylation450k') { 
      fData(x)$COLOR_CHANNEL = mget(featureNames(x),
                                   IlluminaHumanMethylation450kCOLORCHANNEL)
    }
  } # }}}
  cy3.probes <- which(fData(x)[['COLOR_CHANNEL']]=='Grn')
  cy5.probes <- which(fData(x)[['COLOR_CHANNEL']]=='Red')
  d2.probes <- which(fData(x)[['COLOR_CHANNEL']]=='Both') 

  dat = list( Cy3=list( M=as.matrix(methylated(x)[cy3.probes,]), # {{{
                        U=as.matrix(unmethylated(x)[cy3.probes,])), # }}}
              Cy5=list( M=as.matrix(methylated(x)[cy5.probes,]), # {{{
                        U=as.matrix(unmethylated(x)[cy5.probes,]))) # }}}
  if(annotation(x) == 'IlluminaHumanMethylation450k') {  # {{{
    #if(allelic) {
      #dat[['D2']] = list( M=as.matrix(methylated(x)[d2.probes, ]),
                          #U=as.matrix(unmethylated(x)[d2.probes, ]) )
    #} else { 
      dat[['Cy3']][['D2']] = as.matrix(methylated(x)[d2.probes, ])
      dat[['Cy5']][['D2']] = as.matrix(unmethylated(x)[d2.probes, ])
    #}
  } # }}}
  # if(!allelic) { # ranges in corrected data, to update existing matrices {{{
  rows=lapply(dat,function(ch) sapply(names(ch),function(nch)nrow(ch[[nch]])))
  last = lapply(rows, cumsum)
  first = lapply(names(last), function(nch) last[[nch]] - rows[[nch]] + 1)
  names(first) = names(last)
  # }}}

  #if(allelic) { # {{{
    #stop("Don't use this (allelic=TRUE) for the time being.")
    #channels = list(Cy5='Cy5',Cy3='Cy3')
    #if(annotation(x) == 'IlluminaHumanMethylation450k') channels[['D2']] = 'D2'
    #estimates = lapply(channels, function(nch) { 
      #alleles = list(M='M',U='U')
      #lapply(alleles, function(al) {
        #xf = dat[[nch]][[al]]
        #if(!is.null(dat[[nch]][['D2']])) xf = rbind(xf, dat[[nch]][['D2']])
        #xs = get.xs(xf, controls[[nch]][[al]], method=method,
                    #offset=offset, correct=correct)
        #names(xs[['params']]) = paste(names(xs[['params']]), nch, sep='.')
        #names(xs[['meta']]) = paste(names(xs[['meta']]), nch, sep='.')
        #xs
      #})
    #}) # }}}
  #} else { 

  if(parallel) {
    estimates = .mclapply(names(dat), function(nch) { 
      xf = rbind(dat[[nch]][['M']], dat[[nch]][['U']])
      if(!is.null(dat[[nch]][['D2']])) xf = rbind(xf, dat[[nch]][['D2']])
      xs = get.xs(xf, controls[[nch]], method=method,
                  offset=offset, correct=correct, parallel=TRUE)
      names(xs[['params']]) = paste(names(xs[['params']]), nch, sep='.')
      names(xs[['meta']]) = paste(names(xs[['meta']]), nch, sep='.')
      xs
    }) 
  } else { 
    estimates = lapply(names(dat), function(nch) { 
      xf = rbind(dat[[nch]][['M']], dat[[nch]][['U']])
      if(!is.null(dat[[nch]][['D2']])) xf = rbind(xf, dat[[nch]][['D2']])
      xs = get.xs(xf, controls[[nch]], method=method,
                  offset=offset, correct=correct, parallel=FALSE)
      names(xs[['params']]) = paste(names(xs[['params']]), nch, sep='.')
      names(xs[['meta']]) = paste(names(xs[['meta']]), nch, sep='.')
      xs
    }) 
  }
  names(estimates) = names(dat)
                        
  ## fix any Illumina control probes using our chosen background estimates
  if(!is.null(controlData(x)) && !allelic) { # {{{
    ctrl = controlData(x)
    xcs = lapply(names(dat), function(nch) {
      xcf = as.matrix(intensitiesByChannel(ctrl, nch))
      get.xcs(xcf, method, params=estimates[[nch]][['params']], correct=correct)
    })
    names(xcs) = names(dat)
    Cy3(ctrl) = xcs[['Cy3']]
    Cy5(ctrl) = xcs[['Cy5']]
    controlData(x) = ctrl
  } # }}}

  #if( allelic ) { # {{{ 
    #stop("Don't use this for the time being.")
    ## {{{
    #methylated(x)[cy3.probes,] <- estimates[['Cy3']][['M']][['xs']]
    #unmethylated(x)[cy3.probes,] <- estimates[['Cy3']][['U']][['xs']]
    #methylated(x)[cy5.probes,] <- estimates[['Cy5']][['M']][['xs']]
    #unmethylated(x)[cy5.probes,] <- estimates[['Cy5']][['U']][['xs']]
    #if(!is.null(estimates[['D2']])) { # 450k design II
      #message('Consider using the random loci for nonspecific intensity here')
      #methylated(x)[d2.probes,] <- estimates[['D2']][['M']][['xs']]
      #unmethylated(x)[d2.probes,] <- estimates[['D2']][['U']][['xs']]
    #} # }}}
  #} else {  # }}}
  cy3.M = first[['Cy3']][['M']]:last[['Cy3']][['M']]
  methylated(x)[cy3.probes,] <- estimates[['Cy3']][['xs']][cy3.M,]
  cy3.U = first[['Cy3']][['U']]:last[['Cy3']][['U']]
  unmethylated(x)[cy3.probes,] <- estimates[['Cy3']][['xs']][cy3.U,]
  cy5.M = first[['Cy5']][['M']]:last[['Cy5']][['M']]
  methylated(x)[cy5.probes,] <- estimates[['Cy5']][['xs']][cy5.M,]
  cy5.U = first[['Cy5']][['U']]:last[['Cy5']][['U']]
  unmethylated(x)[cy5.probes,] <- estimates[['Cy5']][['xs']][cy5.U,]
  if(!is.null(dat[['Cy3']][['D2']])) { # 450k design II
    d2.M = first[['Cy3']][['D2']]:last[['Cy3']][['D2']]
    d2.U = first[['Cy5']][['D2']]:last[['Cy5']][['D2']]
    methylated(x)[d2.probes,] <- estimates[['Cy3']][['xs']][d2.M,]
    unmethylated(x)[d2.probes,] <- estimates[['Cy5']][['xs']][d2.U,]
  }
  #} 

  for(ch in names(estimates)) { # {{{
    chnames = names(estimates[[ch]][['params']])
    for(nm in chnames) pData(x)[,nm] = estimates[[ch]][['params']][[nm]]
    varMetadata(x)[chnames,] = paste(ch, estimates[[ch]][['meta']]) 
  } # }}}

  pcutoff <- pval.detect(x)
  betas(x) <- pmax(methylated(x), 1) / pmax(total.intensity(x), 2)
  betas(x)[ which(pvals(x) > pcutoff) ] <- NA 

  history.command <- deparse(match.call())
  history.finished <- t.finish()
  x@history<- rbind(x@history,
                    data.frame(submitted=history.submitted,
                               finished=history.finished,
                               command=history.command))
  return(x)
} # }}}

# dispatchers for arbitrary background correction methods
get.xs <- function(xf, controls, method, offset=50, robust=T, correct=T, parallel=F) { # {{{
  if(method == 'normexp') return(normexp.get.xs(xf, controls, offset, robust))
  if(method == 'median') return(median.get.xs(xf, controls, offset))
  if(method == 'illumina') return(illumina.get.xs(xf, controls, offset))
  if(method == 'gamma') return(gamma.get.xs(xf, controls, offset, correct, 
                                            parallel=parallel))
  if(method == 'mode') return(gammaM.get.xs(xf, controls, offset, correct))
  if(method == 'lumi') return(lumi.get.xs(xf, controls, offset))
  else stop(paste('Method',method,'has not been added to get.xs() yet'))
}  # }}}
normexp.get.xs <- function(xf, controls, offset=50, robust=T, ...){#{{{
  cat("Background mean & SD estimated from", nrow(controls), "probes\n")
  if(robust){ # {{{ slower
    require(MASS)
    mu <- sigma <- alpha <- rep(NA, ncol(xf))
    for( i in 1:ncol(xf) ) {
      ests <- huber(controls[, i])
      mu[i] <- ests$mu
      sigma[i] <- ests$s
      alpha[i] <- max(huber(xf[, i])$mu - mu[i], 10)
    } # }}}
  } else { # {{{ faster 
    require(matrixStats) 
    mu = colMeans(controls, na.rm=T)
    sigma = colSds(xf, na.rm=T)
    alpha <- pmax((colMeans(xf, na.rm=T) - mu), 10)
  } # }}}
  pars = data.frame(mu=mu, lsigma=log(sigma), lalpha=log(alpha))
  for(i in 1:ncol(xf)) xf[,i] <- normexp.signal(as.list(pars[i,]), xf[,i])
  return(list(xs=xf+offset, 
              params=data.frame(mu=mu, sigma=sigma, alpha=alpha, offset=offset),
              meta=c('background mean','background SD','signal mean','offset')))
} # }}}
median.get.xs <- function(xf, controls, offset=50, robust=T, ...){#{{{
  cat("Background median estimated from", nrow(controls), "probes\n")
  require(matrixStats)
  bg = colMedians(controls, na.rm=T)
  for(i in 1:ncol(xf)) xf[,i] <- pmax( xf[,i] - bg[i], 1 )
  return(list(xs=xf+offset, 
              params=data.frame(median=bg, offset=offset),
              meta=c('background median','offset')))
} # }}}
illumina.get.xs <- function(xf, controls, offset=50, robust=T, ...){#{{{
  require(matrixStats)
  bg = colQuantiles(controls, 0.05)
  for(i in 1:ncol(xf)) xf[,i] <- pmax( xf[,i] - bg[i], 1 )
  return(list(xs=xf+offset, 
              params=data.frame(bg=bg, offset=offset),
              meta=c('background fifth percentile','offset')))
} # }}}
gamma.get.xs <- function(xf,controls,offset=50,correct=T,parallel=F,...){#{{{

  #require(rGammaGamma)
  bg = sapply(1:ncol(xf), function(i) gamma.mle(controls[,i]))
  if(correct) { # {{{
    bgmu = colMeans(controls, na.rm=T)
    fg = sapply(1:ncol(xf), function(i) gamma.mle(pmax(xf[,i]-bgmu[i], 1))) #}}}
  } else { # {{{
    fg = sapply(1:ncol(xf), function(i) gamma.mle(xf[,i]))
  } # }}}
  params = cbind(t(fg), t(bg))
  colnames(params) = c('gamma','alpha','delta','beta')
  meta = c('signal shape','signal scale','background shape','background scale')
  names(meta) = c('gamma','alpha','delta','beta')

  if( parallel == TRUE ) {
    xs = data.matrix(as.data.frame(.mclapply(1:ncol(xf), function(i) {
           gamma.integral(xf[,i], params[i,], offset=offset)
    })))
  } else { 
    cat('Estimating xs serially (probably not what you want)...', "\n")
    xs = data.matrix(as.data.frame(lapply(1:ncol(xf), function(i) {
           gamma.integral(xf[,i], params[i,], offset=offset)
    })))
  }
  params = cbind(params, c(rep(offset, nrow(params))))
  colnames(params) = c('gamma','alpha','delta','beta','offset')
  meta['offset'] = 'offset'
  return(list(xs=xs, params=as.data.frame(params), meta=meta))
  
} # }}}
gammaM.get.xs <- function(xf,controls,offset=15,correct=T,parallel=F,...){#{{{

  #require(rGammaGamma)
  bg = sapply(1:ncol(xf), function(i) gamma.mode(gamma.mle(controls[,i])))
  xs = sapply(1:ncol(xf), function(s) pmax(xf[,s] - bg[s], offset))
  params = data.frame(mode=bg)
  params$offset = offset
  meta = c('background mode','offset')
  names(meta) = c('mode','offset')
  cat("Background mode estimated from", nrow(controls), "probes\n")
  return(list(xs=xs, params=params, meta=meta))
  
} # }}}
lumi.get.xs <- function(xf, controls, offset=50, robust=T, nbin=1000,...){#{{{

  # from lumi's estimateBG() function
  bg <- apply(controls, 2, function(x) {
    hh.x <- hist(x, nbin, plot = FALSE)
    Th <- hh.x$breaks[which.max(hh.x$counts) + 1] * 2
    dd.x <- density(x[x < Th], na.rm = TRUE)
    bg.x <- dd.x$x[which.max(dd.x$y)]
  })
  cat("Background mode estimated from", nrow(controls), "probes\n")

  for(i in 1:ncol(xf)) { # screams out to be parallelized
    xf[,i] <- pmax((xf[,i] - bg[i]), 1)
  }
  return(list(xs=xf+offset, 
              params=data.frame(mode=bg, offset=offset),
              meta=c('background mode','offset')))
} # }}}

# dispatchers for correcting controls w/o including them in estimates 
get.xcs <- function(xcf, method, params, robust=T, correct=T) { # {{{
  if(method == 'normexp') return(normexp.get.xcs(xcf, params))
  if(method == 'median') return(median.get.xcs(xcf, params))
  if(method == 'illumina') return(illumina.get.xcs(xcf, params,correct=correct))
  if(method == 'gamma') return(gamma.get.xcs(xcf, params))
  if(method == 'mode') return(gammaM.get.xcs(xcf, params))
  if(method == 'lumi') return(lumi.get.xcs(xcf, params))
  else stop(paste('Method',method,'has not been added to get.xcs() yet'))
}  # }}}
normexp.get.xcs <- function(xcf, params, ...){#{{{

  stopifnot(any(grepl('mu', names(params))))
  stopifnot(any(grepl('sigma', names(params))))
  stopifnot(any(grepl('alpha', names(params))))
  stopifnot(any(grepl('offset', names(params))))
  pars = data.frame(mu=params[[grep('mu', names(params), value=T)]],
                    sigma=log(params[[grep('sigma', names(params), value=T)]]),
                    alpha=log(params[[grep('alpha', names(params), value=T)]]) )
  for(i in 1:ncol(xcf)) xcf[,i] = normexp.signal( pars[i,], xcf[,i] )
  return( xcf + params[[grep('offset', names(params), value=T)]][1] )

} # }}}
median.get.xcs <- function(xcf, params, robust=T, ...){#{{{

  stopifnot(any(grepl('median', names(params))))
  stopifnot(any(grepl('offset', names(params))))
  xcs = xcf
  mu = params[[grep('median', names(params), value=T)]]
  offset = params[[grep('offset', names(params), value=T)]]
  for(i in 1:ncol(xcf)) xcs[,i] <- pmax(xcf[,i] - mu[i], 1)
  rm(xcf)
  xcs = xcs + offset
  return(xcs)

} # }}}
illumina.get.xcs <- function(xcf, params, robust=T, ...){#{{{

  stopifnot(any(grepl('bg', names(params))))
  stopifnot(any(grepl('offset', names(params))))
  bg = params[[grep('bg', names(params), value=T)]]
  offset = params[[grep('offset', names(params), value=T)]][1]
  for(i in 1:ncol(xcf)) xcf[,i] <- pmax(xcf[,i] - bg[i], 1)
  return(xcf+offset)

} # }}}
gamma.get.xcs <- function(xcf, params, robust=T, parallel=F,...){#{{{

  #require(rGammaGamma)
  offset = params[[grep('offset', names(params), value=T)]][1]
  params[[grep('offset', names(params), value=T)]] = NULL
  
  xcs = xcf
  for(i in 1:ncol(xcf)) {
    xcs[,i]=gamma.integral(xcf[,i],params=as.numeric(params[i,]),offset=offset)
  }
  # notice that the offset was added during the calculation of xcs|params.
  return(xcs)
  
} # }}}
gammaM.get.xcs <- function(xcf, params, robust=T, parallel=F, ...){ #{{{
  bgmode = params[[grep('mode', names(params), value=T)]]
  sapply(1:ncol(xcf), function(i) pmax((xcf[,i] - bgmode[i]), offset))
} # }}}
lumi.get.xcs <- function(xcf, params, robust=T, ...){#{{{

  stopifnot(any(grepl('mode', names(params))))
  stopifnot(any(grepl('offset', names(params))))
  mu = params[[grep('mode', names(params), value=T)]]
  offset = params[[grep('offset', names(params), value=T)]][1]
  for(i in 1:ncol(xcf)) xcf[,i] <- pmax(xcf[,i] - mu[i], 1)
  return(xcf+offset)

} # }}}

# normal-exponential deconvolution (conditional expectation of xs|xf; WEHI code)
normexp.signal <- function (par, x)  { # {{{
  par = unlist(par)
  mu <- par[1]
  sigma <- exp(par[2])
  sigma2 <- sigma * sigma
  alpha <- exp(par[3])
  if (alpha <= 0) 
    stop("alpha must be positive")
  if (sigma <= 0) 
    stop("sigma must be positive")
  mu.sf <- x - mu - sigma2/alpha
  signal <- mu.sf + sigma2 * exp(dnorm(0, mean = mu.sf, sd = sigma, log = TRUE)-
    pnorm(0, mean = mu.sf, sd = sigma, lower.tail = FALSE, log = TRUE))
  o <- !is.na(signal)
  if (any(signal[o] < 0)) {
    warning("Limit of numerical accuracy reached with very low intensity or very high background:\nsetting adjusted intensities to small value")
    signal[o] <- pmax(signal[o], 1e-06)
  }
  signal
} # }}}

# gamma deconvolution (conditional expectation of xs|xf; my code)
gamma.signal <- function (par, x)  { # {{{
  #require(rGammaGamma)
  par = as.numeric(par)
  gamma.integral(x, par, offset=0)
} # }}}
