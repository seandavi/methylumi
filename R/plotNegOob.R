plotNegOob <- function(x, log2=F) {

  require(ggplot2)
  samples = dim(x)[2]
  if(samples > 4) stop("Too many samples, choose a subset for a decent plot")
  par(mfrow=c(samples,2))
  max.neg = max( max(negctls(x, 'Cy5')), max(negctls(x, 'Cy3')) )
  xl = c(1, ifelse(log2, log2(max.neg+500), max.neg+500)) # x limits
  dat = list(Cy3.ctl=negctls(x, 'Cy3'), # {{{
             Cy3.oob.M=intensities.OOB.allelic(x, 'Cy3', 'methylated'),
             Cy3.oob.U=intensities.OOB.allelic(x, 'Cy3', 'unmethylated'),
             Cy5.ctls=negctls(x, 'Cy5'),
             Cy5.oob.M=intensities.OOB.allelic(x, 'Cy5', 'methylated'),
             Cy5.oob.U=intensities.OOB.allelic(x, 'Cy5', 'unmethylated')) # }}}
  merged = do.call(rbind, lapply(names(dat), function(d) { # {{{
    subdat = melt(dat[[d]])
    names(subdat)[2] = 'chip'
    names(subdat)[1] = 'probe'
    names(subdat)[3] = 'intensity'
    subdat$chip = gsub(' ', '.', subdat$chip)
    subdat$channel = substr(d, 1, 3)
    subdat$probes = substr(d, 5, 7)
    return(subdat)
  })) # }}}
  merged$channel.probes = as.factor(paste(merged$channel,merged$probes,sep='.'))
  merged$channel = as.factor(merged$channel)
  merged$probes = as.factor(merged$probes)
  merged$chip = as.factor(merged$chip)

  ch = levels(merged$probes) = c('Controls','Out-of-band')
  bychannel = as.character(sapply(c('Cy3','Cy5'), function(z) paste(z, ch)))
  levels(merged$channel.probes) = bychannel
  chcolors = c('Cy3 Controls'='darkgreen',
               'Cy3 Out-of-band'='green',
               'Cy5 Controls'='darkred',
               'Cy5 Out-of-band'='red')
  names(chcolors) = bychannel

  ( ggplot(merged, aes(x=intensity, 
                       group=channel.probes, 
                       fill=channel.probes)) +
           geom_histogram(aes(y=..density..),
                          binwidth=25, 
                          alpha=I(0.5), 
                          position=position_identity()) +
           facet_grid(chip ~ channel) + 
           scale_x_continuous(limits=xl) + 
           scale_y_continuous(breaks=NA) + 
           scale_fill_manual(values=chcolors) + 
           opts(title='Negative controls and out-of-band probe intensities',
                legend.title='Probe group') + 
           theme_bw() )

}
