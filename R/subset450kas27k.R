subset.common.probes <- function(x, ...) { # {{{
  if (length(list(...)) > 0) do.call(subset.common.probes(list(...)))
  else if(annotation(x) == 'IlluminaHumanMethylation27k') subset.27k.as.450k(x)
  else if(annotation(x) == 'IlluminaHumanMethylation450k') subset.450k.as.27k(x)
  else stop(paste("Don't know how to handle chip type", annotation(x)))
} # }}}
subset.450k.as.27k <- function(x) { # {{{
  require(IlluminaHumanMethylation450k.db)
  stopifnot(annotation(x) == 'IlluminaHumanMethylation450k')
  if(!('platform' %in% varLabels(x))) x$platform = 'HumanMethylation450'
  return(x[as.character(unlist(IlluminaHumanMethylation450k_get27k())),])
} # }}}
subset.27k.as.450k <- function(x) { # {{{
  require(IlluminaHumanMethylation450k.db)
  stopifnot(annotation(x) == 'IlluminaHumanMethylation27k')
  if(!('platform' %in% varLabels(x))) x$platform = 'HumanMethylation27'
  return(x[as.character(unlist(IlluminaHumanMethylation450k_get27k())),])
} # }}}
