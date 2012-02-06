subset450kas27k <- function(x) { # {{{
  require(IlluminaHumanMethylation450k.db)
  stopifnot(annotation(x) == 'IlluminaHumanMethylation450k')
  if(!('platform' %in% varLabels(x))) x$platform = 'HumanMethylation450'
  return(x[as.character(unlist(IlluminaHumanMethylation450k_get27k())),])
} # }}}
subset27kas450k <- function(x) { # {{{
  require(IlluminaHumanMethylation450k.db)
  stopifnot(annotation(x) == 'IlluminaHumanMethylation27k')
  if(!('platform' %in% varLabels(x))) x$platform = 'HumanMethylation27'
  return(x[as.character(unlist(IlluminaHumanMethylation450k_get27k())),])
} # }}}
subsetCommonProbes <- function(x, ...) { # {{{
  if (length(list(...)) > 0) do.call(subsetCommonProbes(list(...)))
  else if(annotation(x) == 'IlluminaHumanMethylation27k') subset27kas450k(x)
  else if(annotation(x) == 'IlluminaHumanMethylation450k') subset450kas27k(x)
  else stop(paste("Don't know how to handle chip type", annotation(x)))
} # }}}
