### A parallel lapply that also works on Windows.
###
### This used to wrap parallel::mclapply behind require(parallel), falling back
### to lapply, and nagged about options('mc.cores'). mclapply() forks, so on
### Windows it silently degrades to serial anyway. BiocParallel::bplapply picks
### an appropriate backend per platform and is configured the standard
### Bioconductor way, via register() or the BPPARAM argument.
###
### Exported, so the signature is kept: callers pass the same arguments they
### would pass to lapply().
.mclapply <- function(...) {
  BiocParallel::bplapply(...)
}
