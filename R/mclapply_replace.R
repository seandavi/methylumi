### Replace mcapply with reasonable default that works on windows.
.mclapply <- function(...) {
  if(.Platform$OS.type=='unix') {
    require(parallel)
    return(mclapply(...))
  } else {
    return(lapply(...))
  }
}
