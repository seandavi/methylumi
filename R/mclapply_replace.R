### Replace mcapply with reasonable default that works on windows.
.mclapply <- function(...) {
  if(.Platform$OS.type=='unix') {
    require(parallel)
    if(is.null(options("mc.cores")[[1]])) {
      message("Autodetecting CPU cores...")
      options("mc.cores"=detectCores())
    }
    return(mclapply(...))
  } else {
    return(lapply(...))
  }
}
