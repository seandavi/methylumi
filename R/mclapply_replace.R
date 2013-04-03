### Replace mcapply with reasonable default that works on windows.
.mclapply <- function(...) {
  if(require(parallel)) {
    if(is.null(options("mc.cores")[[1]])) {
      message("Remember to set options('mc.cores') if you want multiple cores used where possible")
    }
    return(mclapply(...))
  } else {
    return(lapply(...))
  }
}
