## copy from nsFilter.R from the genefilter package
rowIQRs <- function(eSet) {
  numSamp <- ncol(eSet)
  lowQ <- rowQ(eSet, floor(0.25 * numSamp))
  upQ <- rowQ(eSet, ceiling(0.75 * numSamp))
  upQ - lowQ
}

## warpping the varFilter function from the genefilter package
.varFilter <- function(eset,
                   var.func=IQR, var.cutoff=0.5,
                   filterByQuantile=TRUE, ...)
{
    if (!is.function(var.func))
        stop("'var.func' must be a function")
            
    filter.log <- new.env(parent=emptyenv())
            
    if (deparse(substitute(var.func)) == "IQR") {
        esetIqr <- rowIQRs(exprs(eset))
    } else {
        esetIqr <- apply(exprs(eset), 1, var.func)
    }
 
    if (filterByQuantile) {
        if ( 0 < var.cutoff && var.cutoff < 1 ) {
            var.cutoff = quantile(esetIqr, var.cutoff, na.rm=TRUE)
        } else stop("Cutoff Quantile has to be between 0 and 1.")
    }
    selected <- esetIqr > var.cutoff
    eset <- eset[selected, ]
    logvar <- "numLowVar"
    assign(logvar, sum(!selected), filter.log)
    list(eset=eset, filter.log=as.list(filter.log))
}


          
setMethod("varFilter", "MethyLumiSet",
          function(eset,
                   var.func=IQR, var.cutoff=0.5,
                   filterByQuantile=TRUE, ...)
          {
              .varFilter(eset, var.func, var.cutoff, filterByQuantile)
          }
)

setMethod("varFilter", "MethyLumiM",
          function(eset,
                   var.func=IQR, var.cutoff=0.5,
                   filterByQuantile=TRUE, ...)
          {
              .varFilter(eset, var.func, var.cutoff, filterByQuantile)
          }
)
