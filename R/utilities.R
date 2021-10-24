setClassUnion("methylData", c("MethyLumiSet","MethyLumiM"))
setClassUnion("NameOrDefault", c("character","missing"))
setClassUnion("NumericDefault", c("numeric","missing"))
setClassUnion("LogicalDefault", c("logical","missing"))

# after qsummary() in package "qvalue", but for B-H adjusted p-value summaries
psummary <- function (pvals, cuts=c(0.001,0.01,0.025,0.05,0.1,1), 
                      digits=getOption("digits"), ...) { # {{{
    cat("Cumulative number of significant p-values:\n")
    cat("\n")
    bh.p <- p.adjust(pvals, method="BH")
    counts <- sapply(cuts, function(p) c(`p-value` = sum(bh.p < p)))
    names(counts) <- paste("<", cuts, sep = "")
    print(counts)
    cat("\n")
    invisible(pvals)
} # }}}

# get the URL of a Bioconductor package:
getBiocPkgUrl <- function(pkg) { # {{{
  paste0("http://www.bioconductor.org/packages/release/bioc/html/",
         pkg,".html")
} # }}}

# useful for feature selection if not using shrinkage
sdmax <- function(y) {
  sd(na.omit(y)) / sqrt(mean(na.omit(y)) * (1 - mean(na.omit(y))))
}

# utility functions I tend to use unconsciously
shift <- function(x) if (length(x) > 1) x[[1]] else x
unshift <- function(x,y) append(y,x)
push <- function(x,y) append(x,y)
pop <- function(x) if (length(x) > 1) x[[length(x)]] else x

# utility operators
"%i%" <- intersect
"%d%" <- setdiff
"%u%" <- union

# using this is a bad habit of mine 
"%nin%" <- "%notin%" <- function(x,table) match(x, table, nomatch=0) == 0

# the following are rather handy for constructing MethyLumiSets and such
t.submit <- function() as.character(Sys.time())
t.finish <- function() as.character(format(Sys.time(), "%H:%M:%S"))

## handy for grabbing all IDAT files in an Illumina directory
getBarcodes <- function(path=".") { # {{{
  oldwd <- getwd()
  setwd(path)
  barcds <- unique(gsub("_(Red|Grn).idat","",
                        list.files(path=path,patt="idat")))
  setwd(oldwd)
  names(barcds) <- barcds
  return(barcds)
} # }}}

# handy for both m-values and beta-values (you heard it here first!)
beta.mme <- function(x, w=NULL, ...) { # {{{
  if(any(is.na(x))) stop("Cannot handle NA values!")
  ## a lie: we COULD, but let the user call impute.knn()
  if(is.null(w)) w <- rep(1 / length(x), length(x))
  xb <- weighted.mean(x, w, na.rm=T)
  s2 <- sum(w * ( (x - xb) ** 2)) / sum(w)
  a <- xb * ( ( (xb * (1 - xb)) / s2) - 1)
  b <- (1 - xb) * ( ( (xb * (1 - xb)) / s2) - 1)
  mme <- c(a=max(a,0), b=max(b,0))
  return(mme)
} # }}}

beta.mode <- function(a, b) { # {{{
  if(length(a) > 1) {
    b <- a[2]
    a <- a[1]
  }
  if(all(c(a,b) < 1)) return(NaN)
  if(a < 1) return(1)
  if(b < 1) return(0)
  return( (a - 1) / (a + b - 2) )
} # }}}

## from Smithson & Verkuilen 2006; shrink to 0.5, the mean, or the mode
##
beta.transform <- function(x, w=NULL, to.mean=TRUE, to.mode=FALSE, s=0.5){ #{{{

  stopifnot( (max(x) <= 1) && (min(x) >= 0) )
  n <- length(x)
  if( is.null(w) ) w <- rep(1, n)
  else w <- (w / (sum(w) / n))
  if(to.mean) s <- weighted.mean(x, w, na.rm=T)
  if(to.mode) s <- pmax(0.01, pmin(0.99, beta.mode(beta.mme(x, w))))
  return( ( ( x * ( n - 1 ) ) + s) / n )

} # }}}

## I'm really, really lazy
#' Total convenience function for processing IDATs like tcga
#' 
#' @param IDATs character() of idat files
#' 
tcgaPipeline <- function(IDATs) { 
  as(normalizeMethyLumiSet(
       stripMethyLumiSet(
         methylumi.bgcorr(
           methylumIDAT(IDATs)))),
     "RangedSummarizedExperiment")
}
