\name{methylumIDAT}
\alias{methylumIDAT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
methylumIDAT
}
\description{
Read a directory of methylumi idat files and return a \code{MethylumiSet}.
}
\usage{
methylumIDAT(barcodes = NULL, pdat = NULL, parallel = F, n = T, n.sd = F, oob = T, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{barcodes}{
A vector of barcodes to read.  Either this argument or \code{pdat} must
be specified.
}
\item{pdat}{
  A data.frame describing the samples.  A special column named
  "barcodes" can be used to specify the barcodes to be read.
}
  \item{parallel}{
If TRUE, an attempt will be made to process using multiple cores on a
multicore machine.
}
  \item{n}{
%%     ~~Describe \code{n} here~~
}
  \item{n.sd}{
%%     ~~Describe \code{n.sd} here~~
}
  \item{oob}{
%%     ~~Describe \code{oob} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
Read a set of .idat files and return a \code{MethylumiSet} object.
}
\value{
A \code{MethylumiSet} object.
}
\author{
Tim Triche, Jr.
}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (barcodes = NULL, pdat = NULL, parallel = F, n = T, 
    n.sd = F, oob = T, ...) 
{
    if (is(barcodes, "data.frame")) 
        pdat = barcodes
    if ((is.null(barcodes)) & (is.null(pdat) | (!("barcode" \%in\% 
        names(pdat))))) {
        stop("\"barcodes\" or \"pdat\" (with pdat$barcode defined) must be supplied.")
    }
    if (!is.null(pdat) && "barcode" \%in\% tolower(names(pdat))) {
        names(pdat)[which(tolower(names(pdat)) == "barcode")] = "barcode"
        barcodes = pdat$barcode
        if (any(grepl("idat", ignore.case = T, barcodes))) {
            message("Warning: filtering out raw filenames")
            barcodes = gsub("_(Red|Grn)", "", barcodes, ignore = TRUE)
            barcodes = gsub(".idat", "", barcodes, ignore = TRUE)
        }
        if (any(duplicated(barcodes))) {
            message("Warning: filtering out duplicates")
            pdat = pdat[-which(duplicated(barcodes)), ]
            barcodes = pdat$barcode
        }
    }
    else {
        if (any(grepl("idat", ignore.case = T, barcodes))) {
            message("Warning: filtering out raw filenames")
            barcodes = unique(gsub("_(Red|Grn)", "", barcodes, 
                ignore.case = TRUE))
            barcodes = unique(gsub(".idat", "", barcodes, ignore.case = TRUE))
        }
        if (any(duplicated(barcodes))) {
            message("Warning: filtering out duplicate barcodes")
            barcodes = barcodes[which(!duplicated(barcodes))]
        }
    }
    files.present = rep(TRUE, length(barcodes))
    idats = sapply(barcodes, function(b) paste(b, c("_Red", "_Grn"), 
        ".idat", sep = ""))
    for (i in colnames(idats)) for (j in idats[, i]) if (!(j \%in\% 
        list.files())) {
        message(paste("Error: file", j, "is missing for sample", 
            i))
        files.present = FALSE
    }
    stopifnot(all(files.present))
    hm27 = hm450 = 0
    hm27 = sum(grepl("_[ABCDEFGHIJKL]", barcodes))
    message(paste(hm27, "HumanMethylation27 samples found"))
    hm450 = sum(grepl("_R0[123456]C0[12]", barcodes))
    message(paste(hm450, "HumanMethylation450 samples found"))
    if (hm27 > 0 && hm450 > 0) {
        stop("Cannot process both platforms simultaneously; please run separately.")
    }
    mlumi = NChannelSetToMethyLumiSet(DFsToNChannelSet(IDATsToDFs(barcodes, 
        parallel = parallel), IDAT = TRUE, parallel = parallel), 
        parallel = parallel, n = n, oob = oob, caller = deparse(match.call()))
    if (is.null(pdat)) {
        pdat = data.frame(barcode = as.character(barcodes))
        rownames(pdat) = pdat$barcode
        pData(mlumi) = pdat
    }
    else {
        pData(mlumi) = pdat
    }
    if (!is.null(mlumi@QC)) {
        sampleNames(mlumi@QC) = sampleNames(mlumi)
    }
    return(mlumi[sort(featureNames(mlumi)), ])
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ IO }
