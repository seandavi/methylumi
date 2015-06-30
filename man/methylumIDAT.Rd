\name{methylumIDAT}
\alias{methylumIDAT}
\title{
methylumIDAT
}
\description{
Read a directory of methylumi idat files and return a \code{MethylumiSet}.
}
\usage{
methylumIDAT(barcodes = NULL, pdat = NULL, parallel = F, n = F, n.sd =
F, oob = T, idatPath=getwd(), ...)
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
  Keep the bead numbers?  (Default: no) 
}
  \item{n.sd}{
  Keep the bead-level SD?  (Default: no)  
}
  \item{oob}{
  Keep the out-of-band (OOB) or opposite-channel signals?  (Default: yes)
}
\item{idatPath}{
  The path to the directory containing the idat files.
  }
\item{\dots}{
  Additional arguments to be passed to sub-functions.
}
}
\details{
Read a set of .idat files and return a \code{MethylumiSet} object.  If you use this function to any significant degree in your analysis, we would appreciate your citing the paper describing it, "Low-level processing of Illumina Infinium DNA methylation beadarrays", TJ Triche, DJ Weisenberger, D Van Den Berg, KD Siegmund, and PW Laird, Nucleic acids research, 2013.
}
\value{
A \code{MethylumiSet} object.
}
\author{
Tim Triche, Jr.
}
\seealso{
  The ``methylumi450k`` vignette: 
  vignette("methylumi450k", package="methylumi")
}
\examples{
\dontrun{
if(require('IlluminaHumanMethylation450k.db')) {
  barcodes <- c('6005486014_R04C02',
              '6005486023_R05C01')
  lumi450k <- methylumIDAT(barcodes,idatPath=system.file('extdata/idat',package='methylumi')) # no normalization done
  sampleNames(lumi450k) <- c('TCGA1','TCGA2')
  show(lumi450k)
}
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ IO }
