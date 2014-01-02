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
  Retain bead numbers?  Generally this is a waste of memory.
}
  \item{n.sd}{
  Retain probe-level standard deviations?  This is also usually a waste of RAM.
}
  \item{oob}{
  Retain out-of-band intensities for Infinium I (paired) probes?  This is useful for methylumi.bgcorr, which defaults to a normal-exponential out-of-band model.
}
\item{idatPath}{
  The path to the directory containing the idat files.
  }
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
Read a set of .idat files and return a \code{MethylumiSet} object.

NOTE: The default for computing detection p-values is to use the empirical CDF
of the negative control probes. This method makes no assumptions about the shape
of the control probe distribution, but may result in more aggressive filtering. The granularity of the resulting p-values is limited by the number of control 
probes on the platform being processed (finer on 450k, coarser on 27k).
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
\dontrun{
if(require('IlluminaHumanMethylation450k.db')) {
  barcodes <- c('6005486014_R04C02',
              '6005486023_R05C01')
  mlumi450k <- methylumIDAT(barcodes,idatPath=system.file('extdata/idat',package='methylumi')) # no normalization done
  sampleNames(mlumi450k) <- c('TCGA1','TCGA2')
  show(mlumi450k)
}
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ IO }
