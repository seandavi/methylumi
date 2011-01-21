\name{estimateM}
\alias{estimateM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate methylation M-value matrix 
}
\description{
Estimate methylation M-value matrix from MethyLumiM-class object or eSet-class object, which include methylated and unmethylated probe intensities
}
\usage{
estimateM(methyLumiM, returnType=c("ExpressionSet", "matrix"), offset=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{methyLumiM}{ MethyLumiM-class object or eSet-class object, which include methylated and unmethylated probe intensities }
  \item{returnType}{determine whether return an ExpressionSet (MethyLumiM in this case) or matrix object}
  \item{offset}{ offset added to the  methylated and unmethylated probe intensities when estimating the M-value }
}
\details{
M-value is the log2 ratio between Illumina methylated and unmethylated probe intensities. As variations of small intensities can cause big changes in the ratio estimation, so an offset is added to methylated and unmethylated probe intensities when estimating the M-value.
 
Please check the lumi package for more details of estimateM function.
}
\value{
 A MethyLumiM or matrix object of methylation M-value
}
\references{
Du, P., Zhang, X, Huang, C.C., Jafari, N., Kibbe, W.A., Hou, L., and Lin, S.M.,  (2010) 'Comparison of Beta-value and M-value methods for quantifying methylation levels by microarray analysis', (under review)
}
\author{
Pan DU
}

\seealso{
  
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
