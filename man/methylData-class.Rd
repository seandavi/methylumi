\name{methylData-class}
\Rdversion{1.1}
\docType{class}
\alias{methylData-class}
\alias{diagnostics,methylData-method}
\alias{methylated.N,methylData-method}
\alias{methylated.N<-,methylData,matrix-method}
\alias{plotNAs,methylData-method}
\alias{pval.detect,methylData-method}
\alias{pval.detect<-,methylData,numeric-method}
\alias{unmethylated.N,methylData-method}
\alias{unmethylated.N<-,methylData,matrix-method}

\title{Class \code{"methylData"}, superclass for MethyLumiSet and MethyLumiM}
\description{
  A superclass (virtual) for MethyLumiSet and MethyLumiM.
}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Methods}{
  \describe{
    \item{diagnostics}{\code{signature(x = "methylData")}: diagnostic
      plots of data }
    \item{methylated.N}{\code{signature(object = "methylData")}:
      accessor for assayData element of the same name}
    \item{methylated.N<-}{\code{signature(object = "methylData", value =
	"matrix")}: 
      replace method for assayData element of the same name}
    \item{plotNAs}{\code{signature(object = "methylData")}: ... }
    \item{pval.detect}{\code{signature(object = "methylData")}: 
      accessor for assayData element of the same name}
    \item{pval.detect<-}{\code{signature(object = "methylData", value =
	"numeric")}: 
            replace method for assayData element of the same name}
    \item{unmethylated.N}{\code{signature(object = "methylData")}:
      accessor for assayData element of the same name}
    \item{unmethylated.N<-}{\code{signature(object = "methylData", value
	= "matrix")}: 
      replace method for assayData element of the same name}
  }
}
\author{
  Tim Triche, Jr.
}
\seealso{
  \code{\linkS4class{MethyLumiSet}},\code{\linkS4class{MethyLumiM}}
}
\examples{
showClass("methylData")
}
\keyword{classes}
