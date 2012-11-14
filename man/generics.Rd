\name{MethyLumi-accessors}
\alias{betas}
\alias{corplot}
\alias{betas<-}
\alias{unmethylated}
\alias{methylated}
\alias{unmethylated<-}
\alias{methylated<-}
\alias{getHistory}
\alias{pvals}
\alias{pvals<-}
\alias{QCdata}
\alias{QCdata<-}
\title{methylumi accessors}
\description{
  These functions serve as getters and setters for information in methylumi classes.
}
\usage{
betas(object)
pvals(object)
methylated(object)
unmethylated(object)
getHistory(object)
QCdata(object)
}
\arguments{
  \item{object}{an object of class MethyLumi or a subclass}
}
\details{
  See the methods definitions in \code{\linkS4class{MethyLumiSet}} and \code{\linkS4class{MethyLumiQC}} for details.  
}
\author{Sean Davis <sdavis2@mail.nih.gov>}
\seealso{\code{\link{normalizeMethyLumiSet}},
  \code{\linkS4class{MethyLumiSet}}, \code{\linkS4class{MethyLumiQC}},
  \code{\linkS4class{eSet}} 
}
\keyword{ IO }
