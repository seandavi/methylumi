\name{qcplot}
\alias{qcplot}
\alias{controlTypes}
\title{
  Methods for dealing with control data for Illumina methylation data.
}
\description{
  The \code{qcplot} function simply generates a plot of the control probe information for a given \code{controlType}.  
}
\usage{
qcplot(object,controltype,...)
controlTypes(object,...)
}
\arguments{
  \item{object}{An object of class \code{\linkS4class{MethyLumiSet}} or \code{\linkS4class{MethyLumiQC}}}
  \item{controltype}{A single character value representing the bead type to plot from the quality control data.  The available types are accessible via the \code{\link{controlTypes}} method.}
  \item{...}{passed to plot function}
}
\details{
  The descriptions of the various control types can be obtained from the Illumina methylation user's guides.  
}
\author{Sean Davis <sdavis2@mail.nih.gov>}
\seealso{
  \code{\linkS4class{MethyLumiSet}}, \code{\linkS4class{MethyLumiQC}}
}
\examples{
data(mldat)
controlTypes(mldat)
qcplot(mldat,controlTypes(mldat)[3])
}
\keyword{hplot}
