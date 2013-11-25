\name{methylumiMU}
\alias{methylumiMU}
\title{
methylumiMU
}
\description{
Read a (perhaps gzipped) file of M and U intensities into a \code{MethylumiSet}.
}
\usage{
methylumiMU(filename, annotation='IlluminaHumanMethylation450k', ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{
    The (possibly compressed) file which should be turned into a MethyLumiSet.
  }
  \item{annotation}{
    The platform annotation (defaults to 'IlluminaHumanMethylation450k')
  }
  \item{\dots}{
    Any other arguments to be passed on.
  }
}
\details{
Read a (perhaps gzipped) file of M and U intensities into a \code{MethylumiSet}.
}
\value{
A \code{MethylumiSet} object.
}
\author{
Tim Triche, Jr.
}
\examples{
  \dontrun{
    filename <- "NBMx4.vs.MDSx4.450k.txt.gz"
    lumi450k <- methylumiMU(filename) 
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ IO }
