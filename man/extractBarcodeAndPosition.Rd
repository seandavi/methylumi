\name{extractBarcodeAndPosition}
\alias{extractBarcodeAndPosition}
\title{Extract the Barcode and Position Information from Sentrix ID}
\description{
  The sentrix IDs from an illumina sentrix array contain positional
  information that might be useful.  This function simply extracts that
  information from the ID itself.
}
\usage{
extractBarcodeAndPosition(sentrixids)
}
\arguments{
  \item{sentrixids}{A character vector of sentrix IDs that look
    like: 1632405013\_R001\_C001 }
}
\value{
  A data.frame with three columns:
  \item{sentrix}{numeric, the sentrix ID}
  \item{row}{numeric, the sentrix row}
  \item{column}{numeric, the sentrix column}
}
\author{Sean Davis <sdavis2@mail.nih.gov>}
\seealso{\code{methylumiR}}
\examples{
extractBarcodeAndPosition(c('12341234_R001_C001'))
}
\keyword{ manip }
