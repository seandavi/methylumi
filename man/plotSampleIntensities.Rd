\name{plotSampleIntensities}
\alias{plotSampleIntensities}
\title{
  Plot the sample intensities.
}
\description{
  The Illumina methylation platforms all show a significant dye bias.  The 
  plotSampleIntensities method shows the density plots for the two channels
  allowing direct visualization of the effect.
}
\usage{
plotSampleIntensities(x,beta.cuts,s)
}
\arguments{
  \item{x}{an object of class MethyLumi or a subclass}
  \item{beta.cuts}{cutoffs for low and high beta values}
  \item{s}{sample number to plot}
}
\examples{
data(mldat)
plotSampleIntensities(mldat,s=1)
}
\keyword{hplot}
