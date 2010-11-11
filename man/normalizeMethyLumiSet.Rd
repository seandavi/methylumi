\name{normalizeMethyLumiSet}
\alias{normalizeMethyLumiSet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Normalize a MethyLumiSet, accounting for dye bias}
\description{
  The Illumina GoldenGate methylation platform uses two colors, one to represent
  the unmethylated state and the other to represent the methylated
  state.  This function corrects that dye bias and recalculates the
  betas based on the corrected intensities.  It is probably not optimal for
  Infinium data since the methylated and unmethylated signals are in the
  same channel for that platform.  
}
\usage{
normalizeMethyLumiSet(x, beta.cuts = c(0.2, 0.8), mapfun = c("atan", "ratio"))
}
\arguments{
  \item{x}{A MethyLumiSet object }
  \item{beta.cuts}{Two numeric values with the first less than the
    second and between 0 and 1, representing the beta cutoffs that will
    be used when determining the median intensities to which to
    correct.  See details below.}
  \item{mapfun}{Either "atan" or "ratio".  See details below.}
}
\details{
  The Illumina GoldenGate methylation platform uses two colors, one to represent
  the unmethylated state and the other to represent the methylated
  state.  This function corrects that dye bias and recalculates the
  betas based on the corrected intensities.

  As a first step, the medians for each of Cy3 and Cy5 are calculated at
  high and low betas, representing the (nearly) fully methylated state
  and the (nearly) fully unmethylated states.  Values of Cy3 and Cy5
  that are negative are set to zero for this process.  Then, the Cy5
  medians are adjusted to match those of the Cy3 channel, thereby
  correcting the dye bias.

  To map the new intensities back to betas, one of two map functions can
  be used.  The default is the atan(Cy3/Cy5).  The ratio maps using the
  function (Cy3/Cy3+Cy5).  The differences should be very small, but we
  feel that the atan map function is probably the mathematically
  appropriate way of doing this.

  Note that this normalization method is not optimal for Infinium methylation
  data and should probably not be used for that platform.
}
\value{
  A new "MethyLumiSet" that contains the corrected betas and the
  adjusted intensities.
}
\author{Sean Davis <sdavis2@mail.nih.gov>}
\examples{
## Read in sample information
samps <- read.table(system.file("extdata/samples.txt",
                                package = "methylumi"),sep="\t",header=TRUE)
## Perform the actual data reading
## This is an example of reading data from an 
## Sentrix Array format file (actually two files,
## one for data and one for QC probes)
mldat <- methylumiR(system.file('extdata/exampledata.samples.txt',package='methylumi'),
                    qcfile=system.file('extdata/exampledata.controls.txt',package="methylumi"),
                    sampleDescriptions=samps)
mldatnorm <- normalizeMethyLumiSet(mldat)
}
\keyword{ manip }
