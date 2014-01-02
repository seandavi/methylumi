\name{methylumiR}
\alias{methylumiR}
\title{ Load data from Illumina methylation platform}
\description{
  This function is useful for loading Illumina methylation data into a
  MethyLumiSet object.  Sample information can be supplied and will
  then be incorporated into the resulting phenoData slot.
}
\usage{
methylumiR(filename, qcfile=NULL, sampleDescriptions = NULL, sep = NULL, ...)
}
\arguments{
  \item{filename}{A filename of the excel-like file from BeadStudio }
  \item{qcfile}{A filename of the excel-like file from BeadStudio }
  \item{sampleDescriptions}{A data.frame that contains at least one column, SampleID (case insensitive).  This column MUST match the part of the column headers before the \code{.Avg_Beta}, etc.  Also, if a column called SampleLabel (case insensitive), it is used for sample labels, IF the sampleLabel column contains unique identifiers}
  \item{sep}{seperator used in the BeadStudio (or GenomeStudio) output file. If it is NULL, the function will automatically estimate it.}
  \item{...}{Passed into \code{read.delim()}}
}
\details{
  Used to construct a MethyLumiSet object....
}
\value{
  A MethyLumiSet object
}
\author{Sean Davis <sdavis2@mail.nih.gov>}
\seealso{ \code{\link{MethyLumiSet-class}}, \code{\link{MethyLumiQC-class}}}
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
mldat
}
\keyword{ IO }
