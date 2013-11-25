\name{MethyLumi-class}
\docType{class}
\alias{MethyLumi-class}
\alias{methylated<-,MethyLumi,matrix-method}
\alias{methylated,MethyLumi-method}
\alias{unmethylated<-,MethyLumi,matrix-method}
\alias{unmethylated,MethyLumi-method}
\alias{pvals<-,MethyLumi,matrix-method}
\alias{pvals,MethyLumi-method}
\alias{summary,MethyLumi-method}
\alias{qcplot,Methylumi-method}

\title{The base class for storing Illumina Methylation data}
\description{This class inherits from \code{\linkS4class{eSet}} from the Biobase package and
  is used as a base class for the other two methylumi classes,
  \code{\linkS4class{MethyLumiSet}} and \code{\linkS4class{MethyLumiQC}}.}
\section{Objects from the Class}{
The \code{MethyLumi} class is a virtual class and is not meant to be instantiated.  Instead, one should instantiate a \code{\linkS4class{MethyLumiSet}} or a \code{\linkS4class{MethyLumiQC}} object.
}
\section{Slots}{
	 \describe{
    \item{\code{assayData}:}{Object of class \code{"AssayData"}}
    \item{\code{phenoData}:}{Object of class \code{"AnnotatedDataFrame"}}
    \item{\code{featureData}:}{Object of class
      \code{"AnnotatedDataFrame"} that will hold the annotation columns
      from the Beadstudio output, if they are available.}
    \item{\code{experimentData}:}{Object of class \code{"MIAME"}}
    \item{\code{annotation}:}{Object of class \code{"character"}; note
      that this slot is not currently used, but may be used in the
      future to store the character name of the annotation package, if available.}
    \item{\code{.__classVersion__}:}{Object of class \code{"Versions"}}
  }
}
\section{Extends}{
Class \code{"\linkS4class{eSet}"}, directly.
Class \code{"\linkS4class{VersionedBiobase}"}, by class "eSet", distance 2.
Class \code{"\linkS4class{Versioned}"}, by class "eSet", distance 3.
}
\section{Methods}{
  \describe{
    \item{pvals<-}{\code{signature(object = "MethyLumi", value = "matrix")}: Set the assayData slot of the same name and stores the P-values from BeadStudio or as computed by methylumIDAT (using the empirical CDF of the negative controls)}
    \item{pvals}{\code{signature(object = "MethyLumi")}: Get the assayData slot of the same name}
    \item{betas<-}{\code{signature(object = "MethyLumi", value = "matrix")}: Set the assayData slot of the same name and represents the methylation values for the samples, analogous to exprs() in gene expression data.}
    \item{betas}{\code{signature(object = "MethyLumi")}: Get the assayData slot of the same name}
    \item{methylated<-}{\code{signature(object = "MethyLumi", value = "matrix")}: Set the assayData slot that represents the Methylated single-channel signal}
    \item{methylated}{\code{signature(object = "MethyLumi")}: Get the assayData slot that represents the Methylated single-channel signal}
    \item{unmethylated<-}{\code{signature(object = "MethyLumi", value = "matrix")}: Set the assayData slot that represents the Unmethylated single-channel signal}
    \item{unmethylated}{\code{signature(object = "MethyLumi")}: Get the assayData slot that represents the Unmethylated single-channel signal}
    \item{controlTypes}{\code{signature(object = "MethyLumi"}}: Find the unique control type beeds in the QCdata slot.  
    \item{qcplot}{\code{signature(object = "MethyLumi",what,...)}: Plot of QC data.  This plot can be useful for diagnosing the problems associated with specific samples or arrays.  The value for "what" is one of the control types (which can be found by using \code{controlTypes()} on the object. }
    \item{summary}{\code{signature(object = "MethyLumi",...)}}: summary method for MethyLumi objects. }
}
\author{Sean Davis <sdavis2@mail.nih.gov>}
\seealso{
	\code{\link{methylumiR}},
	\code{\linkS4class{MethyLumiSet}}, \code{\linkS4class{MethyLumiQC}},
	\code{\linkS4class{eSet}}
}
\examples{
## The class structure
showClass("MethyLumi")
## read in some data
## Read in sample information
samps <- read.table(system.file("extdata/samples.txt",
                                package = "methylumi"),sep="\t",header=TRUE)
## Perform the actual data reading
## This is an example of reading data from a 
## Sentrix Array format file (actually two files,
## one for data and one for QC probes)
mldat <- methylumiR(system.file('extdata/exampledata.samples.txt',
        package='methylumi'),
      qcfile=system.file('extdata/exampledata.controls.txt',
        package="methylumi"),
      sampleDescriptions=samps)
mldat
## Get history information
getHistory(mldat)
## Get QC data, which is another eSet-derived object
QCdata(mldat)
}
\keyword{classes}
