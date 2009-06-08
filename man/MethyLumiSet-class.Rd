\name{MethyLumiSet-class}
\docType{class}
\alias{MethyLumiSet-class}
\alias{[,MethyLumiSet-method}
\alias{betas<-,MethyLumiSet,matrix-method}
\alias{betas,MethyLumiSet-method}
\alias{boxplot,MethyLumiSet-method}
\alias{combine,MethyLumiSet,MethyLumiSet-method}
\alias{corplot,MethyLumiSet-method}
\alias{exprs<-,MethyLumiSet,matrix-method}
\alias{exprs,MethyLumiSet-method}
\alias{getHistory,MethyLumiSet-method}
\alias{hist,MethyLumiSet-method}
\alias{initialize,MethyLumiSet-method}
\alias{pairs,MethyLumiSet-method}
\alias{plotSampleIntensities,MethyLumiSet-method}
\alias{QCdata<-,MethyLumiSet,MethyLumiQC-method}
\alias{QCdata,MethyLumiSet-method}
\alias{show,MethyLumiSet-method}
\alias{methylated<-,MethyLumiSet,matrix-method}
\alias{methylated,MethyLumiSet-method}
\alias{unmethylated<-,MethyLumiSet,matrix-method}
\alias{unmethylated,MethyLumiSet-method}
\alias{qcplot,MethyLumiSet-method}
\alias{controlTypes,MethyLumiSet-method}

\title{Class to hold Illumina Methylation data}
\description{This class inherits from the MethyLumi class (and
  therefore, from eSet in Biobase) and is designed to hold both the
  intensities and the calculated betas, as well as pvalues if present.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("MethyLumiSet", assayData, phenoData, featureData, experimentData, annotation, betas)}.  An object of this type is the main storage class for methylation data from Illumina.  Subsetting, etc., works as normal (rows represent genes, columns represent samples).  There is also a rudimentary history tracking system that is modeled after that from the lumi package.
}
\section{Slots}{
	 \describe{
    \item{\code{QC}:}{Object of class \code{"QCDataOrNULL"}, containing
      either the \code{\linkS4class{MethyLumiQC}} object or NULL}
    \item{\code{history}:}{Object of class \code{"data.frame"},
      containing a running history of transforms to the data contained herein}
    \item{\code{assayData}:}{Object of class \code{\linkS4class{AssayData}}}
    \item{\code{phenoData}:}{Object of class \code{\linkS4class{AnnotatedDataFrame} }}
    \item{\code{featureData}:}{Object of class
      \code{\linkS4class{AnnotatedDataFrame}}, containing the annotation columns from the Illumina Beadstudio output}
    \item{\code{experimentData}:}{Object of class \code{\linkS4class{MIAME}}}
    \item{\code{annotation}:}{Object of class \code{"character"}, not
      currently used}
    \item{\code{.__classVersion__}:}{Object of class \code{"Versions"}}
  }
}
\section{Extends}{
Class \code{"\linkS4class{MethyLumi}"}, directly.
Class \code{"\linkS4class{eSet}"}, by class "MethyLumi", distance 2.
Class \code{"\linkS4class{VersionedBiobase}"}, by class "MethyLumi", distance 3.
Class \code{"\linkS4class{Versioned}"}, by class "MethyLumi", distance 4.
}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "MethyLumiSet")}: subsetting, genes as rows, samples as columns}
    \item{betas<-}{\code{signature(object = "MethyLumiSet", value = "matrix")}: Set the assayData slot of the same name}
    \item{betas}{\code{signature(object = "MethyLumiSet")}: Get the assayData slot of the same name}
    \item{boxplot}{\code{signature(x = "MethyLumiSet")}: boxplot of all
      sample betas}
    \item{combine}{\code{signature(x = "MethyLumiSet", y = "MethyLumiSet")}}
    \item{corplot}{\code{signature(x = "MethyLumiSet")}}
    \item{exprs<-}{\code{signature(object = "MethyLumiSet", value =
	"matrix")}: alias for "<-betas"}
    \item{exprs}{\code{signature(object = "MethyLumiSet")}: alias for "betas"}
    \item{getHistory}{\code{signature(object = "MethyLumiSet")}: returns
    a data.frame containing the history for this object}
    \item{hist}{\code{signature(x = "MethyLumiSet")}: histogram of the
      betas for the data}
    \item{initialize}{\code{signature(.Object = "MethyLumiSet")}}
    \item{pairs}{\code{signature(x = "MethyLumiSet")}: pairs plot of the
    betas for the object.  Note that pairs plots of more than a few samples are not very useful.}
    \item{plotSampleIntensities}{\code{signature(x = "MethyLumiSet")}:
      The intensities as output by the Beadstudio software often show a
      considerable amount of dye bias.  This method shows a graphical
      example  of this dye bias.  In short, for each of the Cy3 and Cy5
      channels, a cutoff in beta is used to calculate which Cy3 and Cy5
      values should be plotted at high-methylation and low-methylation
      status.  Any offset between Cy3 and Cy5 when plotted in this way
      likely represents dye bias and will lead to biases in the estimate
      of beta. }
    \item{QCdata<-}{\code{signature(object = "MethyLumiSet", value =
	"MethyLumiQC")}: assign QC data to the QC slot }
    \item{QCdata}{\code{signature(object = "MethyLumiSet")}: retrieve
      the QC data. }
    \item{show}{\code{signature(object = "MethyLumiSet")}}
    \item{methylated<-}{\code{signature(object = "MethyLumiSet", value =
	"matrix")}: Set the assayData slot associated with methylated intensity}
    \item{methylated}{\code{signature(object = "MethyLumiSet")}: Get the assayData slot associated with methylated intensity}
    \item{unmethylated<-}{\code{signature(object = "MethyLumiSet", value = "matrix")}: Set the assayData slot associated with unmethylated intensity}
    \item{unmethylated}{\code{signature(object = "MethyLumiSet")}: Get the assayData slot associated with unmethylated intensity}
    \item{qcplot}{\code{signature(object = "MethyLumiSet", what, ...)}:
      QC plots of various controltypes }
    \item{controlTypes}{\code{signature(object = "MethyLumiSet")}:
      determine the character vector of control types from the QCdata information }
}
}
\author{Sean Davis <sdavis2@mail.nih.gov>}
\seealso{
  \code{\link{methylumiR}}, \code{\link{normalizeMethyLumiSet}},
  \code{\linkS4class{MethyLumiSet}}, \code{\linkS4class{MethyLumiQC}},
  \code{\linkS4class{eSet}} 
}
\examples{
showClass("MethyLumiSet")
}
\keyword{classes}
