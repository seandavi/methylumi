\name{MethyLumiQC-class}
\docType{class}
\alias{MethyLumiQC-class}
\alias{exprs<-,MethyLumiQC,matrix-method}
\alias{exprs,MethyLumiQC-method}
\alias{initialize,MethyLumiQC-method}
\alias{QCDataOrNULL-class}
\alias{controlTypes,MethyLumiQC-method}
\alias{qcplot,MethyLumiQC-method}
\title{Class to hold Illumina Methylation QC data}
\description{This class inherits from the MethyLumi class (and
  therefore, from eSet in Biobase) and is designed to hold QC data from
  Illumina Beadstudio output.  These data can be potentially useful when
determining the cause for quality problems.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("MethyLumiQC", assayData, phenoData, featureData, experimentData, annotation, betas)}.
}
\section{Slots}{
	 \describe{
    \item{\code{assayData}:}{Object of class \code{"AssayData"}}
    \item{\code{phenoData}:}{Object of class \code{"AnnotatedDataFrame"}}
    \item{\code{featureData}:}{Object of class
      \code{"AnnotatedDataFrame"} containing the annotation columns from
    the Illumina Beadstudio output.  In particular, the names of the
    probes describe the types of control probes.}
    \item{\code{experimentData}:}{Object of class \code{"MIAME"}}
    \item{\code{annotation}:}{Object of class \code{"character"}, not
      currently used}
    \item{\code{.__classVersion__}:}{Object of class \code{"Versions"}}
  }
}
\section{Extends}{
Class \code{"\linkS4class{MethyLumi}"}, directly.
Class \code{"\linkS4class{QCDataOrNULL}"}, directly.
Class \code{"\linkS4class{eSet}"}, by class "MethyLumi", distance 2.
Class \code{"\linkS4class{VersionedBiobase}"}, by class "MethyLumi", distance 3.
Class \code{"\linkS4class{Versioned}"}, by class "MethyLumi", distance 4.
}
\section{Methods}{
  \describe{
    \item{qcplot}{\code{signature(object = "MethyLumiQC", what, ...)}:
      QC plots of various controltypes }
    \item{controlTypes}{\code{signature(object = "MethyLumiQC")}:
      determine the character vector of control types from the QCdata information }
    \item{initialize}{\code{signature(.Object = "MethyLumiQC")}}
  }
}
\author{Sean Davis <sdavis2@mail.nih.gov>}
\seealso{
  \code{\link{methylumiR}}, \code{\linkS4class{MethyLumiSet}},
  \code{\linkS4class{MethyLumiQC}}, \code{\linkS4class{eSet}}
}
\examples{
showClass("MethyLumiQC")
}
\keyword{classes}
