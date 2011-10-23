\name{MethyLumiQC-class}
\Rdversion{1.1}
\docType{class}
\alias{MethyLumiQC-class}
\alias{Cy3.N,MethyLumiQC-method}
\alias{Cy3<-,MethyLumiQC,matrix-method}
\alias{Cy5.N,MethyLumiQC-method}
\alias{Cy5<-,MethyLumiQC,matrix-method}
\alias{QCdata<-,MethyLumiSet,MethyLumiQC-method}
\alias{combine,MethyLumiQC,MethyLumiQC-method}
\alias{controlData<-,MethyLumiSet,MethyLumiQC-method}
\alias{controlTypes,MethyLumiQC-method}
\alias{hist,MethyLumiQC-method}
\alias{intensitiesByChannel,MethyLumiQC-method}
\alias{methylated,MethyLumiQC-method}
\alias{negctls.stderr,MethyLumiQC,character-method}
\alias{negctls.stderr,MethyLumiQC,missing-method}
\alias{negctls,MethyLumiQC,character-method}
\alias{negctls,MethyLumiQC,missing-method}
\alias{negnorm,MethyLumiQC,character-method}
\alias{negnorm,MethyLumiQC,missing-method}
\alias{normctls,MethyLumiQC-method}
\alias{qcplot,MethyLumiQC-method}
\alias{unmethylated,MethyLumiQC-method}

\title{Class \code{"MethyLumiQC"} for holding Illumina methylation QC data}
\description{
  This class inherits from the MethyLumi class (and
  therefore, from eSet in Biobase) and is designed to hold QC data from
  Illumina Beadstudio output.  These data can be potentially useful when
  determining the cause for quality problems.
}
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
Class \code{"\linkS4class{eSet}"}, by class "MethyLumi", distance 2.
Class \code{"\linkS4class{VersionedBiobase}"}, by class "MethyLumi", distance 3.
Class \code{"\linkS4class{Versioned}"}, by class "MethyLumi", distance 4.
}
\section{Methods}{
  \describe{
    \item{initialize}{\code{signature(.Object = "MethyLumiQC")}}
    \item{Cy3.N}{\code{signature(object = "MethyLumiQC")}: ... }
    \item{Cy3<-}{\code{signature(object = "MethyLumiQC", value = "matrix")}: ... }
    \item{Cy5.N}{\code{signature(object = "MethyLumiQC")}: ... }
    \item{Cy5<-}{\code{signature(object = "MethyLumiQC", value = "matrix")}: ... }
    \item{QCdata<-}{\code{signature(object = "MethyLumiSet", value = "MethyLumiQC")}: ... }
    \item{combine}{\code{signature(x = "MethyLumiQC", y = "MethyLumiQC")}: ... }
    \item{controlData<-}{\code{signature(object = "MethyLumiSet", value = "MethyLumiQC")}: ... }
    \item{controlTypes}{\code{signature(object = "MethyLumiQC")}:
      determine the character vector of control types from the QCdata information }
    \item{hist}{\code{signature(x = "MethyLumiQC")}: ... }
    \item{intensitiesByChannel}{\code{signature(object = "MethyLumiQC")}: ... }
    \item{methylated}{\code{signature(object = "MethyLumiQC")}: ... }
    \item{negctls.stderr}{\code{signature(object = "MethyLumiQC", channel = "character")}: ... }
    \item{negctls.stderr}{\code{signature(object = "MethyLumiQC", channel = "missing")}: ... }
    \item{negctls}{\code{signature(object = "MethyLumiQC", channel = "character")}: ... }
    \item{negctls}{\code{signature(object = "MethyLumiQC", channel = "missing")}: ... }
    \item{negnorm}{\code{signature(object = "MethyLumiQC", channel = "character")}: ... }
    \item{negnorm}{\code{signature(object = "MethyLumiQC", channel = "missing")}: ... }
    \item{normctls}{\code{signature(object = "MethyLumiQC")}: ... }
    \item{qcplot}{\code{signature(object = "MethyLumiQC", what, ...)}:
      QC plots of various controltypes }
    \item{unmethylated}{\code{signature(object = "MethyLumiQC")}: ... }
	 }
}
\author{Sean Davis <sdavis2@mail.nih.gov>}
\seealso{
  \code{\link{methylumiR}}, \code{\linkS4class{MethyLumiSet}},
  \code{\linkS4class{MethyLumi}}, \code{\linkS4class{eSet}}
}
\examples{
showClass("MethyLumiQC")
}
\keyword{classes}
