\name{MethyLumiSet-class}
\Rdversion{1.1}
\docType{class}
\alias{MethyLumiSet-class}
\alias{Cy3.N,MethyLumiSet-method}
\alias{Cy5.N,MethyLumiSet-method}
\alias{QCdata,MethyLumiSet-method}
\alias{QCdata<-,MethyLumiSet}
\alias{[,MethyLumiSet-method}
\alias{betas,MethyLumiSet-method}
\alias{betas<-,MethyLumiSet,matrix-method}
\alias{boxplot,MethyLumiSet-method}
\alias{combine27k450k,MethyLumiSet,MethyLumiSet-method}
\alias{combine,MethyLumiSet,MethyLumiSet-method}
\alias{controlData,MethyLumiSet-method}
\alias{controlData<-,MethyLumiSet}
\alias{controlTypes,MethyLumiSet-method}
\alias{corplot,MethyLumiSet-method}
\alias{exprs,MethyLumiSet-method}
\alias{getHistory,MethyLumiSet-method}
\alias{hist,MethyLumiSet-method}
\alias{intensities.IB,MethyLumiSet,character-method}
\alias{intensities.IB,MethyLumiSet,missing-method}
\alias{intensities.M,MethyLumiSet,character-method}
\alias{intensities.M,MethyLumiSet,missing-method}
\alias{intensities.OOB.allelic,MethyLumiSet,character,character-method}
\alias{intensities.OOB.allelic,MethyLumiSet,missing,missing-method}
\alias{intensities.OOB,MethyLumiSet,character-method}
\alias{intensities.OOB,MethyLumiSet,missing-method}
\alias{intensities.U,MethyLumiSet,character-method}
\alias{intensities.U,MethyLumiSet,missing-method}
\alias{intensitiesByChannel,MethyLumiSet-method}
\alias{methylated,MethyLumiSet-method}
\alias{methylated<-,MethyLumiSet,matrix-method}
\alias{negctls.stderr,MethyLumiSet,character-method}
\alias{negctls.stderr,MethyLumiSet,missing-method}
\alias{negctls,MethyLumiSet,character-method}
\alias{negctls,MethyLumiSet,missing-method}
\alias{negnorm,MethyLumiSet,character-method}
\alias{negnorm,MethyLumiSet,missing-method}
\alias{normctls,MethyLumiSet-method}
\alias{pairs,MethyLumiSet-method}
\alias{plotSampleIntensities,MethyLumiSet-method}
\alias{probeNAs,MethyLumiSet-method}
\alias{qcplot,MethyLumiSet-method}
\alias{sampleNAs,MethyLumiSet-method}
\alias{show,MethyLumiSet-method}
\alias{total.intensity,MethyLumiSet-method}
\alias{unmethylated,MethyLumiSet-method}
\alias{unmethylated<-,MethyLumiSet,matrix-method}

\title{Class \code{"MethyLumiSet"} for containing Illumina methylation data}
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
    \item{\code{protocolData}:}{Object of class
      \code{"AnnotatedDataFrame"} that contains protocol information,
      including scan date if available}
  }
}
\section{Extends}{
  Class \code{"\linkS4class{MethyLumi}"}, directly.
  Class \code{"\linkS4class{methylData}"}, directly.
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
    \item{exprs}{\code{signature(object = "MethyLumiSet")}: returns m-values}
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
      determine the character vector of control types from the QCdata
      information }
    \item{Cy3.N}{\code{signature(object = "MethyLumiSet")}: ... }
    \item{Cy5.N}{\code{signature(object = "MethyLumiSet")}: ... }
    \item{combine27k450k}{\code{signature(x = "MethyLumiSet", y = "MethyLumiSet")}: ... }
    \item{controlData}{\code{signature(object = "MethyLumiSet")}: ... }
    \item{controlData<-}{\code{signature(object = "MethyLumiSet", value = "MethyLumiQC")}: ... }
    \item{featureFilter}{\code{signature(eset = "MethyLumiSet")}: ... }
    \item{intensities.IB}{\code{signature(x = "MethyLumiSet", channel = "character")}: ... }
    \item{intensities.IB}{\code{signature(x = "MethyLumiSet", channel = "missing")}: ... }
    \item{intensities.M}{\code{signature(x = "MethyLumiSet", channel = "character")}: ... }
    \item{intensities.M}{\code{signature(x = "MethyLumiSet", channel = "missing")}: ... }
    \item{intensities.OOB.allelic}{\code{signature(x = "MethyLumiSet", channel = "character", allele = "character")}: ... }
    \item{intensities.OOB.allelic}{\code{signature(x = "MethyLumiSet", channel = "missing", allele = "missing")}: ... }
    \item{intensities.OOB}{\code{signature(x = "MethyLumiSet", channel = "character")}: ... }
    \item{intensities.OOB}{\code{signature(x = "MethyLumiSet", channel = "missing")}: ... }
    \item{intensities.U}{\code{signature(x = "MethyLumiSet", channel = "character")}: ... }
    \item{intensities.U}{\code{signature(x = "MethyLumiSet", channel = "missing")}: ... }
    \item{intensitiesByChannel}{\code{signature(object = "MethyLumiSet")}: ... }
    \item{negctls.stderr}{\code{signature(object = "MethyLumiSet", channel = "character")}: ... }
    \item{negctls.stderr}{\code{signature(object = "MethyLumiSet", channel = "missing")}: ... }
    \item{negctls}{\code{signature(object = "MethyLumiSet", channel = "character")}: ... }
    \item{negctls}{\code{signature(object = "MethyLumiSet", channel = "missing")}: ... }
    \item{negnorm}{\code{signature(object = "MethyLumiSet", channel = "character")}: ... }
    \item{negnorm}{\code{signature(object = "MethyLumiSet", channel = "missing")}: ... }
    \item{normctls}{\code{signature(object = "MethyLumiSet")}: ... }
    \item{plotSampleIntensities}{\code{signature(x = "MethyLumiSet")}: ... }
    \item{probeNAs}{\code{signature(object = "MethyLumiSet")}: ... }
    \item{sampleNAs}{\code{signature(object = "MethyLumiSet")}: ... }
    \item{total.intensity}{\code{signature(object = "MethyLumiSet")}: ... }
    \item{varFilter}{\code{signature(eset = "MethyLumiSet")}: ... }
	 }

}
\author{Sean Davis & Tim Triche, Jr.}
\seealso{
  \code{\link{methylumiR}}, \code{\link{normalizeMethyLumiSet}},
  \code{\link{methylumIDAT}}, \code{\linkS4class{MethyLumiQC}},
  \code{\linkS4class{eSet}} 
}
\examples{
showClass("MethyLumiSet")
}
\keyword{classes}
