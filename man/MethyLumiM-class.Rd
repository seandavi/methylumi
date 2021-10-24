\name{MethyLumiM-class}
\Rdversion{1.1}
\docType{class}
\alias{MethyLumiM-class}
\alias{MethyLumiM}
\alias{QCdata,MethyLumiM-method}
\alias{[,MethyLumiM,ANY-method}
\alias{betas,MethyLumiM-method}
\alias{combine,MethyLumiM,MethyLumiM-method}
\alias{controlData,MethyLumiM-method}
\alias{controlData<-,MethyLumiM,ANY-method}
\alias{controlTypes,MethyLumiM-method}
\alias{detection,MethyLumiM-method}
\alias{detection<-,MethyLumiM-method}
\alias{getHistory,MethyLumiM-method}
\alias{intensitiesByChannel,MethyLumiM-method}
\alias{methylated.N,MethyLumiM-method}
\alias{methylated.N<-,MethyLumiM,ANY-method}
\alias{methylated,MethyLumiM-method}
\alias{methylated<-,MethyLumiM,ANY-method}
\alias{negctls,MethyLumiM,character-method}
\alias{negctls,MethyLumiM,missing-method}
\alias{negnorm,MethyLumiM,character-method}
\alias{negnorm,MethyLumiM,missing-method}
\alias{normctls,MethyLumiM-method}
\alias{probeNAs,MethyLumiM-method}
\alias{pvals,MethyLumiM-method}
\alias{sampleNAs,MethyLumiM-method}
\alias{total.intensity,MethyLumiM-method}
\alias{unmethylated.N,MethyLumiM-method}
\alias{unmethylated.N<-,MethyLumiM,ANY-method}
\alias{unmethylated,MethyLumiM-method}
\alias{unmethylated<-,MethyLumiM,ANY-method}
\alias{dataType,MethyLumiM-method}
\alias{dataType<-,MethyLumiM,ANY-method}
\alias{dataType<-,MethyLumiM-method}
\alias{[,MethyLumiM,ANY,ANY,ANY-method}

\title{Class "MethyLumiM": for Illumina Methylation microarray data
  using logRatios}
\description{
MethyLumiM is a class inherited from \code{\link[Biobase]{ExpressionSet-class}}. It is designed for Illumina Methylation microarray data. The exprs dataMatrix included in the assayData slot of MethyLumiM object includes a matrix of M-values, which is the log2 ratio of methylated and unmethylated probe intensities. The MethyLumiM class include a boxplot function uniquely designed for two-mode histogram data. It also include a coerce function to map from \code{\link[methylumi]{MethyLumi-class}}, \code{\link[methylumi]{MethyLumiSet-class}} or other \code{\link[Biobase]{eSet-class}} inherited object to MethyLumiM class object. 
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("MethyLumiM", exprs, methylated, unmethylated, detection, methylated.N, unmethylated.N, ..., assayData)}.
The "exprs" is a matrix of M-values, which is the log2 ratio of methylated and unmethylated probe intensities; "methylated" and "unmethylated" are intensity matrix measured by methylated and unmethylalted probes of Illumina Infinium methylation microarray; "detection" is the detection p-value outputted by Illumina GenomeStudio software; "methylated.N" and "unmethylated.N" are bead numbers for methylated and unmethylalted probes. "exprs", "methylated" and "unmethylated" information are required for MethyLumiM class. When creating a new MethyLumiM object, the information of "exprs", "methylated", "unmethylated" and "detection" can also be provided directly through "assayData".
}
\section{Slots}{
  \describe{
    \item{\code{history}:}{Object of class \code{"data.frame"} recording the operation history of the LumiBatch object.}
    \item{\code{controlData}:}{Object of class \code{"MethyLumiQC"} to keep the QC probe measurement information.}
    \item{\code{dataType}:}{The type of data stored in the "exprs" data matrix in "assayData". It can be "M" (M-value), "Beta" (Beta-value") or "Intensity" (Intensity of CpG-site)}
    \item{\code{assayData}:}{Object of class \code{"AssayData"}, which includes "exprs", "methylated", "unmethylated", "detection", "methylated.N" and "unmethylated.N" data matrix}
    \item{\code{phenoData}:}{Object of class \code{"AnnotatedDataFrame"}, See \code{\link[Biobase]{eSet-class}} }
    \item{\code{featureData}:}{Object of class \code{"AnnotatedDataFrame"}, See \code{\link[Biobase]{eSet-class}} }
    \item{\code{experimentData}:}{Object of class \code{"MIAME"}, See \code{\link[Biobase]{eSet-class}} }
    \item{\code{annotation}:}{Object of class \code{"character"}, See \code{\link[Biobase]{eSet-class}} }
    \item{\code{protocolData}:}{Object of class \code{"AnnotatedDataFrame"}, See \code{\link[Biobase]{eSet-class}} }
    \item{\code{.__classVersion__}:}{Object of class \code{"Versions"}, See \code{\link[Biobase]{eSet-class}} }
  }
}
\section{Extends}{
Class \code{"\link{ExpressionSet}"}, directly.
Class \code{"\link{eSet}"}, by class "ExpressionSet", distance 2.
Class \code{"\link{VersionedBiobase}"}, by class "ExpressionSet", distance 3.
Class \code{"\link{Versioned}"}, by class "ExpressionSet", distance 4.
}
\section{Methods}{
  \describe{
    \item{boxplot}{\code{signature(x = "MethyLumiM")}: plot distribution of M-value }
    \item{coerce}{\code{signature(from = "eSet", to = "MethyLumiM")}: map from \code{\link[methylumi]{MethyLumi-class}}, \code{\link[methylumi]{MethyLumiSet-class}} or other \code{\link[Biobase]{eSet-class}} inherited object to MethyLumiM class object. MethyLumiM object will only keep "exprs", "methylated", "unmethylated" and "detection" data matrix in the assayData.}
    \item{getHistory}{\code{signature(object = "MethyLumiM")}: access the operation history of \code{MethyLumiM} object. }
    \item{initialize}{\code{signature(.Object = "MethyLumiM")}: class initialization }
    \item{methylated}{\code{signature(object = "MethyLumiM")}: retrieve the data matrix measured by methylated probes}
    \item{methylated<-}{\code{signature(object = "MethyLumiM")}: set the data matrix measured by methylated probes }
    \item{unmethylated}{\code{signature(object = "MethyLumiM")}: retrieve the data matrix measured by unmethylated probes }
    \item{unmethylated<-}{\code{signature(object = "MethyLumiM")}: set the data matrix measured by unmethylated probes }
    \item{methylated.N}{\code{signature(object = "MethyLumiM")}: retrieve the data matrix keeping the number of beads of methylated probes}
    \item{methylated.N<-}{\code{signature(object = "MethyLumiM")}: set the data matrix keeping the number of beads of methylated probes }
    \item{unmethylated.N}{\code{signature(object = "MethyLumiM")}: retrieve the data matrix keeping the number of beads of unmethylated probes }
    \item{unmethylated.N<-}{\code{signature(object = "MethyLumiM")}: set the data matrix keeping the number of beads of unmethylated probes }
    \item{detection}{\code{signature(object = "MethyLumiM")}: retrieve \code{detection} data matrix in \code{AssayData-class} }
    \item{detection<-}{\code{signature(object = "MethyLumiM")}: set \code{detection} data matrix in \code{AssayData-class}  }
    \item{controlData}{\code{signature(object = "MethyLumiM")}: retrieve the \code{controlData} in \code{MethyLumiQC-class} }
    \item{controlData<-}{\code{signature(object = "MethyLumiM")}: set \code{controlData} in \code{MethyLumiQC-class}  }
    \item{dataType}{\code{signature(object = "MethyLumiM")}: retrieve the dataType, by default it is "M", it can also be "Beta" or "Intensity" }
    \item{dataType<-}{\code{signature(object = "MethyLumiM")}: set \code{dataType} in \code{MethyLumiM-class}, the value can be "M", "Beta" or "Intensity"}
	 }
}
\references{
1. Du, P., Zhang, X, Huang, C.C., Jafari, N., Kibbe, W.A., Hou, L., and Lin, S.M.,  (2010) 'Comparison of Beta-value and M-value methods for quantifying methylation levels by microarray analysis'
}
\author{
Pan DU
}

\seealso{
 \code{\link[methylumi]{MethyLumi-class}} and \code{\link[methylumi]{MethyLumiSet-class}}
}
\examples{
showClass("MethyLumiM")
}
\keyword{classes}
