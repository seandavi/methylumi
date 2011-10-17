\name{MethyLumi-strippers}
\alias{stripOOB}
\alias{stripBeadNs}
\alias{stripBeadSDs}
\alias{stripMethyLumiSet}
\title{Strip excessive probe-level data from MethyLumiSets}
\description{
  450k datasets with probe-level stderrs, out-of-band intensities, and bead numbers can become huge.  These functions help to manage their growth in memory, at least until preprocessing and QC is completed, whereupon the summary data can be exported to a RangedData-based object of some sort for integration.
}
\usage{
stripMethyLumiSet(object)
stripBeadNs(object)
stripBeadSDs(object)
stripOOB(object)
}
\arguments{
  \item{object}{an object of class MethyLumi or a subclass}
}
\author{Tim Triche, Jr. <tim.triche@gmail.com>}
