\name{MethyLumi-strippers}
\alias{stripBeadNs}
\alias{stripBeadSDs}
\title{Strip excessive probe-level data from MethyLumiSets}
\description{
  450k datasets with probe-level stderrs and bead numbers can become huge.  These functions help to manage their growth in memory.
}
\usage{
stripBeadNs(midat)
stripBeadSDs(object)
}
\arguments{
  \item{object}{an object of class MethyLumi or a subclass}
}
\author{Tim Triche, Jr. <tim.triche@gmail.com>}
