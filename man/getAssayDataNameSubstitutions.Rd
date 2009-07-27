\name{getAssayDataNameSubstitutions}
\Rdversion{1.1}
\alias{getAssayDataNameSubstitutions}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Return a data.frame of AssayData name substitutions.
}
\description{
The Illumina methylation platforms use two distinct platforms, the "goldengate" platform and the "infinium" platform.  Each of these uses different file formats as well as different assay techologies.  To make the downstream data handling more straightforward and uniform between the two different systems, a simple mapping from the column names in the output files from the Illumina software is used to convert things from Red/Green or Cy5/Cy3 to unmethylated/methylated.  This function simply returns that mapping.
}
\usage{
getAssayDataNameSubstitutions()
}
\details{
A file in the extdata directory called "substitutions.txt" contains two columns.  The function loads this file and uses the first column as a match against column names in the data file (with the "sample part" removed).  If matched, the second column gives the replacement.  
}
\value{
A data.frame with two columns, regex and replacement. 
}
\author{
Sean Davis <seandavi@gmail.com>
}
\examples{
getAssayDataNameSubstitutions()
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ IO }