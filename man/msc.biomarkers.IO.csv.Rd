\name{msc.biomarkers.read.csv & msc.biomarkers.write.csv}
\alias{msc.biomarkers.read.csv}
\alias{msc.biomarkers.write.csv}
\title{Read and Write biomarker matrix in CSV format }
\description{
 Functions to read and write CSV (comma separated values) text files containing
 biomarkers (aligned peaks) in the format used by Ciphergen's biomarker file, 
 with spectra (samples) as rows, and biomarkers as columns (features).
}

\usage{
  X = msc.biomarkers.read.csv(fname)
  msc.biomarkers.write.csv(X, fname)
}

\arguments{
  \item{fname}{either a character string naming a file or a connection.}
  \item{X}{biomarker data in form of a 2D matrix (nFeatures \eqn{\times}{x} 
    nSamples) or 3D array (nFeatures \eqn{\times}{x} nSamples \eqn{\times}{x} 
    nCopies. Notice that this data is in format which is a transpose 
    of data in CSV file.}
}

\value{
  Function \code{msc.biomarkers.read.csv} returns peak information data frame. 
  See argument \code{X} above. Function \code{msc.biomarkers.write.csv} does 
  not return anything.
}
\author{Jarek Tuszynski (SAIC) \email{jaroslaw.w.tuszynski@saic.com}} 
\seealso{\code{\link{msc.biomarkers.fill}}}
\examples{
  example("msc.peaks.align", verbose=FALSE) # create biomarkers data
  X = Y$Bmrks   #  biomarkers data is stored in variable 'Y$Bmrks'
  msc.biomarkers.write.csv(X, "biomarkers.csv")
  Y = msc.biomarkers.read.csv("biomarkers.csv")
  file.remove("biomarkers.csv")
  stopifnot( all(X==Y, na.rm=TRUE) )
}
\keyword{ts}
