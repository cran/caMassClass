\name{msc.peaks.read.csv}
\alias{msc.peaks.read.csv}
\alias{msc.peaks.write.csv}
\title{Read and Write Mass Spectra Peaks in CSV Format }
\description{
 Functions to read and write CSV (comma separated values) text files containing
 peaks in the format used by Ciphergen's peak file.
}
\usage{
  X = msc.peaks.read.csv(fname)
  msc.peaks.write.csv(fname, X)
}

\arguments{
  \item{fname}{either a character string naming a file or a connection.}
  \item{X}{Peak information. A data-frame in the same format as returned by
    \code{\link{msc.peaks.find}}, containing five components: 
  \itemize{
    \item Spectrum.Tag - sample name of each peak
    \item Spectrum. - sample number of each peak
    \item Peak. - peak number within each sample
    \item Intensity - peak height (intensity)
    \item Substance.Mass - x-axis position, or corresponding mass of 
    the peak measured in M/Z} 
  }
}

\value{
  Function \code{msc.peaks.read.csv} returns peak information data frame. See
  argument \code{X} above. Function \code{msc.peaks.write.csv} does not return 
  anything.
}

\author{Jarek Tuszynski (SAIC) \email{jaroslaw.w.tuszynski@saic.com}} 
\seealso{\code{\link{msc.peaks.find}} and \code{\link{msc.peaks.align}}}
\examples{
  example("msc.peaks.find") # create peak data
  X = Peaks                 # Peak data is stored in variable 'Peaks'
  msc.peaks.write.csv("peaks.csv", X)
  X = msc.peaks.read.csv("peaks.csv")
  file.remove("peaks.csv")
  stopifnot(X==Peaks)
}
\keyword{ts}
