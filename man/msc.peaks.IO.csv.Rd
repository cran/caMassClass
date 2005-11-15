\name{msc.peaks.read.csv & msc.peaks.write.csv}
\alias{msc.peaks.read.csv}
\alias{msc.peaks.write.csv}
\title{Read and Write Mass Spectra Peaks in CSV Format }
\description{
 Functions to read and write CSV (comma separated values) text files containing
 peaks in the format used by Ciphergen's peak file.
}
\usage{
  X = msc.peaks.read.csv(fname, mzXML.record=FALSE)
  msc.peaks.write.csv(X, fname)
}

\arguments{
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
  \item{fname}{either a character string naming a file or a connection.}
  \item{mzXML.record}{should mzXML record be created to store mata-data (input 
    file names)? }
}

\value{
  Function \code{msc.peaks.read.csv} returns peak information data frame. See
  argument \code{X} above.  If \code{mzXML.record} was set to true than mzXML 
  record with input file names will be attached to \code{X} as "\code{mzXML}" 
  attribute.
  
  Function \code{msc.peaks.write.csv} does not return anything. 
}

\author{Jarek Tuszynski (SAIC) \email{jaroslaw.w.tuszynski@saic.com}} 
\seealso{\code{\link{msc.peaks.find}} and \code{\link{msc.peaks.align}}}
\examples{
  # load "Data_IMAC.Rdata" file containing raw MS spectra 'X'  
  if (!file.exists("Data_IMAC.Rdata")) example("msc.project.read")
  load("Data_IMAC.Rdata")

  # Find Peaks and save them
  Peaks = msc.peaks.find(X) # create peak data
  msc.peaks.write.csv(Peaks, "peaks.csv")
  Pks = msc.peaks.read.csv("peaks.csv", mzXML.record=TRUE)
  stopifnot(Pks==Peaks)
  mzXML = attr(Pks,"mzXML")
  strsplit(mzXML$parentFile, '\n')      # show mzXML$parentFile record
  # Suggestion: inspect 'peaks.csv' using any text editor          
  file.remove("peaks.csv")
}
\keyword{ts}
