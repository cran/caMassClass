\name{msc.peaks.read.mzXML & msc.peaks.write.mzXML}
\alias{msc.peaks.read.mzXML}
\alias{msc.peaks.write.mzXML}

\title{Read / write calculated peaks heights and positions to/from mzXML Files}
\description{
 Functions to read and write mzXML files containing peaks stored in the format 
 used by peak finding functions.
}

\usage{
msc.peaks.write.mzXML(scans, filename, mzXML=NULL, ...)
msc.peaks.read.mzXML(filename, wipe=TRUE)
}

\arguments{
  \item{scans}{Peak information to be stored in mzXML file. A data-frame in the
    similar to format as returned by \code{\link{msc.peaks.find}}, containing 
    four components: 
  \itemize{
    \item Spectrum. - sample number of each peak (the same for all peaks from 
      the same sample)
    \item Peak. - peak number within each sample
    \item Intensity - peak height (intensity)
    \item Substance.Mass - x-axis position, or corresponding mass of 
    the peak measured in M/Z} 
  }
  \item{filename}{character string with name of the file.}
  \item{mzXML}{class storing partially parsed mzXML data}
  \item{wipe}{Should all scans that were returned be also deleted (wiped) from 
		mzXML record? Set to \code{TRUE} by default to minimize memory use.}
  \item{...}{additional parameters to be passed to \code{\link{write.mzXML}} 
    function (precision)}
}

\details{
  Functions \code{\link{read.mzXML}} and \code{\link{write.mzXML}} use very 
  general data type to communicate with \code{mzXML} files. Functions 
  \code{msc.rawMS.read.mzXML} and \code{msc.rawMS.write.mzXML} 
  allow passing information using data format specialized for storing 
  peak data.  
}

\value{
  Function \code{msc.peaks.read.mzXML} returns data in the same format as 
  \code{msc.peaks.write.mzXML} input parameter \code{scans}. In addition, 
  object of type \code{mzXML} is attached as "\code{mzXML}" attribute. 
  See \code{\link{read.mzXML}} for details.
  
  Functions \code{msc.peaks.write.mzXML} do not return anything. 
}

\author{Jarek Tuszynski (SAIC) \email{jaroslaw.w.tuszynski@saic.com}} 

\seealso{ 
  \itemize{
  \item \code{\link{read.mzXML}}, \code{\link{write.mzXML}} are more general 
   mzXML file reader/writer.
  \item \code{\link{msc.rawMS.read.mzXML}} & \code{\link{msc.rawMS.write.mzXML}}
   functions also read/write mzXML file, but use different data format.
  \item \code{\link{msc.peaks.read.csv}} & \code{\link{msc.peaks.write.mzXML}} 
  function can read/write peak data using CSV files. 
  }
}

\examples{
  # load "Data_IMAC.Rdata" file containing raw MS spectra 'X'  
  directory  = system.file("Test", package = "caMassClass")
  X = msc.rawMS.read.csv(directory, "IMAC_.*csv", mzXML.record=TRUE) 

  # Find Peaks and save them
  Peaks = msc.peaks.find(X) # create peak data
  msc.peaks.write.mzXML(Peaks, "peaks.mzXML", mzXML=attr(X,"mzXML"), 
                        precision='64')
  Output = msc.peaks.read.mzXML("peaks.mzXML")
  stopifnot(Output$Substance.Mass == Peaks$Substance.Mass,
            Output$Intensity      == Peaks$Intensity, 
            Output$Spectrum.      == Peaks$Spectrum., 
            Output$Peak.          == Peaks$Peak.)
  mzXML = attr(Output,"mzXML")
  strsplit(mzXML$parentFile, '\n')      # show mzXML$parentFile record
  # Suggestion: inspect 'peaks.mzXML' using any text editor          
  file.remove("peaks.mzXML")
}

\keyword{file}
