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
msc.peaks.read.mzXML(filename)
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
  \item{...}{additional parameters to be passed to \code{\link{write.mzXML}} 
    function (precision)}
}

\details{
  Functions \code{\link{read.mzXML}} and \code{\link{write.mzXML}} use very 
  general data type to comunicate with \code{mzXML} files. Functions 
  \code{msc.rawMS.read.mzXML} and \code{msc.rawMS.write.mzXML} 
  allow passing informations using data format specialized for storing 
  peak data.  
}

\value{
  Function \code{msc.peaks.read.mzXML} returns list with folowing objects: 
  \item{scans}{Data in the same format as \code{msc.peaks.write.mzXML} input 
    parameter \code{scans}.}
  \item{mzXML}{Object of type \code{mzXML}. See \code{\link{read.mzXML}} for 
    details }
  
  Functions \code{msc.peaks.write.mzXML} do not return anything. 
}

\author{Jarek Tuszynski (SAIC) \email{jaroslaw.w.tuszynski@saic.com}} 

\seealso{ 
  \itemize{
  \item \code{\link{read.mzXML}}, \code{\link{write.mzXML}} are more genaral 
   mzXML file reader/writer.
  \item \code{\link{msc.rawMS.read.mzXML}} & \code{\link{msc.rawMS.write.mzXML}}
   functions also read/write mzXML file, but use different data format.
  \item \code{\link{msc.peaks.read.csv}} & \code{\link{msc.peaks.write.mzXML}} 
  function can read/write peak data using CSV files. 
  }
}

\examples{
  # load "Data_IMAC.Rdata" file containing raw MS spectra 'X'  
  if (!file.exists("Data_IMAC.Rdata")) example("msc.project.read")
  load("Data_IMAC.Rdata")
  
  # Find Peaks and save them
  Peaks = msc.peaks.find(X) # create peak data
  msc.peaks.write.mzXML(Peaks, "peaks.mzXML", precision='64')
  X = msc.peaks.read.mzXML("peaks.mzXML")
  stopifnot(X$scans$Substance.Mass == Peaks$Substance.Mass,
            X$scans$Intensity == Peaks$Intensity, 
            X$scans$Spectrum.==Peaks$Spectrum., X$scans$Peak.==Peaks$Peak.)
  # Suggestion: inspect 'peaks.mzXML' using any text editor          
  file.remove("peaks.mzXML")
}

\keyword{file}
