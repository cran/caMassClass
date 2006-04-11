\name{msc.rawMS.read.mzXML & msc.rawMS.write.mzXML}
\alias{msc.rawMS.read.mzXML}
\alias{msc.rawMS.write.mzXML}

\title{Read / write raw Mass Spectra to/from mzXML Files}
\description{Read / write raw protein mass spectra to/from mzXML Files}

\usage{
msc.rawMS.write.mzXML(scans, filename, mzXML=NULL, ...)
msc.rawMS.read.mzXML(input, scanIdx=NULL, wipe=TRUE)
}

\arguments{
  \item{scans}{data to be stored in mzXML file, in 
  form of a 2D matrix (nFeatures \eqn{\times}{x} nSamples) or 3D array
  (nFeatures \eqn{\times}{x} nSamples \eqn{\times}{x} nCopies).
  }
  \item{filename}{character string with name of the file (connection)}
  \item{mzXML}{class storing partially parsed mzXML data}
  \item{input}{Either \code{mzXML} object, or character string with name of the 
    file (connection)}
  \item{scanIdx}{List of scans to return. Optional. By default all will be 
    returned, but one can choose only a subset using this argument.}
  \item{wipe}{Should all scans that were returned be also deleted (wiped) from 
		mzXML record? Set to \code{TRUE} by default to minimize memory use.}
  \item{...}{additional parameters to be passed to \code{\link{write.mzXML}} 
    function}
}

\value{
  Function \code{msc.rawMS.read.mzXML} returns data in the matrix format (nFeatures 
  \eqn{\times}{x} nSamples) \code{\link{rownames}(scan)} storing masses (M/Z)
  of each feature.  In addition, object of type \code{mzXML} is attached 
  as "\code{mzXML}" attribute. See \code{\link{read.mzXML}} for details.
  
  Functions \code{msc.rawMS.write.mzXML} do not return anything. 
}

\author{Jarek Tuszynski (SAIC) \email{jaroslaw.w.tuszynski@saic.com}} 

\seealso{ 
  \itemize{
  \item \code{\link{read.mzXML}}, \code{\link{write.mzXML}} are more general 
   mzXML file reader/writer.
  \item \code{\link{msc.peaks.read.mzXML}}, \code{\link{msc.peaks.write.mzXML}}
   functions also read/write mzXML file, but use different data format.
  \item \code{\link{msc.rawMS.read.csv}} function can read raw MS files from 
  CSV files. 
  }
}

\examples{
  # load "Data_IMAC.Rdata" file containing raw MS spectra 'X'  
  if (!file.exists("Data_IMAC.Rdata")) example("msc.project.read")
  load("Data_IMAC.Rdata")
  
  # save raw MS data as mzXML using 32-bit precision
  msc.rawMS.write.mzXML(X, "rawMS32.mzXML",  precision="32")
  Y = msc.rawMS.read.mzXML("rawMS32.mzXML")
  dim(Y) = c(nrow(Y), ncol(Y)/2, 2)
  stopifnot(all.equal(X,Y, tolerance=1e-5, check.attributes = FALSE)) # about the same
  
  # save raw MS data as mzXML using 64-bit precision
  msc.rawMS.write.mzXML(X, "rawMS64.mzXML",  precision="64")
  Y = msc.rawMS.read.mzXML("rawMS64.mzXML")
  dim(Y) = c(nrow(Y), ncol(Y)/2, 2)
  stopifnot(X==Y)                           # exactly the same
  
  # Suggestion: inspect 'rawMS32.mzXML' and 'rawMS64.mzXML' using a text editor
  file.remove("rawMS32.mzXML") # delete temporary files 
  file.remove("rawMS64.mzXML") 
}

\keyword{file}
