\name{msc.baseline.subtract}
\alias{msc.baseline.subtract}
\title{Baseline Subtraction for Mass Spectra Data}
\description{Perform baseline subtraction on batch of mass spectra data}
\usage{msc.baseline.subtract(X, \dots) }

\arguments{
  \item{X}{Spectrum data either in matrix format [nFeatures \eqn{\times}{x} 
    nSamples] or in 3D array format [nFeatures \eqn{\times}{x} nSamples 
    \eqn{\times}{x} nCopies]. Row names \code{(rownames(X)} store M/Z mass of 
    each row/feature. }
  \item{\ldots}{parameters to be passed to \code{\link[PROcess]{bslnoff}} 
    function from \pkg{PROcess} library. See details for explanation of 
    \code{breaks}, \code{qntl}, and \code{bw}. Boolean parameter \code{plot} 
    can be used to plot results.
  }
}

\details{Perform baseline subtraction for every sample in a batch of data, using 
  \code{\link[PROcess]{bslnoff}} function from \pkg{PROcess} library. The 
  \code{\link[PROcess]{bslnoff}} function splits spectrum into \code{breaks}
  number of exponentially growing regions. Baseline is calculated by appling 
  \code{\link{quantile}(\dots,probs=qntl)} to each region and smoothing the 
  results using \code{\link{loess}(\dots, span=bw, degree=1)} function. 
}

\value{ Data in the same format and size as input variable \code{X} but with 
  the subtracted baseline.
} 

\author{Jarek Tuszynski (SAIC) \email{jaroslaw.w.tuszynski@saic.com}} 
\seealso{ 
  \itemize{
  \item Part of \code{\link{msc.preprocess.run}} and 
   \code{\link{msc.project.run}} pipelines.
  \item Previous step in the pipeline was \code{\link{msc.project.read}} 
    and \code{\link{msc.rawMS.read.csv}} 
  \item Next step in the pipeline is \code{\link{msc.mass.cut}} 
  \item This function uses \code{\link[PROcess]{bslnoff}} (from \pkg{PROcess} 
  library) which is a single-spectrum baseline removal function implemented using 
  \code{\link{loess}} function.
  \item Function \code{\link[PROcess]{rmBaseline}} (from \pkg{PROcess} 
  library) can read all CSV files in directory and remove their baselines.
  }
}

\examples{
  # load "Data_IMAC.Rdata" file containing raw MS spectra 'X'  
  if (!file.exists("Data_IMAC.Rdata")) example("msc.project.read")
  load("Data_IMAC.Rdata")
  
  # run msc.baseline.subtract using 3D input
  Y = msc.baseline.subtract(X)
  cat("Size before: ", dim(X), " and after :", dim(Y), "\n")

  # test on data provided in PROcess package (2D input)
  directory  = system.file("Test", package = "PROcess")
  X = msc.rawMS.read.csv(directory) 
  Y = msc.baseline.subtract(X, plot=TRUE)
  cat("Size before: ", dim(X), " and after :", dim(Y), "\n")

}

\keyword{ts}

