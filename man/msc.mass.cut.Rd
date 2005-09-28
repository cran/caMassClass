\name{msc.mass.cut}
\alias{msc.mass.cut}
\title{Remove Low Mass Portion of the Mass Spectra Data.}
\description{Remove low-mass portion of the protein mass spectra (SELDI) data. }
\usage{msc.mass.cut( X, MinMass=3000) }

\arguments{
  \item{X}{Spectrum data either in matrix format [nFeatures \eqn{\times}{x} nSamples] or in 
    3D array format [nFeatures \eqn{\times}{x} nSamples \eqn{\times}{x} nCopies]. Row names 
    \code{(rownames(X))} store M/Z mass of each row. }
  \item{MinMass}{Minimum mass threshold. All data below that mass will be 
    deleted}
 }

\details{Low-mass portion of the protein mass spectra is removed since it is 
  not expected to have any biological information, and it has large enough 
  amplitude variations that can skew normalization process. This function also 
  removes all the masses (features) where the values in all the samples are 
  identical. That happens sometimes when the ends of the samples are set to 
  zero.
}


\value{Data in the similar format as input variable \code{X} but likely with 
  fewer features.} 

\author{Jarek Tuszynski (SAIC) \email{jaroslaw.w.tuszynski@saic.com}} 

\seealso{
  \itemize{
  \item Part of \code{\link{msc.preprocess.run}} and 
   \code{\link{msc.project.run}} pipelines.
  \item Previous step in the pipeline was \code{\link{msc.baseline.subtract}} 
  \item Next step in the pipeline is \code{\link{msc.mass.adjust}} 
  }
}

\examples{
  # load "Data_IMAC.Rdata" file containing raw MS spectra 'X'  
  if (!file.exists("Data_IMAC.Rdata")) example("msc.project.read")
  load("Data_IMAC.Rdata")
  
  # run in 3D input
  Y = msc.mass.cut( X, MinMass=3000)
  cat("Size before: ", dim(X), " and after :", dim(Y), "\n")

  # test on data provided in PROcess package (2D input)
  directory  = system.file("Test", package = "PROcess")
  X = msc.rawMS.read.csv(directory) 
  Y = msc.mass.cut( X, MinMass=4000)
  cat("Size before: ", dim(X), " and after :", dim(Y), "\n")

}

\keyword{ts}
