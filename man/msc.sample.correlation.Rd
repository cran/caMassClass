\name{msc.sample.correlation}
\alias{msc.sample.correlation}
\title{Sample Correlation}
\description{Calculates correlations between different samples and correlations 
  between different copies of the same sample}

\usage{msc.sample.correlation(X, PeaksOnly=FALSE)}

\arguments{ 
  \item{X}{Spectrum data either in matrix format [nFeatures \eqn{\times}{x} nSamples] or in 
    3D array format [nFeatures \eqn{\times}{x} nSamples \eqn{\times}{x} nCopies]. Row names 
    \code{(rownames(X))} store M/Z mass of each row. }
  \item{PeaksOnly}{Should only peaks be used in calculating the correlation? 
    In case of raw mass spectra data it does not make much sense to calculate 
    correlation of "valleys" between peaks so one can set this flag to TRUE and
    only points above sample mean will be used.} 
}

\details{
  Function calculates for each copy of each sample two variables:
  \itemize{
    \item inner correlation - average correlation between multiple copies of 
    the same sample. Inner correlation measures how similar copies are to each 
    other. For example \code{innerCor[iSamp,iCopy]} measures average 
    correlation between \code{X[,iSamp,iCopy]} and all other copies of that 
    sample. In case of one copy of the data \code{innerCor} is set to one. In 
    case of two 
    copies \code{innerCor[iSamp,1]}\code{=innerCor[iSamp,2]}
    \code{=cor(X[,iSamp,1],X[,iSamp,2])}.  
    In case of 3 copies \code{innerCor[iSamp,1]} 
    \code{= (cor(X[,iSamp,1],X[,iSamp,2])} 
    \code{+ cor(X[,iSamp,1],X[,iSamp,3]))/2} ,etc.
    \item outer correlation - average correlation between each sample and every 
    other sample within the same copy. Outer correlation measures how similar 
    each copy is to everybody else.
  }
}


\value{Returns list with two components: \code{innerCor} and \code{outerCor} 
  both of size [nSamples \eqn{\times}{x} nCopies].} 
\author{Jarek Tuszynski (SAIC) \email{jaroslaw.w.tuszynski@saic.com}} 

\seealso{
  \itemize{
  \item Part of \code{\link{msc.preprocess.run}} and 
    \code{\link{msc.project.run}} pipelines.
  \item Used by \code{\link{msc.copies.merge}} function.
  \item Uses \code{\link{cor}} function
  }
}

\examples{
   # load input data 
  if (!file.exists("Data_IMAC.Rdata")) example("msc.project.read")
  load("Data_IMAC.Rdata")
  
  # run in 3D input data using long syntax
  out = msc.mass.adjust.calc (X);
  Y   = msc.mass.adjust.apply(X, out$ShiftX, out$ScaleY, out$ShiftY)
  
  # check what happen to sample correlation
  A = msc.sample.correlation(X, PeaksOnly=TRUE)
  B = msc.sample.correlation(Y, PeaksOnly=TRUE)
  cat("Mean corelation between two copies of the same sample:\n")
  cat(" before: ", mean(A$innerCor)," after: ", mean(B$innerCor), "\n")
  cat("Mean corelation between unrelated samples:\n")
  cat(" before: ", mean(A$outerCor)," after: ", mean(B$outerCor), "\n")
 }

\keyword{ts}

