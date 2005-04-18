\name{msc.peaks.clust}
\alias{msc.peaks.clust}
\title{Clusters Peaks of Mass Spectra}
\description{Clusters peaks from multiple protein mass spectra (SELDI) samples}
\usage{
  msc.peaks.clust(dM, S, BinSize=c(0,sum(dM)), tol=0.97, verbose=FALSE) 
}

\arguments{
  \item{S}{Peak sample number, used to identify the spectrum the peak come from.}
  \item{dM}{Distance between sorted peak positions (masses, m/z).}
  \item{BinSize}{Upper and lower bound of bin-sizes, based on expected 
    experimental variation in the mass (m/z) values. Size of any bin is 
    measured as \code{(R-L)/mean(R,L)} where \code{L} and \code{R} are masses 
    (m/z values) of left and right boundaries.
    All resulting bin sizes will be between \code{BinSize[1]} and 
    \code{BinSize[2]}. Default is \code{c(0,sum(dM))} which ensures that no 
    \code{BinSizes} is not being used.}
  \item{tol}{gaps bigger than \code{tol*max(gap)} are assumed to be the same 
    size as the largest gap. See details.}
  \item{verbose}{boolean flag turns debugging printouts on.} 
}

\details{
  This is a low level function used by \code{msc.peaks.alignment} and not intended to
  be directly used by many users. However it might be usefull for other code 
  developers. It clusters peaks from different samples into bins in 
  such a way as to satisfy constraints in following order:
 \itemize{
   \item bin sizes are in between \code{BinSize[1]} and \code{BinSize[2]}
    \item no two peaks from the same sample are present in the same bin
    \item bins are split in such a way as to minimize bin size and maximize 
    spaces between bins
    \item if there are multiple, equally good, ways to split a bin than bin is 
    split in such a way as to minimize number of repeats on each smaller sub-bin
  }
}

\value{ 
  The output is binary array of the same size as \code{dM} and \code{S} where 
  left boundaries of each clusters-bin (biomarker) are marked
} 

\references{
  The initial version of this function started as implementation of algorithm
  described on webpage of Virginia Prostate Center (at Virginia Medical School)
  documenting their PeakMiner Software. See 
  \url{http://www.evms.edu/vpc/seldi/peakminer.pdf}
} 

\author{Jarek Tuszynski (SAIC) \email{jaroslaw.w.tuszynski@saic.com}} 

\seealso{
  \itemize{
  \item Part of \code{\link{msc.preprocess.run}} and 
    \code{\link{msc.project.run}} pipelines.
  \item Previous step in the pipeline was \code{\link{msc.peaks.find}} 
  \item Next step in the pipeline is \code{\link{msc.peaks.align}} and 
    \code{\link{msc.biomarkers.fill}} 
  \item Part of \code{\link{msc.peaks.align}} function
  }
}

\examples{
  # example with simple made up data (18 peaks, 3 samples)
  M = c(1,5,8,12,17,22, 3,5,7,11,14,25, 1, 5, 7,10,17,21) # peak position/mass
  S = c(1,1,1, 1, 1, 1, 2,2,2, 2, 2, 2, 3, 3, 3, 3, 3, 3) # peak's sample number
  idx = sort(M, index=TRUE)$ix;  # sort peaks by mass
  M   = M[idx];               # sorted mass
  S   = S[idx];               # arrange sample numbers in the same order
  bin = msc.peaks.clust(diff(M), S, verbose=TRUE) 
  rbind(S,M,bin)              # show results
  
  # use the results to align peaks into biomarkers matrix
  Bmrks = matrix(NA,sum(bin),max(S)); # init feature (biomarker) matrix
  bin   = cumsum(bin);                # find bin numbers for each peak in S array
  for (j in 1:length(S))              # Bmrks usually store height H of each peak
    Bmrks[bin[j], S[j]] =  1;         # but in this example it will be "1"
  Bmrks
}

\keyword{ts}
