\name{msc.peaks.align}
\alias{msc.peaks.align}
\alias{msc.peaks.alignment}
\title{Align Peaks of Mass Spectra into a "Biomarker" Matrix}
\description{Align peaks from multiple protein mass spectra (SELDI) samples into 
  a single "biomarker" matrix}

\usage{
  msc.peaks.align(Peaks, SampFrac=0.3, BinSize=c(0.002, 0.008), \dots)
  msc.peaks.alignment(S,  M,  H, Tag=0, SampFrac=0.3, BinSize=c(0.002, 0.008), \dots)
}

\arguments{
  \item{Peaks}{Peak information. Could have two formats: a filename where to 
    find the data, or the data itself. In the first case, \code{Peaks} is string 
    containing path to a file saved by \code{\link{msc.peaks.find}}, 
    \code{\link[PROcess]{getPeaks}} (from \pkg{PROcess} package), or by 
    other software.  In the second case,  it is a data-frame in the same 
    format as returned by
    \code{\link{msc.peaks.find}}. A third way to pass the same input data is through
    use of \code{S, M, H} and \code{Tag} variables (described below) used by 
    \code{msc.peaks.alignment} function. }
  \item{S}{Peak sample number. Unique number of the sample the peak belongs to.
    Likely to come from \code{Peaks$Spectrum.} .}
  \item{M}{Peak center mass. Position of the peak on the x-axis. 
    Likely to come from \code{Peaks$Substance.Mass}.}
  \item{H}{Ppeak height. Likely to come from \code{Peaks$Intensity}.}
  \item{Tag}{Peak sample name. Unique name of the sample the peak belongs to.
    Likely to come from \code{Peaks$Spectrum.Tag}. 
    Optional since is used only to set column-names of output data.}
  \item{SampFrac}{After peak alignment, bins with fewer peaks than 
    \code{SampFrac*nSamp} are removed.}
  \item{BinSize}{Upper and lower bound of bin-sizes, based on expected 
    experimental variation in the mass (m/z) values. Size of any bin is 
    measured as \code{(R-L)/mean(R,L)} where L and R are masses (m/z values) of 
    left and right boundaries. All resulting bin sizes will all be between 
    \code{BinSize[1]} and \code{BinSize[2]}. Since SELDI data is often assumed 
    to have \deqn{\pm}{+-}3\% mass drift than a good bin size is twice that 
    number (0.006). 
    Same as \code{BinSize} variable in \code{\link{msc.peaks.clust}}, except for 
    default. }
  \item{\dots}{Two additional parameters that can be passed to 
    \code{\link{msc.peaks.clust}} are mostly for expert users fine-tuning the code:
    \itemize{
      \item tol - gaps bigger than \code{tol*max(gap)} are assumed to be the 
      same size as the largest gap. See details.
       \item verbose - boolean flag turns debugging printouts on. 
     }
  }
}

\details{
  Two interfaces were provided to the same function:
  \itemize{
    \item \code{msc.peaks.alignment} is a lower level function with more detailed 
    inputs and outputs. Possibly easier to customize for other purposes than 
    processing SELDI data.
    \item \code{msc.peaks.align} is a higher level function with simpler interface 
    customized for processing SELDI data.
  }
  
 This function align peaks from different samples into bins in such a way as 
 to satisfy constraints in following order:
 \itemize{
   \item bin sizes are in between \code{BinSize[1]} and \code{BinSize[2]}
    \item no two peaks from the same sample are present in the same bin
    \item bins are split in such a way as to minimize bin size and maximize 
    spaces between bins
    \item if there are multiple, equally good, ways to split a bin than bin is 
    split in such a way as to minimize number of repeats on each smaller sub-bin
  }

  The algorithm used does the following:
  \itemize{
    \item Store mass and sample number of each peak into an array
    \item Concatenate arrays from all samples and sort them according to mass 
    \item Group sets of peaks into subsets (bins). Each subset will consist of 
    peaks from different spectra that have similar mass. That is done by puting all
    peaks into a single bin and recursivly going through the following steps:
    \itemize{
      \item Check size of the current bin: if it is too small than we are done, 
      if it is too big than it will be split and if it is already in the 
      desired range than it will be split only if multiple peaks from the same 
      sample are present.
      \item If bin needs to be split than find the biggest gap between peaks
      \item If multiple gaps were found with the same size as the largest gap 
      (or within \code{tol} tolerance from it) than minimizes number of 
      multiple peaks from the same sample after cut 
      \item Divide the bin into two sub-bins: to the left and to the right of 
      the biggest gap
      \item Recursively repeat the above four steps for both sub-bins
    }
    \item Store peaks into 2D array (bins by samples) 
    \item Remove bins with fewer peaks than \code{SampFrac*nSamp}
  }
  The algorithm for peak alignment is described as recursive algorithm but the 
  actual implementation uses internal stack, instead in order to increase speed.
}


\value{ 
  \item{Bmrks}{Biomarker matrix containing one sample per column and one 
    biomarker per row. If a given sample does not have a peak in some bin than 
    \code{NA} is inserted.}
  \item{BinBounds}{Mass of left-most and right-most peak in the bin}
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
    \item Input comes most Likely from: \code{\link{msc.peaks.find}}, 
    \code{\link[PROcess]{getPeaks}} (from \pkg{PROcess} package), or Ciphergen's 
    software
    \item Output can be processed further by: \code{\link{msc.biomarkers.fill}} or 
    \code{\link{msc.copies.merge}} 
    \item Part of \code{\link{msc.preprocess.run}} pipeline
    \item Uses \code{\link{msc.peaks.clust}} function to do most of the work
    \item Uses \code{\link{msc.peaks.read.csv}} function to read peak file
    \item Uses \code{\link{msc.biomarkers.write.csv}} function to save results
    \item Function \code{\link[PROcess]{pk2bmkr}} from \pkg{PROcess} package 
     performs similar function.
   }
}

\examples{
  # load "Data_IMAC.Rdata" file containing raw MS spectra 'X'  
  if (!file.exists("Data_IMAC.Rdata")) example("msc.project.read")
  load("Data_IMAC.Rdata")

  # Find and Align peaks
  Peaks = msc.peaks.find(X)
  cat(nrow(Peaks), "peaks were found in", Peaks[nrow(Peaks),2], "files.\n")
  Y = msc.peaks.align(Peaks)
  print( t(Y$Bmrks) , na.print=".",  digits=2)
  stopifnot( dim(Y$Bmrks)==c(22, 40) )
}

\keyword{ts}

