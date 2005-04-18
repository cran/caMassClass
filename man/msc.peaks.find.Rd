\name{msc.peaks.find}
\alias{msc.peaks.find}
\title{Find Peaks of Mass Spectra}
\description{Find Peaks in a Batch of Protein Mass Spectra (SELDI) Data.}
\usage{msc.peaks.find(X, PeakFile=0, SNR=2, span=c(81,11), zerothresh=0.9) }

\arguments{
  \item{X}{Spectrum data either in matrix format [nFeatures \eqn{\times}{x} nSamples] or in 
    3D array format [nFeatures \eqn{\times}{x} nSamples \eqn{\times}{x} nCopies]. Row names 
    \code{(rownames(X))} store M/Z mass of each row. }
  \item{PeakFile}{optional filename. If provided than CSV file will be created 
    in the same format as Ciphergen's peak-info file, with following columns 
    of data: "Spectrum.Tag", "Spectrum.", "Peak.", "Intensity" and 
    "Substance.Mass".   }
  \item{SNR}{signal to noise ratio (z-score) criterion for peak detection. 
    Similar to \code{SoN} variable in \code{\link[PROcess]{isPeak}} from \pkg{PROcess} 
    package.}
  \item{span}{two moving window widths. Smaller one will be used for smoothing 
    and local maxima finding. Larger one will be used for local variance 
    estimation. Similar to \code{span} and \code{sm.span} variables in 
    \code{\link[PROcess]{isPeak}} from \pkg{PROcess} package.}
  \item{zerothresh}{Intensity threshold criterion for peak detection. Positive 
    numbers in range [0,1), 
    like default 0.9, will be used to calculate a single threshold used 
    for all samples using \code{quantile(X,zerothresh)} equation. Negative 
    numbers in range (-1, 0) will be used to calculate threshold for each single 
    sample \code{i} using \code{quantile(X[i,],-zerothresh)}. 
    Similar to \code{zerothrsh} variable in 
    \code{\link[PROcess]{isPeak}} from \pkg{PROcess} package. }
}

\details{
  Peak finding is done using the following algorithm:   
  \tabular{lcl}{
		\code{x}      \tab \code{=} \tab \code{ X[j,] } \cr
    \code{thresh} \tab \code{=} \tab \code{ if(zerothresh>=0) quantile(X,zerothresh) else quantile(x,-zerothresh) } \cr  
    \code{sig}    \tab \code{=} \tab \code{ runmean(x, span[2]) } \cr              
    \code{rMax}   \tab \code{=} \tab \code{ runmax (x, span[2]) } \cr                
    \code{rAvr}   \tab \code{=} \tab \code{ runmed (x, span[1]) } \cr
    \code{rStd}   \tab \code{=} \tab \code{ runmad (x, span[1], center=rAvr) } \cr
    \code{peak}   \tab \code{=} \tab \code{ (rMax == x) & (sig > thresh) & (sig-rAvr > SNR*rStd) }
  }
  
  What means that a peak have to meet the following criteria to be classified 
  as a peak:
  \itemize{
    \item be a local maxima in \code{span[2]} neighborhood
    \item smoothed sample (\code{sig}) is above user defined threshold 
    \code{zerothresh}
    \item locally calculated z-score (see 
    \url{http://mathworld.wolfram.com/z-Score.html}) of the signal is above 
    user defined signal-to-noise ratio 
  }
  It is very similar to the \code{isPeak} and \code{getPeaks} functions from 
  \pkg{PROcess} library (ver 1.3.2) written by Xiaochun Li. For example 
  \code{getPeaks(X, PeakFile, SoN=SNR, span=span[1], sm.span=span[2], 
    zerothrsh=zerothresh, area.w=0.003, ratio=0) } would give very similar 
  results as \code{msc.peaks.find} the differences include: speed ( \code{msc.peaks.find} 
  uses much faster C-level code), different use of signal-to-noise-ratio 
  variable, and \code{msc.peaks.find} does not do or use area calculations.
}

\value{
  A data frame, in the same format as data saved in \code{peakinfofile}, have 
  five components: 
  \item{Spectrum.Tag}{sample name of each peak}
  \item{Spectrum.}{sample number of each peak}
  \item{Peak.}{peak number within each sample}
  \item{Intensity}{peak height (intensity)}
  \item{Substance.Mass}{x-axis position, or corresponding mass of 
  the peak measured in M/Z, which were extracted from row names of the 
  \code{X} matrix.}
}

\author{Jarek Tuszynski (SAIC) \email{jaroslaw.w.tuszynski@saic.com}} 

\seealso{
  \itemize{
  \item Part of \code{\link{msc.preprocess.run}} and 
    \code{\link{msc.project.run}} pipelines.
  \item Previous step in the pipeline was \code{\link{msc.mass.adjust}} 
  \item Functions \code{\link{msc.peaks.align}} or \code{\link[PROcess]{pk2bmkr}} 
    can be used to align peaks from different samples in order to find 
    biomarkers. 
  \item Peak data can be read and writen by \code{\link{msc.peaks.read.csv}} 
        and \code{\link{msc.peaks.write.csv}}. 
  \item Functions \code{\link[PROcess]{isPeak}} and 
  \code{\link[PROcess]{getPeaks}} from \pkg{PROcess} package are very similar.
  \item Uses \code{\link{runmax}}, \code{\link{runmean}}, 
    \code{\link{runmed}}, \code{\link{runmad}} functions.
  }
}

\examples{
  # load input data 
  if (!file.exists("Data_IMAC.Rdata")) example("msc.project.read")
  load("Data_IMAC.Rdata")

  # Find Peaks
  Peaks = msc.peaks.find(X)
  cat(nrow(Peaks), "peaks were found in", Peaks[nrow(Peaks),2], "files.\n")
}

\keyword{ts}

