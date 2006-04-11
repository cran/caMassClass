\name{msc.preprocess.run}
\alias{msc.preprocess.run}
\title{Preprocessing Pipeline of Protein Mass Spectra }
\description{Pipeline for preprocessing protein mass spectra (SELDI) data 
  before classification. }
\usage{
msc.preprocess.run ( X, mzXML=NULL,
    baseline.removal = 0,
      breaks=200, qntl=0, bw=0.005,                    # bslnoff
    min.mass = 3000,                                   # msc.mass.cut
    mass.drift.adjustment = 1,
      shiftPar=0.0005,                                 # msc.mass.adjust
    peak.extraction = 0, 
     PeakFile=0, SNR=2, span=c(81,11), zerothresh=0.9, # msc.peaks.find
     BmrkFile=0, BinSize=c(0.002, 0.008), tol=0.97,    # msc.peaks.align 
     FlBmFile=0, FillType=0.9,                         # msc.biomarkers.fill
    merge.copies = 1,                                  # msc.copies.merge
    verbose = TRUE) 
}

\arguments{
  \item{X}{Spectrum data either in matrix format [nFeatures \eqn{\times}{x} nSamples] or in 
    3D array format [nFeatures \eqn{\times}{x} nSamples \eqn{\times}{x} nCopies]. Row names 
    \code{(rownames(X)} store M/Z mass of each row. }
  \item{mzXML}{optional record of experimental setup and processing so far, 
    beeing prepared for possible output as mzXML file.}
  \item{baseline.removal}{Remove baseline from each spectrum? (boolean or 0/1 
    integer). See function \code{\link{msc.baseline.subtract}} and
    \code{\link[PROcess]{bslnoff}} 
    from \pkg{PROcess} library for other parameters that can be passed:
    \bold{\code{breaks}}, \bold{\code{qntl}} and \bold{\code{bw}}. }
  \item{min.mass}{Cutting place when removing data corresponding to low masses 
    (m/z). See function \code{\link{msc.mass.cut}} for details.}
  \item{mass.drift.adjustment}{Controls mass drift adjustment and scaling. 
    If 0 than no mass adjustment or scaling will be performed; otherwise, it is
    passed to \code{\link{msc.mass.adjust}} function as \code{scalePar}. Because 
    of that: 1 means that afterwards all samples will have the same mean, 2 
    means that afterwards all samples will have the same mean and medium. See 
    function \code{\link{msc.mass.adjust}} for details and additional parameter 
    \bold{\code{shiftPar}} that can be passed.}
  \item{peak.extraction}{Perform peak extraction and alignment, or keep on
    working with the raw spectra? (boolean or 0/1 integer). See following 
    functions for other parameters that can be passed: 
    \itemize{
      \item \code{\link{msc.peaks.find}} - see parameters: \bold{\code{PeakFile}}, 
      \bold{\code{SNR}}, \bold{\code{span}} and \bold{\code{zerothresh}} 
      \item \code{\link{msc.peaks.align}} - see parameters: \bold{\code{BmrkFile}}, 
      \bold{\code{BinSize}} and \bold{\code{tol}} 
      \item \code{\link{msc.biomarkers.fill}} - see parameters: 
      \bold{\code{FlBmFile}} and \bold{\code{FillType}}   
    }
    Especially filenames to store intermediate results. }
  \item{PeakFile}{Optional filename, storing peak finding results. If provided 
    than CSV file will be created 
     in the same format as Ciphergen's peak-info file, with following columns 
    of data: "Spectrum.Tag", "Spectrum.", "Peak.", "Intensity" and 
    "Substance.Mass". }
  \item{BmrkFile}{Optional filename, storing peak alignment results. If provided 
    than CSV file will be created 
    in the same format as Ciphergen's biomarker file, with spectra (samples) 
  	as rows, and biomarkers as columns (features).}
  \item{FlBmFile}{Optional filename, storing results of 
    \code{\link{msc.biomarkers.fill}}. If provided than CSV file will be created 
    in the same format as Ciphergen's biomarker file, with spectra (samples) 
  	as rows, and biomarkers as columns.}

  \item{merge.copies}{In case multiple copies of data exist should they be
    merged and how? Passed to \code{\link{msc.copies.merge}} function as 
    \code{mergeType} variable. See that function for more details. }
  \item{verbose}{Boolean flag turns debugging printouts on.} 
  \item{breaks}{parameter to be passed to \code{\link[PROcess]{bslnoff}} 
    function from \pkg{PROcess} library by \code{\link{msc.baseline.subtract}}  }
  \item{qntl}{parameter to be passed to \code{\link[PROcess]{bslnoff}} 
    function from \pkg{PROcess} library by \code{\link{msc.baseline.subtract}}  }
  \item{bw}{parameter to be passed to \code{\link[PROcess]{bslnoff}} 
    function from \pkg{PROcess} library by \code{\link{msc.baseline.subtract}}  }
  \item{shiftPar}{parameter to be passed to \code{\link{msc.mass.adjust}}  }
  \item{SNR}{parameter to be passed to \code{\link{msc.peaks.find}}  }
  \item{span}{parameter to be passed to \code{\link{msc.peaks.find}}  }
  \item{zerothresh}{parameter to be passed to \code{\link{msc.peaks.find}} }
  \item{BinSize}{parameter to be passed to \code{\link{msc.peaks.align}}  }
  \item{tol}{parameter to be passed to \code{\link{msc.peaks.align}}  }
  \item{FillType}{parameter to be passed to \code{\link{msc.biomarkers.fill}}  }
}

\details{
  Function containing several pre-processing steps preparing protein mass 
  spectra (SELDI) data for classification. This function is a 
  "pipeline" performing several operations, all of which do not need class label
  information. Any and all steps are optional and can be skipped: 
  \itemize{
    \item Remove baseline from each spectrum, using 
      \code{\link{msc.baseline.subtract}} and \code{\link[PROcess]{bslnoff}} 
      from \pkg{PROcess} library.
    \item Remove data corresponding to low masses (m/z), using 
      \code{\link{msc.mass.cut}}.
    \item Adjust for mass drift and normalize data, using 
      \code{\link{msc.mass.adjust}}.
    \item Find peaks and align them into "biomarker" matrix, using 
      \code{\link{msc.peaks.find}}, \code{\link{msc.peaks.align}} and 
      \code{\link{msc.biomarkers.fill}}.
    \item Merge multiple copies of data, using \code{\link{msc.copies.merge}}.
  }
}

\value{Return matrix containing features as rows and samples as columns, unless 
  \code{merge.copies} was 0,4, or 8 when no merging is done and data is 
  returned in same or similar format as the input format 
  [nFeatures \eqn{\times}{x} nSamples \eqn{\times}{x} nCopies]. 
  Row names \code{(rownames(X)} store M/Z mass of each row. 
  If \code{mzXML} input argument was not null than updated version of mzXML 
  record will be outputted as "\code{mzXML}"attribute of \code{X}.
}

\author{Jarek Tuszynski (SAIC) \email{jaroslaw.w.tuszynski@saic.com}} 
\seealso{
  \itemize{
  \item Input data likely come from \code{\link{msc.project.read}} or 
    \code{\link{msc.rawMS.read.csv}} functions
  \item As mentioned above function uses the following lower level functions:
   \code{\link{msc.baseline.subtract}}, \code{\link[PROcess]{bslnoff}} 
   from \pkg{PROcess} library, \code{\link{msc.mass.cut}},
   \code{\link{msc.mass.adjust}}, \code{\link{msc.peaks.find}}, 
   \code{\link{msc.peaks.align}}, \code{\link{msc.biomarkers.fill}}, and
   \code{\link{msc.copies.merge}}.
  \item Output data can be latter used for classification using 
   \code{\link{msc.classifier.test}} function
  }
}

\examples{
  # load "Data_IMAC.Rdata" file containing raw MS spectra 'X'  
  if (!file.exists("Data_IMAC.Rdata")) example("msc.project.read")
  load("Data_IMAC.Rdata") # load data: X & mzXML
  
  # run preprocess with peak extraction
  Y = msc.preprocess.run(X, mzXML=mzXML, peak.extraction=1, 
    PeakFile="peaks_IMAC.mzXML", BmrkFile="bmrk_IMAC.csv")
  cat("Size before: ", dim(X), " and after :", dim(Y), "\n")
  stopifnot( dim(Y)==c(25, 40) )        # make sure it is what's expected
  YmzXML = attr(Y, "mzXML")
  strsplit(YmzXML$dataProcessing, '\n') # show mzXML$dataProcessing record
  # inspect by hand output files: "peaks_IMAC.mzXML" & "bmrk_IMAC.csv"
  
  # run preprocess with no peak extraction
  Y = msc.preprocess.run(X, mzXML=mzXML)
  cat("Size before: ", dim(X), " and after :", dim(Y), "\n")
  stopifnot( dim(Y)==c(9377, 40) )      # make sure it is what's expected
  YmzXML = attr(Y, "mzXML")
  strsplit(YmzXML$dataProcessing, '\n') # show mzXML$dataProcessing record
}

\keyword{ts}
