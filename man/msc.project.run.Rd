\name{msc.project.run}
\alias{msc.project.run}
\title{Read and Preprocess Protein Mass Spectra}
\description{Read and preprocess protein mass spectra (SELDI) files
  where files could contain multiple copies of spectra taken from the same 
  sample, or spectra from multiple experiments performed on the same sample.}
\usage{msc.project.run(ProjectFile, directory.out=NULL, verbose = TRUE, \dots)}

\arguments{
  \item{ProjectFile}{path and name of text file in Excel's CSV format which is 
    used to 
    store information about a batch of Mass Spectra data files.  See details. } 
  \item{directory.out}{Optional character vector with name of directory where 
     output files will be saved. Use "/" slashes 
     in directory name. By default the directory containing \code{ProjectFile}
     and all Mass Spectra files is used, and this argument is provided in case 
     that directory is read-only and user have to choose a different directory.}
  \item{verbose}{boolean flag turns debugging printouts on.} 
  \item{\ldots}{parameters to be passed to \code{\link{msc.preprocess.run}}}
}

\details{
  High level processing of protein mass spectra (SELDI) data. 
  \code{\link{msc.project.read}} supports projects with multiple
    sets of spectra taken under different experimental condition. Those sets 
    will be called \dfn{batches}. With that in mind, following steps 
  are performed:
  \itemize{
    \item \code{\link{msc.project.read}(ProjectFile, directory)} is called 
    which reads and saves different batches of mass spectra (SELDI) data into 
    separate files. 
    List of those files is saved in temporary "RInputFiles.csv" file. 
    In future calls to \code{msc.project.run}, if above file exist than 
    \code{\link{msc.project.read}} is not called again.
    \item Each batch of data is \code{\link{load}}ed and preprocessed by
    calls to \code{\link{msc.preprocess.run}}. All the required parameters 
    have to be passed through "\dots" mechanism.
    \item In case of multiple batches of data results are \code{\link{rbind}}ed
  }
}

\value{ 
  \item{X}{Spectrum data either in matrix format [nFeatures \eqn{\times}{x} 
    nSamples].  Row names 
    \code{(rownames(X)} store M/Z mass of each row merged with batch name } 
   \item{SampleLabels}{Class label for each sample as read from 
    \code{\link{msc.project.read}} }
   \item{SameSample}{array of the same length as number of columns in \code{X}
     indicating samples that are multiple copies that came from the same 
     physical sample, and should be the same, if not for noise. }
   \item{mxXML}{Experiment information in mzXML file format. Included only if
     input data was read from mzXML files, and \code{NULL} otherwise. }
} 

\author{Jarek Tuszynski (SAIC) \email{jaroslaw.w.tuszynski@saic.com}} 
\seealso{\code{\link{msc.project.read}}, \code{\link{msc.preprocess.run}}}

\examples{
  directory  = system.file("Test", package = "caMassClass")
  ProjectFile = file.path(directory,"InputFiles.csv")
  Data = msc.project.run(ProjectFile, '.',  
      baseline.removal=0, mass.drift.adjustment=1, min.mass=3000, 
      peak.extraction=1, merge.copies=7, shiftPar=0.0004)
  stopifnot( dim(Data$X)==c(25,60) )
  strsplit(Data$mzXML$parentFile, '\n')      # show mzXML$parentFile record
  strsplit(Data$mzXML$dataProcessing, '\n')  # show mzXML$dataProcessing record
}

\keyword{ts}
