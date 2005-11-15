\name{msc.rawMS.read.csv}
\alias{msc.rawMS.read.csv}
\title{Read Protein Mass Spectra from CSV files}

\description{
  Read multiple protein mass spectra (SELDI) files, listed in \code{FileList}, 
  from a given directory and combine them into a single data structure. Files 
  are in CSV format, possibly compresses. Data is stored as a matrix one file 
  per column. 
}

\usage{
  msc.rawMS.read.csv(directory=".", FileList="\\.csv", mzXML.record=FALSE)
}

\arguments{
  \item{directory}{a character vector with name of directory where all the 
    files can be found. Use "/" slashes in directory name. The default 
    corresponds to the working directory \code{\link{getwd}}().}
  \item{FileList}{List of files to read. List can be in the following formats: 
     \itemize{
       \item single string - a regular expression (see \code{\link{regex}}) to 
       be used in selecting files to read, for example "\\.csv"
       \item list - list of file names to be read
      }
     The last format also support file zip and gzip file compression. 
     For example if individual file name is in the format:
     \itemize{
       \item "dir/a.csv" - uncompressed file 'a.csv' in directory 'dir'
       \item "dir/b.zip/a.csv" - file 'a.csv' within zipped file 'b.zip'
       \item "dir/a.csv.gz" - gziped individual file
     }
  } 
  \item{mzXML.record}{should mzXML record be created to store mata-data (input 
    file names)? }
}

\details{
  All files should be in Excel's CSV format (table in text format: 1 row per 
  line, comma delaminated columns). Each file is assumed to have two columns, in
  case of SELDI data: column 1 (x-axis) is mass/charge (M/Z), and column 2 
  (y-axis) is spectrum intensity.  All files are assumed to have identical 
  first (M/Z) column. 
}


\value{Data structure containing all the data read from the files, in 
  form of a 2D matrix (nFeatures \eqn{\times}{x} nSamples). If 
  \code{mzXML.record} was set to true than mzXML record with input file names
  will be attached to \code{X} as "\code{mzXML}" attribute.
} 

\author{Jarek Tuszynski (SAIC) \email{jaroslaw.w.tuszynski@saic.com}} 

\seealso{
  \itemize{
  \item Part of  \code{\link{msc.project.run}} pipeline.
  \item \code{\link{msc.project.read}} gives user much more flexibility in 
  defining the meaning of the data to be read.
  \item \code{\link{msc.preprocess.run}} is often used as a next step in 
  the process
  \item \code{\link[PROcess]{read.files}} from \pkg{PROcess} library can 
  read a single SELDI file and \code{\link[PROcess]{rmBaseline}} can read in a 
  directory of files and substract their baselines.
  \item \code{\link[ppc]{ppc.read.raw.batch}} and 
  \code{\link[ppc]{ppc.read.raw.nobatch}} from \pkg{ppc} library can 
  also read SELDI files, assuming correct directory structure.
  }
}

\examples{
  # example of mode "single string" FileList
  directory  = system.file("Test", package = "caMassClass")
  X = msc.rawMS.read.csv(directory, "IMAC_normal_.*csv")
  stopifnot ( dim(X) == c(11883, 20) ) # make sure it is what's expected
  
  # example of explicit 1D FileList
  ProjectFile = file.path(directory,"InputFiles.csv")
  FileList = read.csv(file=ProjectFile, comment.char = "")
  FileList[,3]
  X = msc.rawMS.read.csv(directory, FileList=FileList[,3], mzXML.record=TRUE)
  stopifnot ( dim(X) == c(11883, 20) ) # make sure it is what's expected
  mzXML = attr(X,"mzXML")
  strsplit(mzXML$parentFile, '\n')      # show mzXML$parentFile record
  
  # example using data provided in PROcess package 
  directory  = system.file("Test", package = "PROcess")
  X = msc.rawMS.read.csv(directory) 
  msc.baseline.subtract(X, plot=TRUE) # used here to plot results
  dim(X)
}

\keyword{ts}
\keyword{file}
