\name{msc.msfiles.read.csv}
\alias{msc.msfiles.read.csv}
\title{Read Protein Mass Spectra from CSV files}

\description{
  Read multiple protein mass spectra (SELDI) files, listed in \code{FileList}, 
  from a given directory and combine them into a single data structure. Files 
  are in CSV format, possibly compresses. If \code{FileList} is
  an 1D list than data is stored as a matrix one file per column. If
  \code{FileList} is a 2D data-frame than data is stored in 3D array.
}

\usage{
  msc.msfiles.read.csv(directory=".", FileList="\\.csv", 
                       SampleNames=NULL, CopyNames=NULL)
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
       \item data.frame (multiple lists of file names)- multiple copies of the 
       same samples are present - see details 
     }
     The last two formats also support file zip and gzip file compression. 
     For example if individual file name is in the format:
     \itemize{
       \item "dir/a.csv" - uncompressed file 'a.csv' in directory 'dir'
       \item "dir/b.zip/a.csv" - file 'a.csv' within zipped file 'b.zip'
       \item "dir/a.csv.gz" - gziped individual file
     }
  } 
  \item{SampleNames}{Optional list of names to be used as sample/column names.}
  \item{CopyNames}{Optional list of names to be used as copy/plane names in 
    case FileList is an 2D data frame.}
}

\details{
  All files should be in Excel's CSV format (table in text format: 1 row per 
  line, comma delaminated columns). Each file is assumed to have two columns, in
  case of SELDI data: column 1 (x-axis) is mass/charge (M/Z), and column 2 
  (y-axis) is spectrum intensity.  All files are assumed to have identical 
  first (M/Z) column. 
  
  If multiple copies of the same sample were collected than one can store them 
  in a 3D array (data cube) where each column correspond to a single sample,
  each row is a single mass (M/Z) and each plane is a single copy. To do so one
  has to pass a 2D data frame as \code{FileList} where each column contains file
  names of multiple copies of the same sample and each row contains filenames of 
  a single copy of different samples.
}


\value{Data structure containing all the data read from the files. It can be in 
  form of a 2D matrix (nFeatures \eqn{\times}{x} nSamples) or 3D array
  (nFeatures \eqn{\times}{x} nSamples \eqn{\times}{x} nCopies) depending on input.} 

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
  X = msc.msfiles.read.csv(directory, "IMAC_normal_.*csv")
  dim(X)
  
  # example of explicite 1D FileList
  ProjectFile = file.path(directory,"InputFiles.csv")
  FileList = read.csv(file=ProjectFile, comment.char = "")
  FileList[,3]
  X = msc.msfiles.read.csv(directory, FileList=FileList[,3], SampleNames=FileList[,1])
  dim(X)
  
  # example of explicite 2D FileList
  FileList[,3:4]
  X = msc.msfiles.read.csv(directory, FileList=FileList[,3:4], 
        SampleNames=FileList[,1], CopyNames=c("copy 1", "copy 2"))
  dim(X)
}

\keyword{ts}
\keyword{file}
