\name{msc.project.read}
\alias{msc.project.read}
\title{Read and Manage a Batch of Protein Mass Spectra}
\description{Read and manage a batch of protein mass spectra (SELDI) files 
  where files could contain multiple spectra taken from the same sample, or 
  multiple experiments performed on the same sample.}
\usage{msc.project.read(ProjectFile, directory.out=NULL) }

\arguments{
  \item{ProjectFile}{path and name of text file in Excel's CSV format which is 
    used to 
    store information about a batch of Mass Spectra data files.  See details. } 
  \item{directory.out}{Optional character vector with name of directory where 
     output files will be saved. Use "/" slashes 
     in directory name. By default the directory containing \code{ProjectFile}
     and all Mass Spectra files is used, and this argument is provided in case 
     that directory is read-only and user have to choose a different directory.}
}

\details{
  Function \code{msc.project.read} allows to user to manage large batches of
  Mass Spectra 
  files, especially when multiple copies of each sample are present. The 
  \code{ProjectFile} contains all the information about the project. An example 
  format might be:
 
  \tabular{lcllll}{
  \bold{Name}, \tab \bold{Class}, \tab \bold{IMAC1}, \tab \bold{IMAC2}, \tab \bold{WCX1},	\tab \bold{WCX2} \cr
  r0008, \tab 1, \tab Nr/imac\_r0008.csv, \tab Nr/imac\_r0008(2).csv, \tab Nr/wcx\_r0008.csv, \tab Nr/wcx\_r0008(2).csv \cr
  r0012, \tab 1, \tab Nr/imac\_r0012.csv, \tab Nr/imac\_r0012(2).csv, \tab Nr/wcx\_r0012.csv, \tab Nr/wcx\_r0012(2).csv \cr
  r0014, \tab 1, \tab Nr/imac\_r0014.csv, \tab Nr/imac\_r0014(2).csv, \tab Nr/wcx\_r0014.csv, \tab Nr/wcx\_r0014(2).csv \cr
  r0021, \tab 2, \tab Ca/imac\_r0021.csv, \tab Ca/imac\_r0021(2).csv, \tab Ca/wcx\_r0021.csv, \tab Ca/wcx\_r0021(2).csv \cr
  r0022, \tab 2, \tab Ca/imac\_r0022.csv, \tab Ca/imac\_r0022(2).csv, \tab Ca/wcx\_r0022.csv, \tab Ca/wcx\_r0022(2).csv \cr
  r0024, \tab 2, \tab Ca/imac\_r0024.csv, \tab Ca/imac\_r0024(2).csv, \tab Ca/wcx\_r0024.csv, \tab Ca/wcx\_r0024(2).csv \cr
  r0027, \tab 2, \tab Ca/imac\_r0027.csv, \tab Ca/imac\_r0027(2).csv, \tab Ca/wcx\_r0027.csv, \tab Ca/wcx\_r0027(2).csv
  }
  
  \code{ProjectFile} always has the following format:
  \itemize{
    \item column 1 - unique name for each sample - Those names will be used in 
    the program to identify the samples
    \item column 2 - class label for each sample - in the classification part 
    of the code those labels will be used as a response vector (target values). 
    Usually a factor for classification, but could be a unique number for 
    regression.
    \item columns 3+ - file path (from \code{directory}) for each file in the 
    project. If \code{ProjectFile} has more than 3 columns than multiple copies 
    of the same sample are present. In that case column labels (IMAC1, IMAC2, 
    WCX1, WCX2) become important, since they distinguish between equivalent 
    copies taken under the same conditions and copies taken under different 
    conditions. In our example both kinds of copies exist: files in columns 
    IMAC1 and IMAC2 contain two copies of spectra collected using Ciphergen's  
    IMAC ProteinChip array and files in columns WCX1 and WCX2 used WCX array.  
    The labels of those columns are expected to use letters as labels for 
    different copies and numbers to mark multiple identical copies.
  }
  
  File names in \code{ProjectFile} could be compressed uzing zip and gzip 
  file compression. For example if individual file name is in the format:
   \itemize{
     \item "a.csv" - trivial case - uncompressed file
     \item "b.zip/a.csv" - file within zipped file
     \item "a.csv.gz" - gziped individual file
   }
}


\value{
  List of .Rdata files storing data that was just read. Each file contains 
  either 2D data (if only one copy of the the data existed) or 3D data (if 
  multiple copies of the data existed). Multiple files are produced if multiple
  experiments were performed under different conditions. In above example two
  files will be produced: Data\_IMAC.Rdata and Data\_WCX.Rdata.
} 

\author{Jarek Tuszynski (SAIC) \email{jaroslaw.w.tuszynski@saic.com}} 

\seealso{
  \itemize{
  \item \code{\link{msc.msfiles.read.csv}} is a lower level function useful for data 
  that does not have multiple copies,  
  \item \code{\link{msc.preprocess.run}} is usually used to process output of 
  this function.
  \item \code{\link[PROcess]{read.files}} from \pkg{PROcess} library can 
  read a single SELDI file and \code{\link[PROcess]{rmBaseline}} can read in a 
  directory of files and substract their baselines.
  \item \code{\link[ppc]{ppc.read.raw.batch}} and 
  \code{\link[ppc]{ppc.read.raw.nobatch}} from \pkg{ppc} library can 
  also read SELDI files, assuming correct directory structure.
  }
}

\examples{
  directory = system.file("Test", package = "caMassClass")
  ProjectFile = file.path(directory,"InputFiles.csv")
  FileNames = msc.project.read(ProjectFile, '.')
  cat("File ",FileNames," was created\n")
}

\keyword{ts}
\keyword{file}
