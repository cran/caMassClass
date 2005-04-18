\name{msc.copies.merge}
\alias{msc.copies.merge}
\title{Merge Multiple Copies of Mass Spectra Samples}
\description{Protein mass spectra (SELDI) samples are sometimes scanned 
  multiple times in order to reduce hardware or software based errors. 
  \code{msc.copies.merge} function is used to merge, concatenate, and/or average all
  of those copies together in preparation for classification. }
\usage{msc.copies.merge( X, mergeType, PeaksOnly=TRUE) }

\arguments{
  \item{X}{Spectrum data in 3D array format [nFeatures \eqn{\times}{x} nSamples \eqn{\times}{x} nCopies]. 
    Row names (\code{rownames(X)}) store M/Z mass of each row. 
    If X is in matrix format [nFeatures \eqn{\times}{x} nSamples] nothing will be done. }
  \item{mergeType}{ an integer variable in [0,11] range, telling how to merge 
    samples and what to do with bad copies:
    \itemize{
       \item 0 - do nothing
       \item add 1 - if all original copies are to be concatenated as separate 
       samples
       \item add 2 - if copies are to be averaged and the average added as a 
       separate sample
       \item add 4 - if for each sample the worst copy is to be deleted
       \item add 8 - if for each sample in case of large differences between 
       copies, a single bad copy of a sample is to be replaced with the best 
       copy. Not to be used with previous option. See details.
    }
  }
  \item{PeaksOnly}{This variable is being passed to function 
    \code{\link{msc.sample.correlation}}. Set it to \code{TRUE} in case of raw 
    spectra and switch to \code{FALSE} in case of data where only peaks 
    (biomarkers) are present.}
}

\details{  
   Quality of a sample is 
   measured by calculating for each copy of each sample two variables:
   inner correlation (average correlation between multiple copies of the same 
   sample) and outer correlation (average correlation between each sample and 
   every other sample within the same copy). Inner correlation measures how 
   similar copies are to each other and outer correlation measures how similar 
   each copy is to everybody else. For example in case of experiment using 
   SELDI technology to distinguish cancerous samples and non-cancerous samples 
   one can assume that most of the proteins present in both cancerous and 
   non-cancerous samples will be the same. In that case one will expect high 
   correlation between samples and even higher correlation between copies of 
   the same sample
   
   if \code{mergeType/4} (\code{mergeType \%/\% 4})  is 
   \itemize{
     \item 0 - all copies are kept
     \item 1 - if inner correlation is smaller than outer correlation, or in 
     other words, if a signature is more similar to other signatures than to 
     other copies of the same signature, than there is some problem with that 
     signature. In that case that bad signature can be replaced with the best 
     copy of the signature. 
     \item 2 - rate each copy of each sample using 
     \code{score=outer_correlation + inner_correlation} measure. 
     Delete worst copy. 
   }
   Option 2 is more suitable in case of data with a lot of copies, when we can 
   afford dropping one copy. Option 1 is designed to patch the most serious problems 
   with the data.
   
   There are also four merging options, if \code{mergeType mod 4} 
   (\code{mergeType \%\% 4})  is 
   \itemize{
     \item 0 - no merging is done to the data and it is left as 3D array
     \item 1 - all copies are concatenated 
       \code{X = cbind(X[,,1], X[,,2], \dots, X[,,nCopy])} 
       so they seem as separate samples
     \item 2 - all copies are averaged 
       \code{X = (X[,,1] + X[,,2] + \dots + X[,,nCopy])/nCopy)}
     \item 3 - all copies are first averaged and than concatenated with extra 
     average copy 
       \code{X = cbind(X[,,1], X[,,2], \dots, X[,,nCopy], Xavr)}
   }
   In preparation for classification one can use multiple copies in several 
   ways: option 2 above improves (one hopes) accuracy of each sample, while 
   options 1 and 3 increase number of samples available during classification.  
   So the choice is: do we want a lot of samples during classification or 
   fewer, better samples?
   
   The best option of \code{mergeType} depends on kind of data.      
   \itemize{
     \item 0 if data has single copy.
     \item 1+2+4 will produce the largest number of samples since we will keep 
     all the copies and an average of all the copies 
     \item 2+8 will produce single most accurate sample from multiple copies 
     (usually if more than 2 copies are present) since we will delete outliers 
     before averaging all the copies
   }
}


\value{Return matrix containing features as rows and samples as columns, unless 
  \code{mergeType} is 0,4, or 8 when no merging is done and data is returned in 
  same or similar format as the input format.} 
\author{Jarek Tuszynski (SAIC) \email{jaroslaw.w.tuszynski@saic.com}} 
\seealso{
  \itemize{
  \item Part of \code{\link{msc.preprocess.run}} and 
   \code{\link{msc.project.run}} pipelines.
  \item Previous step in the pipeline was \code{\link{msc.mass.adjust}} or
  peak finding functions: \code{\link{msc.peaks.find}}, \code{\link{msc.peaks.align}}, and
  \code{\link{msc.biomarkers.fill}} 
  \item Next step in the pipeline is data classification
   \code{\link{msc.classifier.test}}
  \item Uses \code{\link{msc.sample.correlation}} 
  }
}

\examples{
  # load input data 
  if (!file.exists("Data_IMAC.Rdata")) example("msc.project.read")
  load("Data_IMAC.Rdata")
  
  # run msc.copies.merge
  Y = msc.copies.merge(X, 1+2+4)
  colnames(Y)
}

\keyword{ts}
