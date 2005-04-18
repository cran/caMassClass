\name{msc.mass.adjust}
\alias{msc.mass.adjust}
\alias{msc.mass.adjust.apply}
\alias{msc.mass.adjust.calc}
\title{Perform Normalization and Mass Drift Adjustment for Mass 
  Spectra Data.}
\description{Perform normalization and mass drift adjustment for protein mass 
  spectra (SELDI) data.}
\usage{
  msc.mass.adjust(X, scalePar=2, shiftPar=0.0005, AvrSamp=0)
  msc.mass.adjust.calc(X, scalePar=2, shiftPar=0.0005, AvrSamp=0)
  msc.mass.adjust.apply(X, shiftX, scaleY, shiftY) 
}

\arguments{
  \item{X}{Spectrum data either in matrix format [nFeatures \eqn{\times}{x} nSamples] or in 
    3D array format [nFeatures \eqn{\times}{x} nSamples \eqn{\times}{x} nCopies]. Row names 
    \code{(rownames(X))} store M/Z mass of each row. }
  \item{scalePar}{Controls scaling (normalization): 1 means that afterwards all 
    samples will have the same mean, 2 means that afterwards all samples will 
    have the same mean and medium (default)}
  \item{shiftPar}{Controls mass adjustment. Shifting sample has to improve 
    correlation by at least that amount to be considered. Designed to prevent 
    shifts based on "improvement" on order of magnitude of machine accuracy. 
    If set to too large will turn off shifting. Default = 0.0005.}
  \item{AvrSamp}{Is used to normalize test set the same way train set was 
    normalized. Test set is processed using \code{AvrSamp} array that was one 
    of the outputs from train-set mass-adjustment. See examples.}
  \item{shiftX}{matrix [nSamp \eqn{\times}{x} nCopy] - integer number of positions a sample 
    should be shifted to the right (+) or left (-). Output from 
    \code{msc.mass.adjust.calc} and input to \code{msc.mass.adjust.apply}.}
  \item{scaleY}{matrix [nSamp \eqn{\times}{x} nCopy] - multiply each sample in order to 
    normalize it. Output from \code{msc.mass.adjust.calc} and input to 
    \code{msc.mass.adjust.apply}.}
  \item{shiftY}{matrix [nSamp \eqn{\times}{x} nCopy] - subtract this number from scaled 
    sample (if matching medians). Output from \code{msc.mass.adjust.calc} and 
    input to \code{msc.mass.adjust.apply}.}
}

\details{
  Mass adjustment assumes that SELDI data has some error associated with 
  inaccuracy of setting the starting point of time measurement (x-axis origin 
  or zero M/Z value). We try to correct this error by allowing the samples to 
  shift a few time-steps to the left or to the right, if that will help with 
  cross-correlation with other samples. The function performs the following 
  steps
   \itemize{
     \item normalize all samples in such a way as to make their means (and 
     optionally medians) the same
     \item if multiple copies exist than
       \itemize{
          \item align multiple copies of each sample to each other
          \item temporarily merge multiple copies of each sample to create a 
          "super-sample" vector with more features
        }
     \item align each sample to the mean of all samples
     \item recalculate mean of all samples and repeat above step
   }
  
   \code{msc.mass.adjust} function 
  was split into two parts (one to calculate parameters and one to apply them) 
  in order to give users more flexibility and information about what is done to 
  the data. This split allows inspection, plotting and/or modification of 
  \code{shiftX, shiftY, scaleY} parameters before data is modified. For example 
  one can set \code{shiftX} to zero to perform normalization without mass 
  adjustment or set \code{shiftY} to zero and \code{scaleY} to one to perform 
  mass adjustment without normalization. Three function provided are:
  \itemize{
    \item \code{msc.mass.adjust.calc} - calculates and returns all the 
    normalization and mass drift adjustment parameters
     \item \code{msc.mass.adjust.apply} - performs normalization 
     and mass drift adjustment using precalculated parameters
     \item \code{msc.mass.adjust} - simple interface version of above 2 functions
  }
}

\value{
  Functions \code{msc.mass.adjust} and \code{msc.mass.adjust.apply} return 
  modified spectra in the same format and size as \code{X}.
  Functions \code{msc.mass.adjust.calc} returns list containing the following:
  \item{shiftX}{matrix [nSamp \eqn{\times}{x} nCopy] - integer number of positions sample 
    should be shifted to the right (+) or left (-) }
  \item{scaleY}{matrix [nSamp \eqn{\times}{x} nCopy] - multiply each sample in order to 
    normalize it}
  \item{shiftY}{matrix [nSamp \eqn{\times}{x} nCopy] - subtract this number from scaled 
    sample (if matching mediums)}
  \item{AvrSamp}{Use AvrSamp returned from train-set mass-adjustment to 
    process test-set}
}

\author{Jarek Tuszynski (SAIC) \email{jaroslaw.w.tuszynski@saic.com}} 

\seealso{
  \itemize{
  \item Part of \code{\link{msc.preprocess.run}} and 
    \code{\link{msc.project.run}} pipelines.
  \item Previous step in the pipeline is \code{\link{msc.mass.cut}}
  \item Next step in the pipeline is either \code{\link{msc.peaks.find}} or 
    \code{\link{msc.copies.merge}}
  }
}

\examples{
  # load input data 
  if (!file.exists("Data_IMAC.Rdata")) example("msc.project.read")
  load("Data_IMAC.Rdata")
  
  # run on 3D input data using long syntax
  out = msc.mass.adjust.calc (X)
  Y   = msc.mass.adjust.apply(X, out$ShiftX, out$ScaleY, out$ShiftY)
  
  # check what happened to means
  Z   = cbind(colMeans(X), colMeans(Y))
  colnames(Z) = c("copy 1 before", "copy 2 before", "copy 1 after", "copy 2 after" )
  cat("Sample means after and after:\n")
  Z
  
  # check what happen to sample correlation
  A = msc.sample.correlation(X, PeaksOnly=TRUE)
  B = msc.sample.correlation(Y, PeaksOnly=TRUE)
  cat("Mean corelation between two copies of the same sample:\n")
  cat(" before: ", mean(A$innerCor)," after: ", mean(B$innerCor), "\n")
  cat("Mean corelation between unrelated samples:\n")
  cat(" before: ", mean(A$outerCor)," after: ", mean(B$outerCor), "\n")
  
  # run on 2D input data using short syntax 
  # check what happened to means and medians
  Y = msc.mass.adjust(X[,,1], scalePar=2)
  Z = cbind(colMeans(X[,,1]), apply(X[,,1],2,median), colMeans(Y), apply(Y,2,median))
  colnames(Z) = c("means before", "medians before", "means after", "medians after" )
  Z
  Y = msc.mass.adjust(X[,,1], scalePar=1)
  Z = cbind(colMeans(X[,,1]), apply(X[,,1],2,median), colMeans(Y), apply(Y,2,median))
  colnames(Z) = c("means before", "medians before", "means after", "medians after" )
  Z
  
  # mass adjustment for train and test sets, where test set is normalized in 
  # the same way as train set was
  Xtrain = X[, 1:10,]
  Xtest  = X[,11:20,]
  out    = msc.mass.adjust.calc (Xtrain);
  Xtrain = msc.mass.adjust.apply(Xtrain, out$ShiftX, out$ScaleY, out$ShiftY)
  out    = msc.mass.adjust.calc (Xtest , AvrSamp=out$AvrSamp);
  Xtest  = msc.mass.adjust.apply(Xtest , out$ShiftX, out$ScaleY, out$ShiftY)  
}

\keyword{ts}
