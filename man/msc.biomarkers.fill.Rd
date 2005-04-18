\name{msc.biomarkers.fill}
\alias{msc.biomarkers.fill}
\title{Fill Empty Spaces in Biomarker Matrix}
\description{Fill empty spaces (NA's) in biomarker matrix created by 
  \code{msc.peaks.align} }
\usage{msc.biomarkers.fill( X, Bmrks, BinBounds, FillType=0.9, BmrkFile=0)}

\arguments{
  \item{X}{Spectrum data either in matrix format [nFeatures \eqn{\times}{x} nSamples] or in 
    3D array format [nFeatures \eqn{\times}{x} nSamples \eqn{\times}{x} nCopies]. Row names 
    \code{(rownames(X))} store M/Z mass of each row. }
  \item{Bmrks}{biomarker matrix containing one sample per column and one 
    biomarker per row}
  \item{BinBounds}{position (mass) of left-most and right-most peak in each bin}
  \item{FillType}{how to fill empty spaces in biomarker data?
    \itemize{
      \item if \code{0<=FillType<=1} than fill spaces with 
      \code{quantile(probs=FillType)}. 
      For example: if \code{FillType=1/2} than medium will be used, 
      if \code{FillType=1} than maximum value will be used, 
      if \code{FillType=0.9} than maximum will be used after discarding 10\% of 
      "outliers"
      \item if \code{FillType<0} than empty spaces will not be filled and NA's
      will remain
      \item if \code{FillType==2} than X value closest to the center of the bin 
      will be used 
      \item if \code{FillType==3} empty spaces will be set to zero 
    }
  }
  \item{BmrkFile}{Optional filename. If provided than CSV file will be created 
    in the same format as Ciphergen's biomarker file, with spectra (samples) 
  	as rows, and biomarkers as columns.}
}

\details{
  This function attepts to correct a problem which is a side-effect of 
  \code{\link{msc.peaks.align}} function. Namely numerous NA's in biomarker data,
  each time when some peak was found only in some of the samples. 
  \code{\link{msc.peaks.align}} already removed the most problematic features using
  \code{SampFrac} variable, but likelly a lot of NA's remain and they can cause
  problem for some classification algorithms. 
}

\value{Data in the same format and size as \code{Bmrks}} 

\note{The whole idea of filling spaces in biomarker matrix is a little bit 
  suspect since we are mixing proverbial apples and oranges. However, it might 
  be better than the other options of filling empty spaces with zeros or 
  keeping NA's.} 
\author{Jarek Tuszynski (SAIC) \email{jaroslaw.w.tuszynski@saic.com}} 

\seealso{
  \itemize{
  \item Part of \code{\link{msc.preprocess.run}} and 
    \code{\link{msc.project.run}} pipelines.
  \item Input comes most likely from: \code{\link{msc.peaks.align}}, 
  or from Ciphergen's software
  \item Output can be processed further by \code{\link{msc.copies.merge}} 
  \item Biomarkers matrix can be read and writen by 
        \code{\link{msc.biomarkers.read.csv}} 
        and \code{\link{msc.biomarkers.write.csv}}. 
  \item Function \code{\link[PROcess]{pk2bmkr}} from \pkg{PROcess} package 
    lso perform similar function.
  }
}

\examples{
  # load input data 
  if (!file.exists("Data_IMAC.Rdata")) example("msc.project.read")
  load("Data_IMAC.Rdata")
  Y = msc.peaks.align(msc.peaks.find(X))
  dim(Y$Bmrks)
  print( Y$Bmrks , na.print=".",  digits=2)
  
  # run msc.biomarkers.fill
  Z = msc.biomarkers.fill( X, Y$Bmrks, Y$BinBounds)
  dim(Z)
  print( Z[,,1] , na.print=".",  digits=2)
  print( Z[,,2] , na.print=".",  digits=2)
  
  # run msc.biomarkers.fill with other FillType
  Z = msc.biomarkers.fill( X, Y$Bmrks, Y$BinBounds, FillType=2)
 }

\keyword{ts}
