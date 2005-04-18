\name{msc.features.select}
\alias{msc.features.select}
\title{Reduce Number of Features Prior to Classification }
\description{Select subset of individual features that are potentially most 
  useful for classification.}
\usage{msc.features.select( x, y, RemCorrCol=0.98, KeepCol=0.6)}

\arguments{
  \item{x}{A matrix or data frame with training data. Rows contain samples 
    and columns contain features/variables}
  \item{y}{Class labels for the training data samples. 
    A response vector with one label for each row/component of \code{x}.
    Can be either a factor, string or a numeric vector.}
  \item{RemCorrCol}{If non-zero than some of the highly correlated columns are 
    removed using \code{\link{msc.features.remove}} function with 
    \code{ccMin=RemCorrCol}.}
  \item{KeepCol}{If non-zero than columns with low AUC are removed. 
     \itemize{
       \item if \code{KeepCol} smaller than 0.5 - do nothing
       \item if \code{KeepCol} in between [0.5, 1] - keep columns with AUC bigger 
       than \code{KeepCol} 
       \item if \code{KeepCol} bigger than one - keep top \code{KeepCol} 
       number of columns
     }
  }
}

\details{
  This function reduces number of features in the data prior to classification,
  using following steps:
  \itemize{
    \item calculate AUC measure for each feature using \code{\link{colAUC}}
    \item remove some of the highly correlated neighboring columns using 
    \code{\link{msc.features.remove}} function. 
    \item remove columns with low AUC
  }

  This function finds subset of individual features that are potentially most 
  useful for classification, and each feature is rated individually.
  However, often set of two or more 
  very poor individual features can produce a superior classifier.  So, this 
  function should be used with care. I found it very useful when classifying 
  raw protein mass spectra (SELDI) data, for reducing dimensionality of the 
  data from 10 000's to 100's prior of classification, instead of peak-finding 
  (see \code{\link{msc.peaks.find}}).
}

\value{Vector of column indexes to be kept.} 
\author{Jarek Tuszynski (SAIC) \email{jaroslaw.w.tuszynski@saic.com}
} 
\seealso{
  \itemize{
  \item Used by \code{\link{msc.classifier.test}} function.
  \item Uses \code{\link{colAUC}}, \code{\link{msc.features.remove}} and 
  \code{\link{msc.features.scale}} functions.
  }
}

\examples{
  # load input data 
  if (!file.exists("Data_IMAC.Rdata")) example("msc.project.read")
  load("Data_IMAC.Rdata")
  
  X = t(X[,,1])
  cidx = msc.features.select(X, SampleLabels, KeepCol=0.7)
  cat(length(cidx),"features were selected out of",ncol(X),"min(auc)=0.7\n")
  cidx = msc.features.select(X, SampleLabels, KeepCol=400)
  cat(length(cidx),"features were selected out of",ncol(X),"\n")
  cat(" min(auc)=", min(colAUC(X[,cidx], SampleLabels)),"\n")
  Y = X[,cidx]
}

\keyword{classif}
