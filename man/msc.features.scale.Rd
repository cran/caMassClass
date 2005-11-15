\name{msc.features.scale}
\alias{msc.features.scale}
\title{Scale Classification Data}
\description{Scale features of the data to be used for classification. Scaling 
  factors are extracted from each column/feature of the train data-set and
  applied to both train and test sets.}
\usage{msc.features.scale( xtrain, xtest, type = c("min-max", "avr-std", "med-mad"))}

\arguments{
  \item{xtrain}{A matrix or data frame with train data. Rows contain samples 
    and columns contain features/variables}
  \item{xtest}{A matrix or data frame with test data. Rows contain samples 
    and columns contain features/variables}
    
  \item{type}{ Following types are recognized 
    \itemize{
       \item "min-max" - data minimum is mapped to 0 and maximum is mapped to 1
       \item "avr-std" - data is mapped to zero \code{\link{mean}} and unit variance 
       \item "med-mad" - data is mapped to zero \code{\link{median}} and unit 
       \code{\link{mad}} (median absolute deviation) 
     }
  }
}

\details{
  Many classification algorithms perform better if input data is scaled 
  beforehand. Some of them perform scaling internally (for example 
  \code{\link[e1071]{svm}}), but many don't. For some it makes no difference
  (for example \code{\link[rpart]{rpart}} or \code{\link[caTools]{LogitBoost}}).
  
  In case \code{xtrain} contains \code{NA} values or infinities all non-finite 
  numbers are omited from scaling parameter calculations.
}

\value{ 
  \item{xtrain}{A matrix or data frame with scaled train data.}
  \item{xtest}{A matrix or data frame with scaled test data. }
} 

\author{Jarek Tuszynski (SAIC) \email{jaroslaw.w.tuszynski@saic.com}} 

\seealso{
  Used by \code{\link{msc.classifier.test}} and \code{\link{msc.features.select}} 
  functions.
}

\examples{
  library(e1071)
  data(iris)
  mask   = sample.split(iris[,5], SplitRatio=1/4) # very few points to train
  xtrain = iris[ mask,-5]  # use output of sample.split to ...
  xtest  = iris[!mask,-5]  # create train and test subsets
  ytrain = iris[ mask, 5] 
  ytest  = iris[!mask, 5] 
  x = msc.features.scale(xtrain, xtest)  
  model = svm(x$xtrain, ytrain, scale=FALSE)
  print(a <- table(predict(model, x$xtest), ytest) )
  model = svm(xtrain, ytrain, scale=FALSE)
  print(b <- table(predict(model, xtest), ytest) )
  stopifnot( sum(diag(a))<sum(diag(b)) )
}

\keyword{classif}

