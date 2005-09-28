\name{msc.classifier.run}
\alias{msc.classifier.run}
\title{Train and Test Chosen Classifier.}
\description{Common interface for training and testing several standard 
  classifiers. Includes feature selection and feature scaling steps. Allows to 
  specify that some test samples are multiple copies of the same sample, and 
  should return the same label.}
\usage{
msc.classifier.run( xtrain, ytrain, xtest, ret.prob=FALSE, 
    RemCorrCol=0, KeepCol=0, prior=1, same.sample=NULL,
    ScaleType=c("none", "min-max", "avr-std", "med-mad"),
    method=c("svm", "nnet", "lda", "qda", "LogitBoost", "rpart"), \dots) 
}

\arguments{
  \item{xtrain}{A matrix or data frame with training data. Rows contain samples 
    and columns contain features/variables}
  \item{ytrain}{Class labels for the training data samples. 
    A response vector with one label for each row/component of \code{x}.
    Can be either a factor, string or a numeric vector. }
  \item{xtest}{A matrix or data frame with test data. Rows contain samples 
    and columns contain features/variables}
  \item{ret.prob}{if set to TRUE than the a-posterior probabilities for 
   each class are returned as attribute called "probabilities".}
  \item{same.sample}{optional parameter which allows to specify that some 
    (or all) test samples have multiple copies which should be used to predict a
    single label for all of them. Can be either a factor, string or a numeric 
    vector, with unique values for different samples and identical values for 
    copies of the same sample.}
  \item{RemCorrCol}{If non-zero than some of the highly correlated columns are 
    removed using \code{\link{msc.features.remove}} function with 
    \code{ccMin=RemCorrCol}.}
  \item{KeepCol}{If non-zero than columns with low AUC are removed. 
     \itemize{
       \item if \code{KeepCol} smaller than 0.5 - do nothing
       \item if \code{KeepCol} in between [0.5, 1] - keep columns with AUC 
       bigger than \code{KeepCol} 
       \item if \code{KeepCol} bigger than one - keep top "\code{KeepCol}" 
       number of columns
     }
  }
  \item{ScaleType}{Optional parameter, if provided than following types are 
    recognized 
    \itemize{
       \item "none"    - no scaling is performed
       \item "min-max" - data minimum is mapped to 0 and maximum is mapped to 1
       \item "avr-std" - data is mapped to zero mean and unit variance 
       \item "med-mad" - data is mapped to zero median and unit mad (median 
       absolute deviation) 
     }
  }
  
  \item{prior}{class weights. following types are recognized 
    \itemize{ 
      \item \code{prior==1} - all samples in all classes have equal weight 
      (default)
      \item \code{prior==2} - all classes have equal weight
      \item \code{prior} is a vector - a named vector of weights for the 
      different classes, used for asymetric class sizes.
    }
  } 
  
  \item{method}{classifier to be used. Following ones are recognized (followed 
    by some parameters that could be passed through \ldots :
    \itemize{
      \item "svm" - see \code{\link[e1071]{svm}} from \pkg{e1071} package. 
        Possible parameters: \code{cost}, \code{gamma}
       \item "nnet" - see \code{\link[nnet]{nnet}} from \pkg{nnet} package. 
        Possible parameters: \code{size}, \code{decay}, \code{maxit}
       \item "LogitBoost" - see \code{\link[caTools]{LogitBoost}} from 
        \pkg{caTools} package
        Possible parameter: \code{nIter}
       \item "lda" - see \code{\link[MASS]{lda}} from \pkg{MASS} package. 
        Possible parameters: \code{method}
       \item "qda" - see \code{\link[MASS]{qda}} from \pkg{MASS} package. 
        Possible parameters: \code{method}
       \item "rpart" - see \code{\link[rpart]{rpart}} from \pkg{rpart} package. 
        Possible parameters: \code{minsplit}, \code{cp}, \code{maxdepth}
     }
  }
  
  \item{\ldots}{Additional parameters to be passed to classifiers. See 
    \code{method} for suggestions. } 
}

\details{
 This function performs the following steps:
 \itemize{
   \item Remove highly correlated columns and columns with low AUC with 
     \code{\link{msc.features.select}} function
   \item Scale each feature separatly using 
     \code{\link{msc.features.scale}} function
   \item Train chosen classifier using \code{xtrain} and \code{ytrain}
   \item Predict labels of \code{xtest} using trained model
   \item If \code{same.sample} variable is given than synchronise predicted 
    labels in such a way that all copies of the same sample return the same 
    label.
    \item Return labels. If \code{ret.prob=TRUE} then return a-posterior 
    probabilities as well.
 }
}

\value{ 
  Predicted class labels for each sample in \code{xtest}. 
  If \code{ret.prob=TRUE} than the a-posterior probabilities of each sample 
  belonging to each class are returned as attribute called "probabilities".
  The returned probabilities do not take into account \code{same.sample} 
  variable, used to synchronize predicted labels.   
} 

\author{Jarek Tuszynski (SAIC) \email{jaroslaw.w.tuszynski@saic.com}} 
\note{This function is not fully tested and might be changed in future versions}

\references{
  \itemize{
   \item "A Practical Guide to Support Vector Classification"
   by Chih-Wei Hsu, Chih-Chung Chang, and Chih-Jen Lin
   (http://www.csie.ntu.edu.tw/~cjlin/papers/guide/guide.pdf)
  }
}

\seealso{
  \itemize{
  \item Used by \code{\link{msc.classifier.test}} function.
  \item Best classifier parameter set can be found by 
   \code{\link[e1071]{tune}} function from \pkg{e1071} package.
  \item Uses \code{\link{msc.features.select}} and 
   \code{\link{msc.features.scale}} functions.
  \item Uses variety of classification algorithms: \code{\link[e1071]{svm}},
   \code{\link[nnet]{nnet}}, \code{\link[caTools]{LogitBoost}}, \code{\link[MASS]{lda}},
   \code{\link[MASS]{qda}}, \code{\link[rpart]{rpart}}
  }
}

\examples{
  data(iris)
  mask  = sample.split(iris[,5], SplitRatio=1/4) # very few points to train
  xtrain = iris[ mask,-5]  # use output of sample.split to ...
  xtest  = iris[!mask,-5]  # create train and test subsets
  ytrain = iris[ mask, 5] 
  ytest  = iris[!mask, 5] 
  table(ytrain, msc.classifier.run(xtrain,ytrain,xtrain, method="svm") )
  table(ytrain, msc.classifier.run(xtrain,ytrain,xtrain, method="LogitBoost") )
  table(ytrain, msc.classifier.run(xtrain,ytrain,xtrain, method="nnet") )
  table(ytrain, msc.classifier.run(xtrain,ytrain,xtrain, method="lda") )
  table(ytrain, msc.classifier.run(xtrain,ytrain,xtrain, method="qda") )
}

\keyword{classif}
