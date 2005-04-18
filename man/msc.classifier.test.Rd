\name{msc.classifier.test}
\alias{msc.classifier.test}
\title{Test a Classifier through Cross-validation}
\description{Test Classifier through Cross-validation. Common interface for 
  Cross-validation of several standard 
  classifiers. Includes feature selection and feature scaling steps. Allows to 
  specify that some test samples are multiple copies of the same sample, and 
  should return the same label.}
\usage{
 msc.classifier.test( X, Y, iters=50, SplitRatio=2/3, verbose=FALSE,
                  RemCorrCol=0, KeepCol=0, prior=1, same.sample=NULL,
                  ScaleType=c("none", "min-max", "avr-std", "med-mad"),
                  method=c("svm", "nnet", "lda", "qda", "LogitBoost", "rpart"), \dots) 
}

\arguments{
  \item{X}{A matrix or data frame with training data. Rows contain samples 
    and columns contain features/variables}
  \item{Y}{Class labels for the training data samples. 
    A response vector with one label for each row/component of \code{x}.
    Can be either a factor, string or a numeric vector. Labels with 'NA' value
    signify test data-set.}
  \item{iters}{Number of iterations. Each iteration consist of splitting the 
   data into train and test sets, performing the classification and storing 
   results}
  \item{SplitRatio}{ Splitting ratio: 
    \itemize{
      \item if \code{(0<=SplitRatio<1)} then \code{SplitRatio} fraction of 
      points from Y will be set toTRUE
      \item if \code{(SplitRatio==1)} then one random point from Y will be set 
      to TRUE
      \item if \code{(SplitRatio>1)} then \code{SplitRatio} number of points 
      from Y will be set to TRUE
    }
  }
  \item{RemCorrCol}{See \code{\link{msc.classifier.run}}.}
  \item{KeepCol}{See \code{\link{msc.classifier.run}}.}
  \item{ScaleType}{See \code{\link{msc.classifier.run}}.}
  \item{prior}{See \code{\link{msc.classifier.run}}.} 
  \item{same.sample}{See \code{\link{msc.classifier.run}}.}
  \item{method}{See \code{\link{msc.classifier.run}}.}
  \item{verbose}{boolean flag turns debugging printouts on.}
  \item{\ldots}{Additional parameters to be passed to classifiers. See 
    \code{method} for suggestions. } 
}

\details{
 This function follows standard cross-validation steps:
 \itemize{
   \item Class labels \code{Y} are used to divide data  \code{X} into train 
    set (with 
    known labels) and test set (labels are unknown and will be calculated)
   \item For number of iterations repeat the following steps:
    \itemize{ 
      \item split train data into temporary train and test sets using 
       \code{msc.sample.split} function
      \item train and test the chosen classifier using temporary train and test 
      data sets and \code{msc.classifier.run} function
    }
  \item Calculate the overall performance of the calassifer
  \item Train the classifier using the whole train data set (all labaled samples)
  \item Use this classifier to predict values of the whole test data set (all
   samples without label - NA.)
 }
}

\value{ 
  \item{Y}{Predicted class labels. If there were any unknown samples in input 
    data, marked by \code{NA}'s in input \code{Y}, than output \code{Y} will
    only hold prediction of those samples, otherwise prediction will be made 
    for all samples. }
  \item{Res}{Holds fraction of correct prediction during cross-validation for
    each iteration. \code{mean(Res)} will give you average acuracy. }
  \item{Tabl}{Contingency table of predictions shows all the input label 
    compared to output labels}
} 

\author{Jarek Tuszynski (SAIC) \email{jaroslaw.w.tuszynski@saic.com}} 
\note{This function is not fully tested and might be changed in future versions}

\seealso{
  \itemize{
  \item Input comes most likely from \code{\link{msc.preprocess.run}} 
  and/or \code{\link{msc.project.run}} functions.
  \item Uses \code{\link{msc.classifier.run}}, \code{\link{msc.features.select}}
   and \code{\link{msc.features.scale}} functions.
  \item Best classifier parameter set can be found by 
  \code{\link[e1071]{tune}} function from \pkg{e1071} package.
  \item Uses variety of classification algorithms: \code{\link[e1071]{svm}},
   \code{\link[nnet]{nnet}}, \code{\link{LogitBoost}}, \code{\link[MASS]{lda}},
   \code{\link[MASS]{qda}}, \code{\link[rpart]{rpart}}
  }
}

\examples{
  data(iris)
  A = msc.classifier.test(iris[,-5],iris[,5], method="LogitBoost", nIter=2) 
  A
  cat("correct classification in",100*mean(A$Res),"+-",100*sd(A$Res),"percent of cases\n")
}

\keyword{classif}