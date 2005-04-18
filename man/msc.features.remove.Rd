\name{msc.features.remove}
\alias{msc.features.remove}
\title{Remove Highly Correlated Features}
\description{Remove Highly Correlated Features. The function checks neighbor
  features looking for highly correlated ones and removes one of them. Used in 
  order to drop dimensionality of the data.}

\usage{msc.features.remove(Data, Auc, ccMin=0.9, verbose=FALSE)}

\arguments{
  \item{Data}{Data containing one sample per row and one feature per column.} 
  \item{Auc}{A measure of usefulness of each column/feature, used to choose
    which one of two highly correlated columns to remove. Usually a measure of
    discrimination power of each feature as measured by 
    \code{\link{colAUC}}, student t-test or other method. See details.}
  \item{ccMin}{Minimum correlation coefficient of "highly correlated" columns.} 
  \item{verbose}{Boolean flag turns debugging printouts on.} 
}

\details{
 If \code{\link{colAUC}} was used and there were more than two classes present
 than \code{Auc} is a matrix with multiple measurments for each feature. In 
 such a case \code{Auc = apply(Auc, 2, mean)} is run in order to extract a single 
 measure per feature. If other measures are desired, like  
 \code{Auc = apply(Auc, 2, max)}, than they should be called beforehand.
}

\value{Vector of column indexes to be kept.} 
\author{Jarek Tuszynski (SAIC) \email{jaroslaw.w.tuszynski@saic.com}} 
\seealso{ 
  \itemize{
  \item Used by \code{\link{msc.classifier.test}} and \code{\link{msc.features.select}} 
  functions.
  \item Uses \code{\link{colAUC}} 
  }
}

\examples{
  # load input data 
  if (!file.exists("Data_IMAC.Rdata")) example("msc.project.read")
  load("Data_IMAC.Rdata")
  
  X = t(X[,,1])
  auc = colAUC(X,SampleLabels)
  quantile(auc)
  cidx = msc.features.remove(X, auc, verbose=TRUE)
  Y = X[,cidx]
}

\keyword{ts}
\keyword{classif}
