#===========================================================================#
# Written by Jarek Tuszynski. Copyright 2001-2003 SAIC.                     #
# Software developed in conjunction with the National Cancer Institute      #
# Distributed under GNU General Public License version 3                    #
#===========================================================================#

msc.sample.correlation = function(X, PeaksOnly=FALSE)
# function [innerCor, outerCor] = msc.sample.correlation(X)
# Input: X:  the data in one ( when X is a matrix [nFeat x nSamp]) or multiple copies  (X is an array [nFeat x nSamp x nCopy])
# Output: innerCor - average corelation between multiple copies of the same sample
#         outerCor - average correlation between each sample and every other sample within the same copy
#
# find correlation between rows in multiple copies of the data and average correlation between 
# different samples within the same data set
{
  d = dim(X)
  nFeat = d[1]
  nSamp = d[2]
  nCopy = if(length(d)==3) d[3] else 1
  dNames = dimnames(X)
  if (length(d)==2) X = as.matrix(X)
  dim(X) = c(nFeat, nSamp, nCopy)
  innerCor = matrix(1, nSamp, nCopy)
  outerCor = matrix(1, nSamp, nCopy)
  rownames(innerCor) = dNames[[2]]
  rownames(outerCor) = dNames[[2]]
  if (length(d)==3) {
    colnames(innerCor) = dNames[[3]]
    colnames(outerCor) = dNames[[3]]
  }
  
  # average over samples and copies to calculate range of numbers
  # that will be used for correlation
  if (PeaksOnly) {
    avr   = apply(X,1,mean)
    Range = which(avr>mean(avr))
    X = X[Range,,]
  }
  if (nCopy>1) 
   for (iSamp in 1:nSamp) {
     x = cor(X[,iSamp,])
     innerCor[iSamp, ] = t((colSums(x)-1) / (nCopy-1))
   }
  for (iCopy in 1:nCopy) {
    x = cor(X[,,iCopy])
    outerCor[,iCopy] = (colSums(x)-1) / (nSamp-1)
  }
  #if (nCopy==2) innerCor = innerCor[,1]
  return (list(innerCor=innerCor, outerCor=outerCor))
}

