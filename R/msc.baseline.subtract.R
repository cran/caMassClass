#===========================================================================#
# Written by Jarek Tuszynski. Copyright 2001-2003 SAIC.                     #
# Software developed in conjunction with the National Cancer Institute      #
# Distributed under GNU General Public License version 3                    #
#===========================================================================#

msc.baseline.subtract = function(X, ...) 
{
  library(PROcess)
  d      = dim(X)                 # dimentions of the data
  dNames = dimnames(X)            # dimention names of the data
  dim(X) = c(d[1], prod(d)/d[1])  # cast into a matrix form
  nFeat  = nrow(X)                # number of rows/features
  nSamp  = ncol(X)                # number of columnss/samples
  f      = matrix(0, nFeat, 2)    # data format needed by bslnoff
  f[,1]  = as.numeric(dNames[[1]])# extract masses (x-axis data) stored as row names         
  for (iSamp in 1:nSamp) {        # for each sample call bslnoff
    f[,2] = X[,iSamp]
    g = bslnoff(f, method="loess", ...)
    X[,iSamp] = g[,2]
  }
  dim(X) = d                      # restore original dimentionality of the data
  dimnames(X) = dNames            # restore original names
  return(X)
}

