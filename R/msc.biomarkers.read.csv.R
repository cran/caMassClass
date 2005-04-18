#===========================================================================#
# Written by Jarek Tuszynski. Copyright 2001-2003 SAIC.                     #
# Software developed in conjunction with the National Cancer Institute, and #
# distributed under "caBIO Software License" included in "COPYING" file.    #
#===========================================================================#

msc.biomarkers.write.csv = function(fname, X)
{
  A = data.frame(t(X))
  A[is.na(A)] = 0
  write.table(A, file=fname, sep=",", col.names=NA )
}

msc.biomarkers.read.csv = function(fname) 
{
  X = read.csv(fname, comment.char = "", header=TRUE)  
  SampleName = X[,1]
  X = t(X[,-1])
  X[X==0] = NA
  colnames(X) = as.character(SampleName)
  return (X)
}

