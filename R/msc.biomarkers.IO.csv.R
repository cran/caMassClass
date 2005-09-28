#===========================================================================#
# Written by Jarek Tuszynski. Copyright 2001-2003 SAIC.                     #
# Software developed in conjunction with the National Cancer Institute, and #
# distributed under "caBIO Software License" included in "COPYING" file.    #
#===========================================================================#

msc.biomarkers.write.csv = function(X, fname)
{
  d = dim(X)                       # dimentions of the data
  if (length(d)==3) {              # if this is 3D data than convert it to 2D
	  rNames  = rownames(X)          # extract mass from input data 
	  Samples = colnames(X)          # extract sample names
     cNames=NULL
    for (iCopy in 1:d[3])          # create new row names
      cNames = c(cNames, paste(Samples,"(", as.character(iCopy), ")", sep=""))
    dim(X) = c(d[1], prod(d)/d[1]) # and resize the data
    dimnames(x) = list(rNames, cNames)
  } 
  A = data.frame(t(X))
  A[is.na(A)] = 0                  # convert NA to zeros
  write.csv(A, file=fname, col.names=NA)
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

