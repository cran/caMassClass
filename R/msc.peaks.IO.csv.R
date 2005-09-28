#===========================================================================#
# Written by Jarek Tuszynski. Copyright 2001-2003 SAIC.                     #
# Software developed in conjunction with the National Cancer Institute, and #
# distributed under "caBIO Software License" included in "COPYING" file.    #
#===========================================================================#

msc.peaks.write.csv = function(X, fname)
{ # save the table into csv file (this is the same format as returned from ciphergen software)
  CNames = c("Spectrum.Tag", "Spectrum.", "Intensity", "Substance.Mass") # column names
  if (!all(CNames %in% colnames(X))) 
    stop("Input variable 'X' is not in proper format") 
	#write.csv(X, file=fname, row.names=FALSE)	
	write.table(X, file=fname, sep=",", row.names=FALSE)	
}

msc.peaks.read.csv = function(fname)
{
  X = read.csv(fname, as.is=TRUE) # read input file with peaks
  CNames = c("Spectrum.Tag", "Spectrum.", "Intensity", "Substance.Mass") # column names
  if (!all(CNames %in% colnames(X))) 
    stop("Input variable 'X' is not in proper format") 
  return(X)
}

