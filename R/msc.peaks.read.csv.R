#===========================================================================#
# Written by Jarek Tuszynski. Copyright 2001-2003 SAIC.                     #
# Software developed in conjunction with the National Cancer Institute, and #
# distributed under "caBIO Software License" included in "COPYING" file.    #
#===========================================================================#

msc.peaks.write.csv = function(fname, X)
{ # save the table into csv file (this is the same format as returned from ciphergen software)
  colnames(X) = c("Spectrum.Tag", "Spectrum.", "Peak.", "Intensity", "Substance.Mass") # column names
	write.table(X, file=fname, sep=",", row.names=FALSE)	
}

msc.peaks.read.csv = function(fname)
{
   read.csv(fname, as.is=TRUE) # read input file with peaks
}

