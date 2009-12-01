#===========================================================================#
# Written by Jarek Tuszynski. Copyright 2001-2003 SAIC.                     #
# Software developed in conjunction with the National Cancer Institute      #
# Distributed under GNU General Public License version 3                    #
#===========================================================================#

msc.peaks.write.csv = function(X, fname)
{ # save the table into csv file (this is the same format as returned from ciphergen software)
  CNames = c("Spectrum.Tag", "Spectrum.", "Intensity", "Substance.Mass") # column names
  if (!all(CNames %in% colnames(X))) 
    stop("Input variable 'X' is not in proper format") 
  #write.csv(X, file=fname, row.names=FALSE)  
  write.table(X, file=fname, sep=",", row.names=FALSE)  
}

msc.peaks.read.csv = function(fname, mzXML.record=FALSE)
{
  X = read.csv(fname, as.is=TRUE) # read input file with peaks
  CNames = c("Spectrum.Tag", "Spectrum.", "Intensity", "Substance.Mass") # column names
  if (!all(CNames %in% colnames(X))) 
    stop("Input variable 'X' is not in proper format") 
  if (mzXML.record) { # create mzXML record 
    S =  unique(X$Spectrum.)
    sha1  = digest(fname, alg="sha1", file=TRUE)
    mzXML$parentFile = paste(mzXML$parentFile, 
      "    <parentFile filename='file:///", fname, 
      "' fileType='processedData' fileSha1='",sha1,"'/>\n", sep="")
    for (i in 1:length(S)) {
      j = S[i]
      scanOrigin = paste("<scanOrigin parentFileID='",sha1,"' num='",j,"'/>", sep="");
      mzXML$scan[[i]] = list(mass=NULL, peaks=NULL, num=j, parentNum=j, 
        msLevel=1, header=NULL, maldi=NULL, scanOrigin=scanOrigin,
        precursorMz=NULL, nameValue=NULL)
    }
    attr(X,"mzXML") = mzXML
  }
  
  return(X)
}

