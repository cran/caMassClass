#===========================================================================#
# Written by Jarek Tuszynski. Copyright 2001-2003 SAIC.                     #
# Software developed in conjunction with the National Cancer Institute, and #
# distributed under "caBIO Software License" included in "COPYING" file.    #
#===========================================================================#

msc.peaks.read.mzXML = function(filename)
{
  a = attr(filename,"class")
  flag = FALSE
  if (!is.null(a))
    if (a=="mzXML") flag=TRUE 
  mzXML = if(flag) filename else read.mzXML(filename)
  nScan = length(mzXML$scan)
  nPeak = 0
  for (i in 1:nScan)
    nPeak = nPeak + length(mzXML$scan[[i]]$mass)

  # save as list of peaks
  sNumb <- pNumb <- pIntn <- pMass <- numeric(nPeak) # initialize arrays for storing: peak number, peak intensity (height) and peak position (mass)
  j = 1
  for (i in 1:nScan) {
    out   = mzXML$scan[[i]]
    n     = length(out$peaks)     # read number of peaks
    if(n>0) {
      idx = j : (j+n-1)
      sNumb[idx] = rep(out$num,n) # create array of sample numbers coresponding to each peak
      pNumb[idx] = 1:n            # generate and store peak numbers
      pIntn[idx] = out$peaks      # store peak intensities (heights)
      pMass[idx] = out$mass       # store peak mass
      j = j+n
    }
  }                               # end of the loop over all the samples
  Data  = data.frame(sNumb, pNumb, pIntn, pMass) # concatinate all this data into a single table
  colnames(Data) = c("Spectrum.", "Peak.", "Intensity", "Substance.Mass") # column names
  return( list(scans=Data, mzXML=mzXML) )
}

#==============================================================================

msc.peaks.write.mzXML = function(scans, filename, mzXML=NULL, ...)
{
  if (is.null(mzXML)) mzXML = new.mzXML()
  
  # add caMassClass signature to mzXML data
  if (is.null(mzXML$dataProcessing) || regexpr("caMassClass", mzXML$dataProcessing)<0) {
    Version = packageDescription("caMassClass")$Version
    Time    = format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
    mzXML$dataProcessing = paste(mzXML$dataProcessing, 
      "    <dataProcessing>\n      <software type='processing'",
      " name='cran.r-project.org/caMassClass' version='", Version,
      "' completionTime='",Time,"'/>\n    </dataProcessing>\n", sep="")
  }

  # Convert list of peaks 'X' to mzXML$scan format
  CNames = c("Spectrum.Tag", "Spectrum.", "Intensity", "Substance.Mass") # column names
  if (!all(CNames %in% colnames(scans))) 
    stop("Input variable 'X' is not in proper format") 
  S   = scans$Spectrum.
  idx = which(!duplicated(S)) # remove duplicated sample numbers and tags
  nIdx = length(idx)
  mzXML$scan  = vector(mode="list", length=nIdx)
  for(i in 1:nIdx) {
    left  = idx[i]
    right = ifelse (i<nIdx, idx[i+1]-1, length(S))
    mzXML$scan[[i]] = list(mass =scans$Substance.Mass[left:right], 
                           peaks=scans$Intensity[left:right], 
                           num  =scans$Spectrum.[left],
                           parentNum=scans$Spectrum.[left], 
                           msLevel=1, maldi=NULL, scanOrigin=NULL, 
                           nameValue=NULL, precursorMz=NULL) 
    }
  # save the file
  write.mzXML(mzXML, filename,  ...)
}

