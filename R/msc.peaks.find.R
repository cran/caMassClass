#===========================================================================#
# Written by Jarek Tuszynski. Copyright 2001-2003 SAIC.                     #
# Software developed in conjunction with the National Cancer Institute, and #
# distributed under "caBIO Software License" included in "COPYING" file.    #
#===========================================================================#

msc.peaks.find=function(X, SNR=2, span=c(81,11), zerothresh=0.9) 
# Input Parameters:
#  X:  the data in one ( when X is a matrix [nFeat x nSamp]) or multiple copies  (X is an array [nFeat x nSamp x nCopy])
{
  if (abs(zerothresh)>1) 
    stop("Zerothresh (",zerothresh,") smaller than minus one or larger than one is incorrect - using default ", zerothresh=0.9);
  GlobalThresh = (zerothresh>=0 && zerothresh<1)
  if (GlobalThresh) zerothresh=quantile(X,zerothresh) 
  span = span+1-span%%2              # make span numbers odd by adding one if they are even
  if (span[1]<span[2]) {tmp=span[1]; span[1]=span[2]; span[2]=tmp; } # swap if needed
	mass  = as.numeric(rownames(X))    # extract mass from input data 
	SampleNames = colnames(X)          # extract sample names
  d     = dim(X)                     # dimentions of the data
  if (length(d)==3) {                # if this is 3D data than convert it to 2D
    nCopy = d[3]
    cNames=NULL
    for (iCopy in 1:nCopy)           # create new row names
      cNames = c(cNames, paste(SampleNames,"(", as.character(iCopy), ")", sep=""))
    dim(X) = c(d[1], prod(d)/d[1])   # and resize the data
    SampleNames = cNames
  } 
  nSamp = ncol(X)                    # number of columns/samples
	cnts  = numeric(nSamp)             # initialize array for storing # peaks per sample
  pNumb=NULL; pIntn=NULL; pMass=NULL;# initialize arrays for storing: peak number, peak intensity (height) and peak position (mass)
	for (j in 1:nSamp) {               # loop over all the samples
		x      = X[,j]                   # put data of next sample into agreed format  
    thresh = if(GlobalThresh) zerothresh else quantile(x,-zerothresh)   # choose between absolute threshold and one calculated from this sample
    sig    = runmean(x, span[2])     # smoothed signal using rectangular kernel                          
    rMax   = runmax (x, span[2])     # calculate running max  - finds local maxima                           
    rAvr   = runmed (x, span[1])     # calculate running median - needed for MAD calculations
    rStd   = runmad (x, span[1], center=rAvr) # calculate running mad  - local median average diviation - robust diviation measure scaled to mimic sd
    peak   = (rMax == x) & (sig > thresh) & (sig-rAvr > SNR*rStd) # the final peak selection
    idx    = which(peak)             # convert from boolean mask to index
		cnts[j]= length(idx)             # count peaks and store their number
		pNumb  = c(pNumb, 1:cnts[j])     # generate and store peak numbers
		pIntn  = c(pIntn, rAvr[idx])     # store peak intensities (heights)
		pMass  = c(pMass, mass[idx])     # store peak mass
 	}                                  # end of the loop over all the samples
	sName = rep(SampleNames, cnts)     # create array of sample names coresponding to each peak
	sNumb = rep(1:nSamp    , cnts)     # create array of sample numbers coresponding to each peak
  Data  = data.frame(sName, sNumb, pNumb, pIntn, pMass) # concatinate all this data into a single table
  colnames(Data) = c("Spectrum.Tag", "Spectrum.", "Peak.", "Intensity", "Substance.Mass") # column names
  return (Data)
}

