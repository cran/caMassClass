#===========================================================================#
# Written by Jarek Tuszynski. Copyright 2001-2003 SAIC.                     #
# Software developed in conjunction with the National Cancer Institute, and #
# distributed under "caBIO Software License" included in "COPYING" file.    #
#===========================================================================#

msc.peaks.alignment = function(S, M, H, Tag=0, SampFrac=0.3, 
                               BinSize=c(0.002, 0.008), ...)
# given :
#  S - peak sample number
#  M - peak center mass
#  H - peak height
#  BinSize - upper and lower bound of bin size measured as (R-L)/mean(R,L)
#            Based on expected experimental variation in the m/z values. size 
#            of any bin is measured as (R-L)/mean(R,L) where L and R are masses 
#            (m/z values) of left and right boundaries if length of BinNumber 
#            is 2 - resulting bin sizes will all be between  BinSize[1] and 
#            BinSize[2] since SELDI data is often asumed to have 3% mass drift 
#            than a good bin size is twice that number (0.006) of each peak 
#            align peaks heights from different samples together into bins 
#            (features) stored in 2D array feat (sample by feature).
{
  #=======================================
  # extract sample names from peak info file
  #=======================================
  nPeak = length(M)                  # how many peaks?
  if (nPeak!=length(S) | nPeak!=length(H)) 
    stop("msc.peaks.alignment: Unequal size arrays in msc.peaks.align")
  if (length(Tag)==nPeak) {
    mask = !duplicated(S)            # remove duplicated sample numbers and tags
    Tg   = Tag[mask]                 # get tags related to unique sample numbers
    x    = sort(S[mask], index=TRUE) # sort sample numbers (they are most ...
    sNames = Tg[x[[2]]]              # ... likely already sorted)
  } else sNames=0
  
  #=========================
  # Prepare for calculations
  #=========================
  x   = sort(M, index=TRUE)          # sort peaks by mass
  M   = x[[1]]                       # sorted mass
  idx = x[[2]]                       # peak order after sorting
  S   = S[idx]                       # arrange sample numbers in the same order
  H   = H[idx]                       # arrange peak heights in the same order
  nPeak = length(M)                  # how many peaks?
  M1 = c(M[2:nPeak], M[nPeak])       # shift M one to the right
  dM = 2*(M1-M)/(M+M1)               # find relative spacing between peaks
  
  #====================================================
  # Perform peak alignment by Divide & conquer approach 
  #====================================================
  bin   = msc.peaks.clust(dM, S, BinSize=BinSize, ...) # bin is boolean array marking left boundary of each bin
  Left  = which(bin==1)              # left boundary of each bin
  bin   = cumsum(bin)                # find bin numbers for each peak in S array
  nSamp = max(S)                     # number of samples
  nFeat = length(Left)               # number of bins
  nul   = -10000
  Bmrks = matrix(nul,nFeat,nSamp)    # init feature (biomarker) matrix
  for (j in 1:nPeak)                 # Bmrks will store height H of each peak
    Bmrks[bin[j], S[j]] = max(Bmrks[bin[j], S[j]], H[j]) # 2 peaks from the ...
  Bmrks[Bmrks==nul] = NA   #... same spectrum in the same bin should be a rare 

  #========================================
  # find boundaries of each bin / biomarker
  #========================================
  BinBounds = c(Left, Left[2:nFeat]-1, nPeak )  # index of left -most and Right-most peak in the cluster 
  BinBounds = M[BinBounds]           # mass of left -most and Right-most peak in the cluster 
  dim(BinBounds) = c(nFeat,2)
  colnames(BinBounds) = c("Left","Right")
  bin = (BinBounds[,1]+BinBounds[,2])/2 # find center of each bin
  rownames(Bmrks) = paste("M",signif(bin,6),sep="")
  if (length(sNames)==ncol(Bmrks)) colnames(Bmrks) = sNames

  #=======================================
  # remove feature/columns with not enough peaks
  #=======================================
  if (SampFrac>0) {
    numCol = apply(!is.na(Bmrks), 1, sum) # count nonempty spaces in each column
    keep   = (numCol > SampFrac*ncol(Bmrks)) 
    Bmrks  = Bmrks    [keep,]     
    BinBounds = BinBounds[keep,]
  }

  return( list(Bmrks=Bmrks, BinBounds=BinBounds) )  
}

#===========================================================================================

msc.peaks.align = function(Peaks, SampFrac=0.3, BinSize=c(0.002, 0.008), ...)
# Peaks - could have 2 formats (a filename where to find the data or data itself:   
#       - a data frame storing peak information
#       - a name of `.csv' file in the same format as Ciphergen's peakinfo file 
#         has to contain columns : Spectrum., Intensity and Substance.Mass.
#  BinSize - upper and lower bound of bin size measured as (R-L)/mean(R,L)
#            Based on expected experimental variation in the m/z values. size 
#            of any bin is measured as (R-L)/mean(R,L) where L and R are masses 
#            (m/z values) of left and right boundaries if length of BinNumber 
#            is 2 - resulting bin sizes will all be between  BinSize[1] and 
#            BinSize[2] since SELDI data is often asumed to have 3% mass drift 
#            than a good bin size is twice that number (0.006) of each peak 
#            align peaks heights from different samples together into bins 
#            (features) stored in 2D array feat (sample by feature).
{
  if (is.character(Peaks)) 
    Peaks = msc.peaks.read.csv(Peaks) # read input file with peaks
    
  out = msc.peaks.alignment(Peaks$Spectrum., Peaks$Substance.Mass, Peaks$Intensity, 
         Tag=Peaks$Spectrum.Tag, BinSize=BinSize, SampFrac=SampFrac, ...)

  return(out)
}

