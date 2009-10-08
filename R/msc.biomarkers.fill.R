#===========================================================================#
# Written by Jarek Tuszynski. Copyright 2001-2003 SAIC.                     #
# Software developed in conjunction with the National Cancer Institute, and #
# distributed under "caBIO Software License" included in "COPYING" file.    #
#===========================================================================#

msc.biomarkers.fill = function( X, Bmrks, BinBounds, FillType=0.9)
# X - a matrix holding the baseline-substracted spectra, where each spectrum is in a separate column. The mass m/z values are stored as row names 
#     and column-names are stored as the spectrum names. 
# Bmrks - a feature matrix created out of peaks extracted from X
# BinBounds - column numbers corresponding to left and right border of each biomarker range. 
# FillType  - what to do with spaces in output matrix without any peaks? 
#            if 0<=FillType<=1 - fill spaces with  quantile(FillType). For example if FillType=1/2 than medium will be used, 
#                                if FillType = 1 than maximum value will be used, is FillType=0.9 than maximum will be used 
#                                after discarding 10% of "outliers", if FillType = 0 than minimum value will be used (not recomended)
#            if FillType<0 than empty spaces will not be filled
#            if FillType>1 than X value closest to the center of the bin will be used 
{
  d      = dim(X)
  dNames = dimnames(X)
  nFeat  = nrow(Bmrks)

  #=======================================
  # filling empty positions in feat matrix
  #=======================================
  if (is.array(X) & FillType>=0) {
    dim(X) = c(d[1], prod(d)/d[1])                                 # if X is a 3D data cube than make it into a matrix
    if (ncol(Bmrks)!=ncol(X)) 
      stop("msc.biomarkers.fill: Unequal number of columns/samples in X(",ncol(Bmrks),") and Bmrks arrays(",ncol(X),")")
    mass   = as.numeric(dNames[[1]])                               # masses in X

    for (i in 1:nFeat) {                                           # for each biomarker bin....
      cidx = which(is.na(Bmrks[i,]))                               # where are the holes to fill?
      n = length(cidx)
      if(n==0) next
      if (FillType>=0 & FillType<=1) {                             # find quantile of the bin
        ridx = which(mass>=BinBounds[i,1] & mass<=BinBounds[i,2])  # find all points of this bin
        if(length(ridx)>1) {
          if(n==1) Bmrks[i,cidx[1]] = quantile( X[ridx, cidx[1]], probs=FillType)
          else Bmrks[i,cidx] = apply( X[ridx, cidx], 2, quantile, probs=FillType)
        } else {
          if(n==1) Bmrks[i,cidx[1]] = X[ridx[1], cidx[1]]
          else Bmrks[i,cidx] = X[ridx[1], cidx]
        }
      } else if (FillType==2) {
         dis  = abs(mass - (BinBounds[i,1]+BinBounds[i,2])/2)
        ridx = which.min(dis)                                      # find nearest point to bin center
        Bmrks[i,cidx] = X[ridx,cidx]
      } else if (FillType==3) {
         A[is.na(A)] = 0
      }
     }
    dim(X) = d
  }
    
  #=======================================
  # cast Bmrks in to the same format as X is
  #=======================================
  if (length(d)==3 && d[3]>1) {
    mass =  rownames(Bmrks)
    dim(Bmrks) = c(nrow(Bmrks), d[2], d[3])   # convert 2D matrix back to a 3D data cube
    dimnames(Bmrks) = list( mass, dNames[[2]], dNames[[3]] )
  } else colnames(Bmrks) = dNames[[2]]
  return(Bmrks)
}

