#===========================================================================#
# Written by Jarek Tuszynski. Copyright 2001-2003 SAIC.                     #
# Software developed in conjunction with the National Cancer Institute, and #
# distributed under "caBIO Software License" included in "COPYING" file.    #
#===========================================================================#

msc.mass.adjust.calc = function(X, scalePar=2, shiftPar=0.0005, AvrSamp=0)
{
  # [ShiftX, ScaleY, ShiftY, AvrSamp] = msc.mass.adjust(data, scalePar=2 shiftPar=0.0005, AvrSamp=0)
  # Scaling and mass adjustment function. It can take one ( when X is a matrix [nSamp x nFeat]) 
  # or multiple copies of data (X is an array [nSamp x nFeat x nCopy]).
  #
  # Input Parameters :
  #  scalePar - is used to control scaling: 0 means no scaling; 1 means that only mean will be matched, 
  #             2 means that both mean and medium will be matched (default)
  #  shiftPar - is used to control mass adjustment shifting sample has to improve correlation by 
  #             at least that amount to be considered. Default = 0.0005.
  #  AvrSamp  - is used if one is trying to normalize test set the same way train set was normalized.
  #             Than test set is using AvrSamp that was outputed from train set mass adjustment
  #
  # Output Parameters:
  #  ShiftX  - matrix [nSamp x nCopy] - integer number of positions sample should be shifted to the right (+) or left (-) 
  #  ScaleY  - matrix [nSamp x nCopy] - multiply each sample in order to normalize it
  #  ShiftY  - matrix [nSamp x nCopy] - substract this number from scaled sample (if matching mediums)
  #  AvrSamp - Use AvrSamp outputed from train set mass adjustment to process test set
 
  align = function(x, y, idx, scalePar, shiftPar)
  # Helper function to align and normalize 2 spectra "x" and "y".
  # full spectra are used for normalization, but only points listed in "idx" 
  # are used for mass drift adjustment (they should be as continous as possible).
  # See msc.mass.adjust help for scalePar and shiftPar parameters.
  {
    if (length(x)!=length(y)) 
      stop("msc.mass.adjust.calc: x and y have to have the same length")
    if (scalePar>0) {  
      if (scalePar==1) {  # scalePar==1 match only mean of x & y arrays
        medy = 0
        medx = 0
      } else {            # scalePar==2 match mean and medium of x & y arrays
        medy = median(y)
        medx = median(x)
      } 
      avrx = mean(x)-medx
      avry = mean(y)-medy
      scaleY = avrx / avry     
      shiftY = scaleY*medy-medx
      y = scaleY*(y-shiftY)     # scale y to match mean and medium of x
    } else {
      scaleY = 1
      shiftY = 0
    }
    x = x[idx]
    y = y[idx]
    n = length(idx)
    dim(x) = c(n,1)
    dim(y) = c(1,n)
    shr = (1:n-1) # shift right
    shl = (2:n)   # shift left
    x3 = t(rbind( c(0,x[shr]), t(x), c(x[shl],0) ))
    Cor = y %*% x3
    i  =  which.max(Cor)
    CC = Cor[2]
    cc = Cor[i]
    shiftX=0
    if (i==3) {
      y  = c(0, y[shr]) 
      while (cc-CC>shiftPar) {
        y  = c(0, y[shr]) 
        CC = cc
        cc = drop(y %*% x)
        shiftX = shiftX+1
      }
    } else if (i==1) {
      y  = c(y[shl], 0) 
      while (cc-CC>shiftPar) {
        y  = c(y[shl], 0)
        CC = cc
        cc = drop(y %*% x)
        shiftX = shiftX-1
      }
    }
    return( c(shiftX, scaleY, shiftY) )
  }

  # ========================
  # memory allocation
  # ========================
  nFeat = nrow(X)
  nSamp = ncol(X)
  nCopy = 1
  if (length(dim(X))==3) nCopy = dim(X)[3]
  ShiftX = matrix(0, nSamp, nCopy)
  ShiftY = matrix(0, nSamp, nCopy)
  ScaleY = matrix(1, nSamp, nCopy)
  
  # =================================================
  # align multiple copies (replicates) to each sample
  # =================================================
  if (nCopy>1) {              
    # average over samples and copies to calculate range of numbers
    # that will be used for mass adjustment (points above average)
    avr    = apply(X,1,mean)
    Range  = which(avr>mean(avr))
    # allign all copies of the sample to each other
    for (iSamp in 1:nSamp) {
      avr = X[,iSamp,1]
      for (iCopy in 2:nCopy) {
        out = align(avr, X[,iSamp,iCopy], Range, scalePar, shiftPar)
        a = out[1]
        ShiftX[iSamp,iCopy] = a
        ScaleY[iSamp,iCopy] = out[2]
        ShiftY[iSamp,iCopy] = out[3]
        Y = out[2]*X[max(1,1-a):min(nFeat,nFeat-a), iSamp, iCopy] - out[3]
        dim(Y) = c(1,nFeat-abs(a))
        if (a==0) XX = Y else
        if (a> 0) XX = cbind(matrix(0,1, a), Y) else
        if (a< 0) XX = cbind(Y, matrix(0,1,-a))
        avr = ((iCopy-1)*avr + XX)/iCopy
        X[,iSamp,iCopy] = XX
        #print(sprintf('%i %1.0f %f %f\n', iSamp,a,out[2],out[3]))
      }
    }
  }
  # =============================================
  # align each sample to the mean of all samples
  # =============================================
  if (length(dim(X))==3) {
    nCopy = dim(X)[3] 
    X=aperm(X,c(1,3,2))
    dim(X) = c(nFeat*nCopy, nSamp) # same as cbind of all copies
  }
  nR = nrow(X)
  nIter = 4
  if (length(AvrSamp)>1) nIter=1  # matching to existing signature - 1 iter is enough
  for (iter in 1:nIter) {         # try it 4 times since mean will be changing
    if (nIter==1) avr = AvrSamp   # Test set processing: use Avr from Train Set
    else avr = rowMeans(X)        # calculate average signal
    Range = which(avr>mean(avr))  # match peaks not valleys
    change = 0
    for (iSamp in 1:nSamp) {
      out = align(avr, X[,iSamp], Range, scalePar, shiftPar)
      a = out[1]                  # shift whole signature to the right or left
      change = change+abs(a)      # check if changing
      ShiftX[iSamp,] = ShiftX[iSamp,]+a
      ScaleY[iSamp,] = ScaleY[iSamp,]*out[2]
      ShiftY[iSamp,] = ShiftY[iSamp,]*out[2]+out[3]
      Y = X[max(1,1-a):min(nR,nR-a),iSamp]  # mass adjustment
      Y = out[2]*Y - out[3]                 # normalization
      if (a==0) X[,iSamp] = Y else
      if (a> 0) X[,iSamp] = c(rep(0,a), Y) else
      if (a< 0) X[,iSamp] = c(Y, rep(0,-a))
      #print(sprintf('%i %1.0f %f %f\n', iSamp,a,out[2],out[3]))
    }
    if (!change) break # stopped changing -> we are done
    scalePar = 0       # scaling need to be done only once at first iteration
  }
  return( list(ShiftX=ShiftX, ScaleY=ScaleY, ShiftY=ShiftY, AvrSamp=rowMeans(X)) )
}

#========================================================================================  

msc.mass.adjust.apply = function(X, shiftX, scaleY, shiftY)
  # X = ShiftScale(X, ShiftX, ScaleY, ShiftY)
  # Apply shifting and scaling calculated by msc.mass.adjust to the data "X".
  # It can take one ( when X is a matrix [nSamp x nFeat]) 
  # or multiple copies of data (X is an array [nSamp x nFeat x nCopy]).
  #
  # Input Parameters :
  #  ShiftX  - matrix [nSamp x nCopy] - integer number of positions sample will be shifted to the right (+) or left (-) 
  #  ScaleY  - matrix [nSamp x nCopy] - multiply each sample in order to normalize it
  #  ShiftY  - matrix [nSamp x nCopy] - substract this number from scaled sample (if matching mediums)
{
  d      = dim(X)          # dimentions of the data
  dNames = dimnames(X)
  nR     = d[1]
  nSamp  = prod(d)/nR
  dim(X) = c(nR, nSamp)    # and resize the data
  for (iSamp in 1:nSamp) {
    a = shiftX[iSamp]              
    Y = X[max(1,1-a):min(nR,nR-a),iSamp] # mass adjustment
    Y = scaleY[iSamp]*Y-shiftY[iSamp]    # normalization
    if (a==0) X[,iSamp] = Y else
    if (a> 0) X[,iSamp] = c(rep(0,a), Y) else
    if (a< 0) X[,iSamp] = c(Y, rep(0,-a))
  }
  dim(X) = d               # dimentions of the data
  dimnames(X) = dNames
  return(X)
}
 
#========================================================================================  

msc.mass.adjust = function(X, scalePar=2, shiftPar=0.0005, AvrSamp=0)
{
  out = msc.mass.adjust.calc (X, scalePar=scalePar, shiftPar=shiftPar, AvrSamp=AvrSamp)
  X   = msc.mass.adjust.apply(X, out$ShiftX, out$ScaleY, out$ShiftY)
  return(X)
}
 
