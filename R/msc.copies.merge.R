#===========================================================================#
# Written by Jarek Tuszynski. Copyright 2001-2003 SAIC.                     #
# Software developed in conjunction with the National Cancer Institute, and #
# distributed under "caBIO Software License" included in "COPYING" file.    #
#===========================================================================#

msc.copies.merge = function( X, mergeType, PeaksOnly=TRUE)              
# mergeType : if copies are to be concatinated (as separate samples) set mergeType to 1 otherwise set it to 0
#             if copies are to be averaged add 2
#             if worst copy is to be deleted add 8 or 
#             if in case of a single bad copy of a sample one wants to replace it with the best copy add 4
{
  d = dim(X)
  if ( length(d)==3 ) { 
    nFeat = nrow(X)
    nSamp = ncol(X)
    nCopy = dim(X)[3]
    out = msc.sample.correlation(X, PeaksOnly=PeaksOnly)    
    outerCor = out$outerCor
    innerCor = out$innerCor
    a = mergeType %/% 4
    b = mergeType %%  4
  
    score = outerCor+innerCor
    # what to do with bad copies ?
    if (a==2) {               # drop worse copy
      worst = max.col(-score) # worse copy of each sample
      for (iSamp in 1:nSamp) X[,iSamp,worst[iSamp]] =  X[,iSamp,nCopy] # overwrite the worst sample by the last one
      nCopy = nCopy-1
      X = X[,,1:nCopy]        # drop last copy
    } else if (a==1) {        # if needed overwrite bad copies by better versions
      idx   = which( apply(innerCor,1,min) < apply(outerCor,1,max))    # samples that have problems
      best  = max.col( score) # best  copy
      worst = max.col(-score) # worse copy
      for (i in 1:length(idx)) {
        iSamp = idx[i]
        X[,iSamp,worst[iSamp]] = X[,iSamp,best[iSamp]]
      }
    }

    # how to merge copies
    if (b>0 && nCopy>1) { # concatinate all copies and/or mean of all copies
      SampleNames = colnames(X)
      XX=NULL
      if (b!=1) { # average copies
        XX = apply(X, c(1, 2), mean)
        colnames(XX) = paste("avr(",SampleNames,")", sep="")
      }
      if (b!=2) 
       for (iCopy in 1:nCopy) { # concatinate all copies 
         colnames(X) = paste(SampleNames,"(", as.character(iCopy), ")", sep="")
         XX = cbind(XX, X[,,iCopy])
       }
      X = XX
      rm(XX)
    }
  } 

  return (X)
}
