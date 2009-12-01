#===========================================================================#
# Written by Jarek Tuszynski. Copyright 2001-2003 SAIC.                     #
# Software developed in conjunction with the National Cancer Institute      #
# Distributed under GNU General Public License version 3                    #
#===========================================================================#

msc.features.scale = function( xtrain, xtest, 
      type=c("min-max", "avr-std", "med-mad")) 
{ 
  type = match.arg(type)  
  nFeat = ncol(xtrain)
  for( iFeat in 1:nFeat) {
    X = xtrain[,iFeat]
    mask = is.finite(X)
    x = X[mask]                  # remove na & inf
    if (type=="min-max") {       # linear scaling: min->0, max->1
      cDif   = min(x)
      cDiv   = max(x)-cDif
    } else if (type=="avr-std") {# linear scaling: mean->0, mean+sd->1
      cDif   = mean(x)
      cDiv   = sd(x)
    } else if (type=="med-mad") {# linear scaling: median->0, median+mad->1
      cDif   = median(x)
      cDiv   = mad(x, center=cDif)
    }
    xtrain[mask,iFeat] =  (x -cDif)/cDiv
    X = xtest[,iFeat]
    mask = is.finite(X)
    xtest [mask,iFeat] =  (X[mask]-cDif)/cDiv
  }
  return(list(xtrain=xtrain, xtest=xtest))
}
