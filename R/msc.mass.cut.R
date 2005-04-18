#===========================================================================#
# Written by Jarek Tuszynski. Copyright 2001-2003 SAIC.                     #
# Software developed in conjunction with the National Cancer Institute, and #
# distributed under "caBIO Software License" included in "COPYING" file.    #
#===========================================================================#

msc.mass.cut = function( X, MinMass=3000) 
# Remove features/columns with no data and Cut Low Masses
# Input / Output Parameters:
#  X:  the data in one ( when X is a matrix [nFeat x nSamp]) 
# or multiple copies  (X is an array [nFeat x nSamp x nCopy])
{  
  d     = dim(X)             # dimentions of the data
  nFeat = d[1]
  Mass  = as.numeric(rownames(X))
  mask  = logical(nFeat)     # rows to delete
  if (length(d)==1) {        # 1D case
    X    = X[Mass>MinMass]
  } else if (length(d)==2) { # 2D case
    for (iFeat in 1:nFeat) mask[iFeat] = !all(X[iFeat, ]==X[iFeat,1  ])
    mask = (mask & Mass>MinMass)
    X    = X[mask,]
  } else if (length(d)==3) { # 3D case
    for (iFeat in 1:nFeat) mask[iFeat] = !all(X[iFeat,,]==X[iFeat,1,1])
    mask = (mask & Mass>MinMass)
    X    = X[mask,,]
  }
  return( X )
}