#===========================================================================#
# Written by Jarek Tuszynski. Copyright 2001-2003 SAIC.                     #
# Software developed in conjunction with the National Cancer Institute, and #
# distributed under "caBIO Software License" included in "COPYING" file.    #
#===========================================================================#

msc.features.select = function( x, y, RemCorrCol=0.98, KeepCol=0.6) 
{
  #=============================================
  # Remove some of the highly correlated columns
  #=============================================
  if (RemCorrCol || KeepCol) {
    Auc = colAUC(x, y)  # measure of ability to distinguish between classes based on the single feature
    if (min(dim(Auc))>1) Auc = apply(Auc, 2, max)
  }
  if (RemCorrCol) {
    idx = msc.features.remove(x, Auc, ccMin=RemCorrCol)
    Auc = Auc[idx]
  } else idx = 1:ncol(x)

  #=============================================
  # Remove columns with low AUC
  #=============================================
  if (KeepCol) {
    if (KeepCol<=1) 
      idx2 = which(Auc>KeepCol)      # keep columns with AUC > KeepCol
    else { 
      out  = sort(Auc, decreasing=TRUE, index.return=TRUE)
      idx2 = sort(out$ix[1:KeepCol]) # keep top "KeepCol" number of columns
    }  
    idx = idx[idx2]
  }
  return (idx)
}                                                                     
