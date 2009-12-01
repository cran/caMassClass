#===========================================================================#
# Written by Jarek Tuszynski. Copyright 2001-2003 SAIC.                     #
# Software developed in conjunction with the National Cancer Institute      #
# Distributed under GNU General Public License version 3                    #
#===========================================================================#

msc.features.remove = function(Data, Auc, ccMin=0.9, verbose=FALSE)
{
# remove worse one of any two higly correlated columns
# Inputs:
#   Data  - 
#   auc   - a measure of how good is each column (usually output of colAUC function)
#   ccMin - minimum correlation coefficient of "higly correlated" columns
  if (is.matrix(Auc) & min(dim(Auc))>1 ) Auc = apply(Auc, 2, mean)
  nC    = ncol(Data)
  mask  = logical(nC)                # global mask
  mask  = !mask                      # set all elements in the mast to TRUE
  ccoef = numeric(nC)                # correlation coefficient between column and column+1
  ccoef[nC] = -1  
  clist = (1:nC)                     # list of remaining columns 
  idx = which(ccoef==0)              # idx is a list of columns to be looked at
  iter = 1
  if (verbose) {
    print("iteration  nCol  nLeads")
    print(sprintf("%9i %5i %7i" ,as.integer(iter), nC, length(idx) ))
  }
  while (length(idx)>0) {            # loop until done
    iter = iter+1
    for (i in 1:length(idx)) {
      k = idx[i]                     # column number
      cc = cor(Data[,k], Data[,k+1]) # find correlation coeficient between column k and k1
      if (is.na(cc)) cc=1
      ccoef[k] = cc                  # ... and store it so we do not have to calculate it again
      if (cc>ccMin) {                # if bigger than ccMin than columns are "higly correlated"
        j = k
        if (Auc[k+1]<Auc[k]) j=k+1   # column to be removed will have lower auc
        mask [j] = FALSE             # mark to be removed
        ccoef[j] = 0                 # since a column changed 2 ccoefs will have to be recalculated
        if (j>1) ccoef[j-1]=0 
      }   
    }
    Data  = Data[,mask]
    ccoef = ccoef[mask]
    Auc   = Auc  [mask]
    clist = clist[mask]
    mask  = mask [mask]
    nC  = length(ccoef)
    ccoef[nC] = -1
    idx = which(ccoef==0)            # idx is a list of columns to be looked at
    if (verbose)
     print(sprintf("%9i %5i %7i" ,as.integer(iter), nC, length(idx) ))
  }
  return(clist)
}
