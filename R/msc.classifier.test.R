#===========================================================================#
# Written by Jarek Tuszynski. Copyright 2001-2003 SAIC.                     #
# Software developed in conjunction with the National Cancer Institute      #
# Distributed under GNU General Public License version 3                    #
#===========================================================================#
               
msc.classifier.test = function( X, Y, iters=50, SplitRatio=2/3, verbose=FALSE,
                  RemCorrCol=0, KeepCol=0, prior=1, same.sample=NULL,
                  ScaleType=c("none", "min-max", "avr-std", "med-mad"),
                  method=c("svm", "nnet", "lda", "qda", "LogitBoost", "rpart"), ...) 
{ 
  if (length(Y)!=nrow(X)) 
    stop("msc.classifier.test: Number of Samples (nrow(X)) does not mach number of labels in Y")
  n = length(same.sample)
  if (n>0 & n!=nrow(X))
    stop("msc.classifier.test: Number of Samples (nrow(X)) does not mach number of elements in 'same.sample' array")
  method    = match.arg(method)
  ScaleType = match.arg(ScaleType)
  nFeat     = ncol(X)
  mask      = is.na(Y)
  UnknIdx   = which( mask) # which samples are unknown to be classified
  TrainIdx  = which(!mask) # which samples are known to be used to train a classifier
  Res       = numeric(iters)
  if(n>0) group=same.sample[TrainIdx] else group=NULL
  N <- Tabl <- S <- 0
  for( iIter in 1:iters) {
    #====================================
    # Split data into test and train sets
    #====================================
    mask = sample.split(Y[TrainIdx], SplitRatio, group=group) 
    trainIdx = TrainIdx[ mask]  # split known samples into temporary train set
    testIdx  = TrainIdx[!mask]  # split known samples into temporary test  set
    xtrain   = X[trainIdx,]
    ytrain   = Y[trainIdx ]
    xtest    = X[testIdx,]
    ytest    = Y[testIdx ]
    if(n>0) ssample=same.sample[testIdx] else ssample=NULL
    if (verbose) cat(iIter,") ")
    #====================================
    # Run Classifier and calculate performance
    #====================================
    y = msc.classifier.run(xtrain, ytrain, xtest, RemCorrCol=RemCorrCol, 
        KeepCol=KeepCol, ScaleType=ScaleType, prior=prior, method=method,
        same.sample=ssample, ...)
    res = table(y, ytest, dnn=c("predicted", "true"))
    Res[iIter] = sum(diag(res))
    Tabl = Tabl + res
    N = N + diag(table(ytest, ytest))
    if (verbose) cat(Res[iIter]/length(testIdx),"\n");
  } # end of iters loop
  Res  = Res/length(testIdx)
  Tabl = Tabl/(matrix(1,nrow(Tabl),1) %*% N)
  
  #=========================================================================
  # Create Classifier based on all labeled data and run it on unknown data
  #=========================================================================
  xtrain = X[TrainIdx,]
  ytrain = Y[TrainIdx ]
  if(length(UnknIdx)>0) {
    xunkn = X[UnknIdx,] 
    if(n>0) ssample=same.sample[UnknIdx] else ssample=NULL
  } else {
    xunkn   = xtrain
    ssample = group
  }
  y = msc.classifier.run(xtrain, ytrain, xunkn, RemCorrCol=RemCorrCol, 
      KeepCol=KeepCol, ScaleType=ScaleType, prior=prior, method=method,
      same.sample=ssample, ret.prob=TRUE, ...)
  Prob = attr(y,"probabilities")
  return(list(Y=y, Res=Res, Tabl=Tabl, Prob=Prob))
}


