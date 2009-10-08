#===========================================================================#
# Written by Jarek Tuszynski. Copyright 2001-2003 SAIC.                     #
# Software developed in conjunction with the National Cancer Institute, and #
# distributed under "caBIO Software License" included in "COPYING" file.    #
#===========================================================================#

msc.classifier.run = function( xtrain, ytrain, xtest, ret.prob=FALSE, 
                  RemCorrCol=0, KeepCol=0, prior=1, same.sample=NULL,
                  ScaleType=c("none", "min-max", "avr-std", "med-mad"),
                  method=c("svm", "nnet", "lda", "qda", "LogitBoost", "rpart"), ...) 
{ 
  if (length(ytrain)!=nrow(xtrain)) 
    stop("msc.classifier.run: Number of Samples (nrow(xtrain)) does not mach number of labels in ytrain")
  if (ncol(xtest)!=ncol(xtrain)) 
    stop("msc.classifier.run: Number of features/columns in  xtrain does not match the number in xtest")
  method    = match.arg(method)
  ScaleType = match.arg(ScaleType)
  #=============================================
  # Remove highly correlated columns and columns with low AUC
  #=============================================
  if (RemCorrCol | KeepCol) {
    idx = msc.features.select( xtrain, ytrain, 
                               RemCorrCol=RemCorrCol, KeepCol=KeepCol)
    xtrain = xtrain[,idx]
    xtest  = xtest [,idx]
  }
  #=============================================
  # Scale each feature separatly
  #=============================================
  if (ScaleType!="none") {
    x = msc.features.scale(xtrain, xtest, type=ScaleType)
    xtrain = x$xtrain
    xtest  = x$xtest
  }
  #=============================================
  # Input/Output probabilities
  #=============================================
  if (length(prior)==1) {
    p = table(ytrain)
    PriorProb = if(prior==1) p/sum(p) else p/p
  } else PriorProb=prior
  nFeat = ncol(xtrain)
  get.prob = (ret.prob || length(same.sample)==nrow(xtest))
  ytrain = as.factor(ytrain)
  #=============================================
  # Classification algorithms
  #=============================================
  if (method=="LogitBoost") {
    model = LogitBoost(xtrain, ytrain, ...)
    predY = predict(model, xtest)
    if(get.prob) Prob = predict(model, xtest, type="raw")
  } else if (method=="svm") {
    library(e1071)
    model = svm(xtrain, y=ytrain, probability=get.prob, ...)
    predY = predict(model, xtest, probability=get.prob)
    if(get.prob) Prob = attr(predY,"probabilities")
  } else if (method=="nnet") {
    library(nnet)
    if (!exists("size" , mode="numeric")) size =2 # add default
    if (!exists("trace", mode="numeric")) trace=FALSE # change default
    ytrain = class.ind(ytrain)
    model = nnet(xtrain, ytrain, trace=trace, size=size, ...)
    out   = predict(model, xtest)
    predY = colnames(ytrain) [max.col(out)]
    if(get.prob) Prob =  predict(model, xtest, type="raw")
  } else if (method=="lda") {
    library(MASS)
    model = lda(xtrain, ytrain, ...)
    out   = predict(model, xtest)
    predY = out$class
    if(get.prob) Prob = out$posterior
  } else if (method=="qda") {
    library(MASS)
    model = qda(xtrain, ytrain, ...)
    out   = predict(model, xtest)
    predY = out$class
    if(get.prob) Prob = out$posterior
  } else if (method=="rpart") {
    library(rpart)
    model = rpart( ytrain~., data = data.frame(cbind(ytrain,xtrain)), ...)
    out   = predict(model, newdata=xtest, type="class")
    if(get.prob) Prob =  predict(model, newdata=xtest, type="prob")
  }
  #=================================================================
  # in case of multiple copies of the same sample - synchronize them
  #=================================================================
  if (length(same.sample)==nrow(xtest)) {
    PostProb = Prob
    lablist = sort(unique(ytrain))
    nC = ncol(Prob)
    rgroup = split(seq(nrow(Prob)),same.sample) # split row idx into groups
    for(i in 1:length(rgroup)) { # for each group containing same.sample rows
      ridx = rgroup[[i]]  # find all rows attributed to the same sample
      nR = length(ridx)
      if (nR>1) {
        meanC = colMeans(Prob[ridx,]) # average all rows
        Prob[ridx,] = matrix(rep(meanC,each=nR),nR,nC) # copy back
      }
    }
    ord   = t(apply(-Prob, 1, order))# find order of sorted Probs 
    predY  = lablist[ord[,1]]        # find label with highest Prob
    Prob = -t(apply(-Prob,1,sort))   # sort Probs
    predY[Prob[,1]==Prob[,2]] = NA   # in case of ties return NA's
    Prob = PostProb                  # restore probabilities
  }
  if (ret.prob) attr(predY,"probabilities") = Prob
  return(predY)
}

