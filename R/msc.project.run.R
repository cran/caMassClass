#===========================================================================#
# Written by Jarek Tuszynski. Copyright 2001-2003 SAIC.                     #
# Software developed in conjunction with the National Cancer Institute      #
# Distributed under GNU General Public License version 3                    #
#===========================================================================#

msc.project.run = function(ProjectFile, directory.out=NULL, verbose=TRUE,  ...)
{
  #============================================================================
  # check if DumpFile exist. If so than load it. It will contain names of all the 
  # Files Storing data in R format. If DumpFile does not exist than samples have 
  # to be converted from CSV files to R files. to do so read file named SampleFile 
  # where names of all the sample files are stored.
  #============================================================================
  directory.in = dirname(ProjectFile) # extract directory out of project file
  if (is.null(directory.out)) directory.out = directory.in;
  DumpFile   = file.path(directory.out, "RInputFiles.csv")
  if (file.exists(DumpFile)) {
    if(verbose) cat("Load precomputed input data")
    FileNames = read.csv(file=DumpFile, comment.char = "", header=FALSE)
    FileNames = as.character(FileNames[,1])
  } else {
    if(verbose) cat("Read CSV files and save them in R format\n")
    FileNames = msc.project.read(ProjectFile, directory.out)
    write.table(FileNames, file=DumpFile, sep=",", row.names=FALSE, col.names=FALSE)
  }
  
  nSet = length(FileNames)  # number of sets
  XX   = NULL    
  for (iSet in 1:nSet) {
    #============================================================================
    # load new file - 3 new variables will be loaded:
    # X - 2D matrix or 3D data cube containing all the samples (features x samples x copies )
    # SampleLabels - for each sample (row) stored in X this array 
    #                contain class labels (cancer? , Normal?, BPH?, test set?)
    # mzXML - info about experiment, (if present)
    load(FileNames[iSet]) 
    Name = gsub(sprintf("%s/Data_", directory.out),"", FileNames[iSet])
    Name = gsub(".Rdata","", Name)
    if(length(dim(X))==3) nCopy=dim(X)[3] else nCopy=1
    if(verbose) cat("Preprocess", Name,"data\n")
    
    #============================================================================
    # preprocess the data, or perform all the processing that does not need SampleLabels
    #============================================================================
    if (!is.null(mzXML) & iSet==1) {
      X = msc.preprocess.run(X, mzXML=mzXML, verbose=verbose, ...)
      mzXML = attr(X,"mzXML")
    } else X = msc.preprocess.run(X, verbose=verbose, ...)
    
    #============================================================================
    # extract unique name of the data from the filename and use it in col names
    #============================================================================
    mass = as.integer(gsub("M","", rownames(X)))  
    rownames(X) = paste(Name, mass, sep="_")
    if (nSet>1) XX = rbind(XX, X)
  }
  #============================================================================
  # return preprocessed data
  #============================================================================
  if (nSet>1) X = XX
  rm(XX)
     
  if(verbose) cat("Update sample labels\n")
  if (length(dim(X))==3) X = msc.copies.merge(X, 1)
  nC = ncol(X)
  nS = length(SampleLabels)
  SameSample   = rep(seq(nS),nC/nS)
  SampleLabels = rep(SampleLabels,nC/nS) 
  # Args = as.list(match.call()) 
  # if("merge.copies" %in% names(Args)) merge.copies=eval(Args$merge.copies)
  # else merge.copies = 1;
  # a = merge.copies %/% 4
  # b = merge.copies %%  4
  # if (nCopy==2 && a==2) b=0; 
  # if (b==1) SampleLabels = c(SampleLabels, SampleLabels)
  # if (b==3) SampleLabels = c(SampleLabels, SampleLabels, SampleLabels)
  # if (b==1) SameSample   = c(SameSample, SameSample)
  # if (b==3) SameSample   = c(SameSample, SameSample, SameSample)
  
  return (list(X=X, SampleLabels=SampleLabels, SameSample=SameSample, mzXML=mzXML))
}
