#===========================================================================#
# Written by Jarek Tuszynski. Copyright 2001-2003 SAIC.                     #
# Software developed in conjunction with the National Cancer Institute, and #
# distributed under "caBIO Software License" included in "COPYING" file.    #
#===========================================================================#

msc.project.read = function(ProjectFile, directory.out = NULL) 
{
  trim = function(Str) 
   {sub('^[[:space:]]+', '', sub('[[:space:]]+$','',as.character(Str)))} 
  #===================
  # Read Index file
  #===================
  directory.in = dirname(ProjectFile) # extract directory out of project file
  if (is.null(directory.out)) directory.out = directory.in;
  if (!file.exists(ProjectFile)) stop(sprintf("Can not find %s file.", ProjectFile))
  FileList = read.csv(file=ProjectFile, comment.char = "")

  nSamp = nrow(FileList)
  nCopy = ncol(FileList)-2
  SampleNames  = trim(FileList[,1])
  SampleLabels = trim(FileList[,2])
  FileList = FileList[,(1:nCopy)+2]
  ColumnNames  = trim(colnames(FileList))

  if (nCopy==1) {
    X = msc.msfiles.read.csv(directory.in, FileList, SampleNames)
    FileNames = file.path(directory.out, sprintf("Data_%s.Rdata", ColumnNames[1]))
    save(X, SampleLabels, file=FileNames, compress=TRUE)
  } else {
    #=======================
    # interpret column names
    #=======================  
    ColStr = gsub("[0-9]","", ColumnNames)          # delete all numbers
    ColNum = gsub("[a-z]","", tolower(ColumnNames)) # delete all letters
    UniqueStr = unique(ColStr)
    UniqueNum = unique(ColNum)
    nSets  = length(UniqueStr)
    M = matrix(0, nSets, length(ColNum))
    for (i in 1:nCopy) {
      Row = pmatch(ColStr[i], UniqueStr)
      Col = pmatch(ColNum[i], UniqueNum)
      M[Row, Col] = i
    }

    #===============================
    # Read and save each set of data
    #===============================
    FileNames = UniqueStr
    for (i in 1:nSets) {
      Cols = M[i,]
      Cols = Cols[Cols>0]
      X = msc.msfiles.read.csv(directory.in, FileList[,Cols], SampleNames, UniqueNum)
      FileNames[i] = file.path(directory.out, sprintf("Data_%s.Rdata", UniqueStr[i]))
      save(X, SampleLabels, file=FileNames[i], compress=TRUE)
    }
    return (FileNames)    
  }
}

