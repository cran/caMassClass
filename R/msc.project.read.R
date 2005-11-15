#===========================================================================#
# Written by Jarek Tuszynski. Copyright 2001-2003 SAIC.                     #
# Software developed in conjunction with the National Cancer Institute, and #
# distributed under "caBIO Software License" included in "COPYING" file.    #
#===========================================================================#

msc.project.read = function(ProjectFile, directory.out = NULL ) 
{
  
    trim = function(Str) # trim white spaces
      sub('^[[:space:]]+', '', sub('[[:space:]]+$','',as.character(Str)))
   
    rawMS.read = function(path, FileList, SampleNames=NULL)  
    { # read multiple MS files listed in FileList from given 
      # directory and save each file as a column in matrix.
      d = dim(FileList)
      nCopy = if (length(d)>1) ncol(FileList) else 1
      
      if (nCopy==1) { # single column: read it from XML or CSV files
        nSamp = length(FileList);
        FileList = trim(FileList)
        nXML  = length(grep("xml/", FileList, ignore.case=TRUE ))
        nCSV  = length(grep(".csv", FileList, ignore.case=TRUE ))
        if (nSamp==nXML) {# all samples will come from a mzXML file
          sg = unlist(strsplit(FileList, "/")) # split into filename  and scan number
          dim(sg) = c(2, nSamp) # filenames will be in raw 1; scan #'s will be in raw 2
          if (length(unique(sg[1,]))>1)
            stop("msc.project.read: Problem with 'ProjectFile': all samples have to come from the same mzXML file") 
          path = file.path(path,sg[1,1]) # get full filename of mzXML file
          X = msc.rawMS.read.mzXML(path, as.integer(sg[2,]))
        } else if (nSamp==nCSV) { # all samples will come from CSV files
          X = msc.rawMS.read.csv(path, FileList, mzXML.record=TRUE) # read them all
        } else 
          stop("msc.project.read: Problem with 'ProjectFile': all samples have to come either from CSV or from mzXML file") 
        colnames(X) = FileList # store file names as sample names, most likelly it will be over written by SampleNames
      } else { # more than one column: call recursivly to work on one column at a time
        Y = rawMS.read(path, FileList[,1])
        X = array(0, c(nrow(Y), ncol(Y), nCopy))
        X[,,1] = Y                # save first copy of data
        mzXML = attr(Y,"mzXML")   # save corresponding metadata
        for (i in 2:nCopy) {
           Y = rawMS.read(path, FileList[,i])
           X[,,i] = Y             # save other copies
           mzXML1 = attr(Y,"mzXML") # this metadata is lost in the current version
           mzXML$parentFile = paste(mzXML$parentFile, mzXML1$parentFile)
           # TO DO: merge mzXML scan's
        }
        mzXML$scan = vector(mode="list") # erase scan data
        rownames(X) = rownames(Y) # mass (m/z)
        if (length(SampleNames)==ncol(X)) colnames(X) = SampleNames
        attr(X,"mzXML") = mzXML
      }
      return(X)
    }
   
  #===================
  # Read Index file
  #===================
  if (is.character(ProjectFile) & length(ProjectFile)==1) {
    directory.in = dirname(ProjectFile) # extract directory out of project file
    if (!file.exists(ProjectFile))
      stop(sprintf("Can not find %s file.", ProjectFile))
    FileList = read.csv(file=ProjectFile, comment.char = "")
  } else {
    FileList = ProjectFile
    directory.in = '.'
  }
  if (is.null(directory.out)) directory.out = directory.in;

  nSamp = nrow(FileList)
  nCopy = ncol(FileList)-2
  SampleNames  = trim(FileList[,1]) # name attached to each sample
  SampleLabels = trim(FileList[,2]) # labels like: cancer vs. normal
  FileList = FileList[,(1:nCopy)+2]
  ColumnNames  = trim(colnames(FileList))
  mzXML = NULL
  
  if (nCopy==1) {
    X = rawMS.read(directory.in, trim(FileList), SampleNames)
    mzXML = attr(X,"mzXML");  attr(X,"mzXML") = NULL;
    FileNames = file.path(directory.out, sprintf("Data_%s.Rdata", ColumnNames[1]))
    save(X, SampleLabels, mzXML, file=FileNames, compress=TRUE)
  } else {
    #=======================
    # interpret column names
    #=======================  
    ColStr = gsub("[0-9]","", ColumnNames)          # delete all numbers
    ColNum = gsub("[a-z]","", tolower(ColumnNames)) # delete all letters
    UniqueStr = unique(ColStr) # signify different batches (sets)
    UniqueNum = unique(ColNum) # signify different copies in the same batch
    nSets  = length(UniqueStr) # number of sets
    M = matrix(0, nSets, length(ColNum))
    for (i in 1:nCopy) {
      Row = pmatch(ColStr[i], UniqueStr) # which set ?
      Col = pmatch(ColNum[i], UniqueNum) # which copy?
      M[Row, Col] = i
    }

    #===============================
    # Read and save each set of data
    #===============================
    FileNames = UniqueStr
    for (i in 1:nSets) {
      Cols = M[i,]        # all columns belonging to set number 'i'
      #Cols = Cols[Cols>0]
      X = rawMS.read(directory.in, FileList[,Cols], SampleNames)
      mzXML = attr(X,"mzXML");  attr(X,"mzXML") = NULL;
      dimnames(X) = list(rownames(X), colnames(X), UniqueNum)
      FileNames[i] = file.path(directory.out, sprintf("Data_%s.Rdata", UniqueStr[i]))
      save(X, SampleLabels, mzXML, file=FileNames[i], compress=TRUE)
    }
  } # end if else
  return (FileNames)    
}

