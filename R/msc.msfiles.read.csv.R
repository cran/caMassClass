#===========================================================================#
# Written by Jarek Tuszynski. Copyright 2001-2003 SAIC.                     #
# Software developed in conjunction with the National Cancer Institute, and #
# distributed under "caBIO Software License" included in "COPYING" file.    #
#===========================================================================#

msc.msfiles.read.csv = function(directory=".", FileList="\\.csv", 
                 SampleNames=NULL, CopyNames=NULL) 
{ # read multiple MS files in multiple copies listed in FileList 
  # from given directory and save each file as a column in matrix. 
  # If multiple copies are present (when FileList have multiple columns)
  # then save them as "planes" in 3D data cube. Use SampNames
  # as row names, and CopyNames as names of planes of the matrix
  
    FileName = function(path, fname)
    { # extract file name. path is strictly a file directory name. fname can be
      # a file name ( "a.csv" - trivial case), file within zipped file 
      # (b.zip/a.csv), gziped file (a.csv.gz)
      fname = as.character(fname)
      fname = sub('^[[:space:]]+', '', sub('[[:space:]]+$','',fname)) #no spaces
      if (regexpr("zip/",fname)>0) {
        sg = strsplit(fname, "/")[[1]]
        fn = unz( file.path(path,sg[1]), file=sg[2])
      } else {
        fn = file.path(path,fname)
        if ("gz" %in% strsplit(fn,split="\\.")[[1]]) fn = gzfile(fn)
      }
      return (fn)
    }
  
    ReadMSFiles = function(path, FileList) 
    { # read multiple MS files listed in FileList from given 
      # directory and save each file as a column in matrix.
      fn    = FileName(path,FileList[1])
      out   = read.csv(file=fn, comment.char = "")
      Mass  = out[,1]
      mask  = (Mass>0)
      nFeat = sum(mask)
      nSamp = length(FileList)
      X     = matrix(0, nFeat, nSamp)
      X[,1] = out[mask,2]
    
      for (j in 2:nSamp) {
         fn = FileName(path,FileList[j])
         out = read.csv(file=fn, comment.char = "")
         X[,j] = out[mask,2]
      }
      rownames(X) = signif(Mass[mask], 6)
      return( X )
    }

  directory = as.character(directory)
  if (is.data.frame(FileList)) nCopy = ncol(FileList)
  else {
    nCopy=1
    if (length(FileList)==1) {
      FileList = list.files(directory, pattern=FileList)
      if (length(FileList)==0) 
        stop("msc.msfiles.read.csv: No files in the directory matches given pattern.")
    }
  }
  if (nCopy==1) {
    X = ReadMSFiles(directory, FileList)
    if (length(SampleNames)==ncol(X)) colnames(X) = SampleNames
    else colnames(X) = FileList
  } else {
    Y = ReadMSFiles(directory, FileList[,1])
    X = array(0, c(nrow(Y), ncol(Y), nCopy))
    X[,,1] = Y
    for (i in 2:nCopy) 
      X[,,i] = ReadMSFiles(directory, FileList[,i])
    rownames(X) = rownames(Y)
    if (length(CopyNames)==nCopy) dimnames(X) = list(rownames(X), colnames(X), CopyNames)
    if (length(SampleNames)==ncol(X)) colnames(X) = SampleNames
  }
  return( X )
}

