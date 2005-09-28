#===========================================================================#
# Written by Jarek Tuszynski. Copyright 2001-2003 SAIC.                     #
# Software developed in conjunction with the National Cancer Institute, and #
# distributed under "caBIO Software License" included in "COPYING" file.    #
#===========================================================================#

msc.rawMS.read.csv = function(directory=".", FileList="\\.csv") 
{ # read multiple MS files listed in FileList from given directory and save 
  # each file as a column in matrix. Use SampNames as row names of the matrix
  
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
  
  #=============================================
  # convert input parameters into uniform format
  #=============================================
  directory = as.character(directory)
  if (length(FileList)==1) { # FileList is a single string
    FileList = list.files(directory, pattern=FileList)
    if (length(FileList)==0) 
      stop("msc.rawMS.read.csv: No files in the directory matches given pattern.")
  } # else FileList is a list - list of file names to be read   
   
  #=============================================
  # Read input files
  #=============================================
  fn    = FileName(directory, FileList[1])
  out   = read.csv(file=fn, comment.char = "")
  Mass  = out[,1]
  mask  = (Mass>0)
  nFeat = sum(mask)
  nSamp = length(FileList)
  X     = matrix(0, nFeat, nSamp)
  X[,1] = out[mask,2]

  for (j in 2:nSamp) {
     fn = FileName(directory, FileList[j])
     out = read.csv(file=fn, comment.char = "")
     X[,j] = out[mask,2]
  }
  rownames(X) = signif(Mass[mask], 6)
  colnames(X) = FileList
  return( X )
}

