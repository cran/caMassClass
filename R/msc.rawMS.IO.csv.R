#===========================================================================#
# Written by Jarek Tuszynski. Copyright 2001-2003 SAIC.                     #
# Software developed in conjunction with the National Cancer Institute, and #
# distributed under "caBIO Software License" included in "COPYING" file.    #
#===========================================================================#

msc.rawMS.read.csv = function(directory=".", FileList="\\.csv", 
                     mzXML.record=FALSE) 
{ # read multiple MS files listed in FileList from given directory and save 
  # each file as a column in matrix. Use SampNames as row names of the matrix
  
    trim = function(Str) # trim white spaces
      sub('^[[:space:]]+', '', sub('[[:space:]]+$','',as.character(Str))) 

    file.conn = function(path, fname)
    { # extract file name. path is strictly a file directory name. fname can be
      # a file name ( "a.csv" - trivial case), file within zipped file 
      # (b.zip/a.csv), gziped file (a.csv.gz)
      fname = trim(fname) #no spaces
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
  if (mzXML.record) mzXML = new.mzXML() 
  fconn = file.conn(directory, FileList[1])
  out   = read.csv (file=fconn, comment.char = "")
  Mass  = out[,1]
  mask  = (Mass>0)
  nFeat = sum(mask)
  nSamp = length(FileList)
  X     = matrix(0, nFeat, nSamp)

  for (j in 1:nSamp) {
    fn    = trim(FileList[j])
    fconn = file.conn(directory, fn)
    fname = file.path(directory, fn)
    out   = read.csv(file=fconn, comment.char = "")
    X[,j] = out[mask,2]
    if (mzXML.record) {
      if (regexpr("zip/",fn)>0) sha1 = digest(out, alg="sha1")
      else sha1 = digest(fname, alg="sha1", file=TRUE)
      mzXML$parentFile = paste(mzXML$parentFile, 
        "    <parentFile filename='file:///", fname, 
        "' fileType='RAWData' fileSha1='",sha1,"'/>\n", sep="")
      scanOrigin = paste("<scanOrigin parentFileID='",sha1,"' num='1'/>", sep="");
      mzXML$scan[[j]] = list(mass=NULL, peaks=NULL, num=j, parentNum=j, msLevel=1,
         header=NULL, maldi=NULL, scanOrigin=scanOrigin,precursorMz=NULL, nameValue=NULL)
    }
  }
  rownames(X) = signif(Mass[mask], 6)
  colnames(X) = FileList
  if (mzXML.record) attr(X,"mzXML") = mzXML
  return( X )
}

