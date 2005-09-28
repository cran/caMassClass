#===========================================================================#
# Written by Jarek Tuszynski. Copyright 2001-2003 SAIC.                     #
# Software developed in conjunction with the National Cancer Institute, and #
# distributed under "caBIO Software License" included in "COPYING" file.    #
#===========================================================================#

msc.rawMS.read.mzXML = function(input, scanIdx=NULL)
{
  a = attr(input,"class")
  flag = FALSE
  if (!is.null(a))
    if (a=="mzXML") flag=TRUE 
  mzXML = if(flag) input else read.mzXML(input) 
  nScan = length(mzXML$scan)
  if (is.null(scanIdx)) { # scanIdx not provided - return all scans
    j = 1
    mScan = nScan 
  } else {                # scanIdx provided - return selected scans
    j = scanIdx[1]
    mScan = length(scanIdx) 
  }
  if (j<=0 || j>nScan) stop("msc.rawMS.read.mzXML: wrong number in 'scanIdx'") 
  mass = mzXML$scan[[j]]$mass
  nMass = length(mass)
  
  Data = array(dim=c(length(mass), mScan)) # store as 2D array
  cNames = integer(mScan)
  for (i in 1:mScan) {
    if (!is.null(scanIdx)) { # scanIdx provided - return selected scans
      j = scanIdx[i]
      if (mzXML$scan[[j]]$num!=j) {# most likelly this is true, but if ...
        m = 0                      # ... is not than look for scan number 'j'
        for (k in 1:nScan) 
          if (j==mzXML$scan[[k]]$num) { m=k; break; }
        j = m   
      }  
      if (j<=0 || j>nScan) stop("msc.rawMS.read.mzXML: wrong number in 'scanIdx'")
    } else j=i # scanIdx not provided - return all scans
    SameMass = (nMass==length(mzXML$scan[[j]]$mass)) # is mass the same?
    if ( SameMass) SameMass = (sum(mass-mzXML$scan[[j]]$mass)==0)
    if (!SameMass) 
     stop("msc.rawMS.read.mzXML: spectras in mzXML files are of different lengths") 
    Data[,i]  = mzXML$scan[[j]]$peaks # store peak intensities (heights)
    cNames[i] = mzXML$scan[[j]]$num
  }
  colnames(Data) = as.character(cNames)
  rownames(Data) = signif(mass, 6)
  if (!is.ordered(cNames)) Data = Data[,order(cNames)] # order samples according to sample number
  return( list(scans=Data, mzXML=mzXML) )
}

#==============================================================================

msc.rawMS.write.mzXML = function(scans, filename, mzXML=NULL, ...)
{
  if (is.null(mzXML)) mzXML = new.mzXML()
  
  # add caMassClass signature to mzXML data
  if (is.null(mzXML$dataProcessing) || regexpr("caMassClass", mzXML$dataProcessing)<0) {
    Version = packageDescription("caMassClass")$Version
    Time    = format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
    mzXML$dataProcessing = paste(mzXML$dataProcessing, 
      "    <dataProcessing>\n      <software type='processing'",
      " name='cran.r-project.org/caMassClass' version='", Version,
      "' completionTime='",Time,"'/>\n    </dataProcessing>\n", sep="")
  }
  
  Mass  = as.numeric(rownames(scans))
  scans[is.na(scans)] = 0                  # convert NA to zeros
  # Convert data stored as 2D or 3D array 'X' to mzXML$scan format
  d = dim(scans)                     # dimentions of the data
  if (length(d)==3)                 # if this is 3D data than convert it to 2D
    dim(scans) = c(d[1], prod(d)/d[1])   # and resize the data 

  stopifnot(length(Mass)== nrow(scans)) 
  nScan = ncol(scans)
  mzXML$scan = vector(mode="list", length=nScan)
  for(i in 1:nScan) {
     mzXML$scan[[i]] = list(mass=Mass, peaks=scans[,i], num=i, parentNum=i,
      msLevel=1, maldi=NULL, scanOrigin=NULL, nameValue=NULL, precursorMz=NULL) 
  } 
  # save the file
  write.mzXML(mzXML, filename, ...)
}
