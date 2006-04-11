#===========================================================================#
# Written by Jarek Tuszynski. Copyright 2001-2003 SAIC.                     #
# Software developed in conjunction with the National Cancer Institute, and #
# distributed under "caBIO Software License" included in "COPYING" file.    #
#===========================================================================#

#TODO: add scanTime and retentionTime
# http://tools.proteomecenter.org/mzXMLschema.php
#source("C:/programs/R/R-2.2.1/src/library/caMassClass/R/mzXML.R")

new.mzXML = function()
{
  object = list(
    header         = NULL, # required - list    - Path to all the ancestor files (up to the native acquisition file) used to generate the current XML instance document.
    parentFile     = NULL, # required - list    - Path to all the ancestor files (up to the native acquisition file) used to generate the current XML instance document.
    dataProcessing = NULL, # required - list    - Description of any manipulation (from the first conversion to mzXML format until the creation of the current mzXML instance document) applied to the data.
    msInstrument   = NULL, # optional - element - General information about the MS instrument.
    separation     = NULL, # optional - element - Information about the separation technique, if any, used right before the acquisition.
    spotting       = NULL, # optional - element - Acquisition independent properties of a MALDI experiment.
    scan           = vector(mode="list")
  )
  class(object) <- "mzXML"
  return(object)
}

#==============================================================================

read.mzXML = function(filename)
{
  Paste = function(...) paste(..., sep="", collapse="")

  strtrunc = function(Str,Sub) {
    lp = attr(regexpr(paste(".*",Sub,sep=""),Str),'match.length')
    return( substring(Str, 1, lp) )
    #y = unlist(strsplit(Str,Sub)) # other way of doing it
    #return( paste(y[-length(y)], sub,  sep="", collapse="") )
  }

  fregexpr = function(pattern, filename)
  { # similar to gregexpr but operating on files not strings
    buf.size=1024
    n  = file.info(filename)$size
    pos = NULL
    fp = file(filename, "rb")
    for (d in seq(1,n,by=buf.size)) {
      m = if (n-d>buf.size) buf.size else n-d
      p = gregexpr(pattern, readChar(fp, m))[[1]]
      if(p[1]>0) pos=c(pos, p+d-1)
    }
    close(fp)
    if (is.null(pos)) pos=-1
    return (pos)
  }

  #-------------------------------
  # define XML handler function
  #-------------------------------
  mzXMLhandlers <- function()
  {
    #---------------------------------------------------------
    # local variables
    #---------------------------------------------------------
    obj      = new.mzXML() # create new mzXML object
    iScan    = 0
    ParentID = vector(mode="integer")
    sha1     = vector(mode="list", length=2) # optional - element - sha-1 sums
    sha1[1] <- sha1[2] <- 0
    # Optional attributes that might come with a scan that will be stored
    OptScanAttr = c("polarity", "scanType", "centroided", "deisotoped", 
      "chargeDeconvoluted", "retentionTime", "ionisationEnergy", 
      "collisionEnergy", "cidGasPressure", "totIonCurrent") 
      
    #-------------------------------
    # local functions
    #-------------------------------
    ToString = function(x, indent = 1)
    { # converts content of a node to a string
      if (is.null(x)) return(NULL);
      spaces = if (indent>0) Paste(rep("  ", indent)) else ""
      Name = xmlName(x, TRUE)
      val  = xmlValue(x)
      if (Name=="text") return( Paste(spaces, val, "\n") )
      if (!is.null(xmlAttrs(x))) {
        att = paste(names(xmlAttrs(x)), paste("\"", xmlAttrs(x),
            "\"", sep = ""), sep = "=", collapse = " ")
        att = paste(" ", att, sep="")
      } else att = ""
      chl = ""
      for (i in xmlChildren(x)) chl = Paste(chl, ToString(i, indent+1))
      if (chl=="") Paste(spaces, "<" , Name, att, "/>\n")
      else Paste(spaces, "<" , Name, att, ">\n", chl, spaces, "</", Name, ">\n")
    }

    CatNodes = function(x,Name, indent = 2)
    { # concatinate strings of several nodes
      Str=NULL
      for (y in xmlElementsByTagName(x, Name))
        Str = paste(Str, ToString(y,indent), sep="")
      return(Str)
    }

    read.mzXML.scan = function(x)
    { # process scan section of mzXML file
      if (is.null(x)) return(NULL)
      if (xmlName(x) != "scan") return(NULL)
      scanOrigin <- precursorMz <- nameValue <- maldi <- mass <- peaks <- NULL
      att         = xmlAttrs(x)
      num         = as.integer(att["num"])
      msLevel     = as.integer(att["msLevel"])
      peaksCount  = as.integer(att["peaksCount"]) # Total number of m/z-intensity pairs in the scan            
      msk         = names(att) %in% OptScanAttr
      if (sum(msk)==0) scanAttr = ""
      else {
        scanAttr = paste( names(att[msk]), paste("\"", att[msk],
                           "\"", sep = ""), sep = "=", collapse = " ")
        scanAttr = paste(" ", scanAttr, sep="")
      }
      maldi       = ToString(x[["maldi"]])
      scanOrigin  = CatNodes(x, "scanOrigin", 3)
      nameValue   = CatNodes(x, "nameValue", 3)
      precursorMz = CatNodes(x, "precursorMz", 3)
      precursorMz = gsub("\n      " , " " , precursorMz)
      for (y in xmlElementsByTagName(x, "scan"))
        ParentID[as.integer(xmlAttrs(y)["num"])] <<- num
      y           = x[["peaks"]]
      att         = xmlAttrs(y)
      peaks       = xmlValue(y) # This is the actual data encoded using base64
      precision   = att["precision"] # nr of bits used by each component (32 or 64)
      byteOrder   = att["byteOrder"] # Byte order of the encoded binary information (must be network)
      pairOrder   = att["pairOrder"] # Order of the m/z intensity pairs (must be m/z-int
      endian      = if(byteOrder=="network") "big" else "little"
      if(precision=="32") size=4
      else if(precision=="64") size=8
      else stop("read.mzXML.scan: incorrect precision attribute of peaks field")
      if (pairOrder!="m/z-int")
        warning("read.mzXML.scan: incorrect pairOrder attribute of peaks field")
      if (peaksCount>0) {
        p = base64decode(peaks, "double", endian=endian, size=size)
        np = length(p) %/% 2
        if (np != peaksCount)
          warning("read.mzXML.scan: incorrect 'peakCount' attribute of 'peaks' field: expected ",
          peaksCount, ", found ", np, "  ",(3*((nchar(peaks)*size)/4))/2, " (scan #",num,")")
        dim(p)=c(2, np)
        mass =p[1,]
        peaks=p[2,]
      }
      #x$children=NULL; # needed to capture the header
      #header <<- toString(x)
      return( list(mass=mass, peaks=peaks, num=num, parentNum=num,
        msLevel=msLevel, scanAttr=scanAttr, maldi=maldi,
        scanOrigin=scanOrigin, precursorMz=precursorMz, nameValue=nameValue) )
    }

    #---------------------------------------------------------
    # the instructions how to parse each section of mzXML file
    #---------------------------------------------------------
    list(
      mzXML  = function(x, ...) {
        y = x[["sha1"]]
        sha1[1]    <<- if (!is.null(y)) xmlValue(y) else 0
        x$children =  NULL
        obj$header <<- toString(x)
        NULL
      },

      msRun = function(x, ...) {
        y = x[["sha1"]]
        sha1[2]            <<- if (!is.null(y)) xmlValue(y) else 0
        obj$msInstrument   <<- ToString(x[["msInstrument"]],2)
        obj$separation     <<- ToString(x[["separation"]],2)
        obj$spotting       <<- ToString(x[["spotting"]],2)
        obj$parentFile     <<- CatNodes(x,"parentFile")
        obj$dataProcessing <<- CatNodes(x,"dataProcessing")
        NULL
      },

      scan  = function(x, ...) {
        iScan <<- iScan+1
        obj$scan[[iScan]] <<- read.mzXML.scan(x)
        x$children=NULL
        x
      },

      data = function() {
        if (is.null(obj$header)) NULL
        else list(mzXML=obj, ParentID=ParentID, sha1=sha1)
      }
    ) #end of list of handler functions
  } # done with local functions

  #---------------------------------
  # begining of read.mzXML function
  #---------------------------------
  require(XML)
  require(digest)
  require(caTools)
  if (!is.character(filename)) stop("read.mzXML: 'filename' has to be a string")
  if (length(filename)>1) filename = paste(filename, collapse = "")  # combine characters into a string

  sha1File = digest(filename, algo="sha1", file=TRUE)
  x = xmlTreeParse(file=filename, handlers=mzXMLhandlers(),
          addAttributeNamespaces=TRUE) $ data()
  if (is.null(x)) # is this file a mzXML file ?
    stop("read.mzXML: This is not mzXML file");
  mzXML    = x$mzXML
  sha1Read = x$sha1

  # sort scans into correct order; find parent numbers of recursive nodes
  n = length(mzXML$scan)
  NumID = integer(n)
  for (i in 1:n) {
    NumID[i] = mzXML$scan[[i]]$num
    mzXML$scan[[i]]$scanOrigin = paste("<scanOrigin parentFileID='",sha1File,
                                       "' num='",NumID[i],"'/>\n", sep="");
  }
  mzXML$scan = mzXML$scan[ order(NumID) ]
  for (i in 1:n)
    if (!is.na(x$ParentID[i])) mzXML$scan[[i]]$parentNum = x$ParentID[i]
    else x$ParentID[i] = mzXML$scan[[i]]$parentNum
  mzXML$scan = mzXML$scan[ order(x$ParentID) ]

  ## read sha1 section
  n = sum(as.integer(lapply(sha1Read, is.character))) # how many sha1 were found
  if( n>0 ) {
    ## sha1 - sha-1 sum for this file (from the beginning of the file up to
    ## (and including) the opening tag of sha1
    if (is.null(sha1Read[[1]])) sha1Read[[1]]=sha1Read[[2]]
    sha1Pos = fregexpr("<sha1>", filename) + 6 # 6 = length("<sha1>")
    for(i in n) { # multiple sha1 sections are possible
      sha1Calc = digest(filename, algo="sha1", file=TRUE, length=sha1Pos[i]-1)
      if (sha1Read[[i]]!=sha1Calc)
        warning("Stored and calculated Sha-1 sums do not match (stored '",
          sha1Read[[i]],"'; calculated '", sha1Calc,"')")
    }
  }

  # strip mzXML terminator from header section
  mzXML$header = gsub("</mzXML>", "", mzXML$header)
  mzXML$header = gsub("^ +", "", mzXML$header)
  # add info about parent file (the file we just read)
  mzXML$parentFile = Paste(mzXML$parentFile, "    <parentFile filename='file://",
    filename, "' fileType='processedData' fileSha1='", sha1File, "'/>\n")
  return( mzXML )
}

#=============================================================================

write.mzXML = function(mzXML, filename, precision=c('32', '64'))
{
  Paste  = function(...) paste(..., sep="", collapse="")

  fprintf = function(fp, level, ..., append=TRUE)
  { # helper function
    x = paste(..., sep="")
    if (length(x)==0 || is.null(x)) return(NULL)
    spaces = if (level>0) Paste(rep("  ", level)) else ""
    x = gsub("'", "\"", x)
    cat(spaces, x, file=fp, sep="")
    NULL
  }  # done with local functions

  require(XML)
  require(digest)
  precision = match.arg(precision)
  if (!is.character(filename)) stop("read.mzXML: 'filename' has to be a string")
  if (length(filename)>1) filename = paste(filename, collapse = "")  # combine characters into a string
  fp  = file(filename, "w")

  if (is.null(mzXML))
    stop("write.mzXML: Variable mzXML has to be an instance of class mzXML");
  if (attr(mzXML, "class")!="mzXML")
    stop("write.mzXML: Variable mzXML has to be an instance of class mzXML");

  #------------------------------------
  # Fill-in required sections if empty
  #------------------------------------
  if (is.null(mzXML$header)) {
    Str = "http://sashimi.sourceforge.net/schema_revision/"
    mzXML$header = Paste( "<mzXML xmlns='",Str,"mzXML_2.1'\n  ",
       "xmlns:xsi='http://www.w3.org/2001/XMLSchema-instance'\n  ",
       "xsi:schemaLocation='",Str,"mzXML_2.1 ",
       Str,"mzXML_2.1/mzXML_idx_2.1.xsd'>\n")
  }
  if (is.null(mzXML$parentFile)) {
    mzXML$parentFile = Paste( "    <parentFile filename='file://unknown' ",
    "fileType='RAWData' fileSha1='0000000000000000000000000000000000000000'/>\n")
  }
  if (is.null(mzXML$dataProcessing)) {
    Version = packageDescription("caMassClass")$Version
    Time    = format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
    mzXML$dataProcessing = Paste("    <dataProcessing>\n",
      "      <software type='processing' name='cran.r-project.org/caMassClass' ",
      "version='",Version,"' completionTime='",Time,"'/>\n    </dataProcessing>")
  }

  #-----------------------------
  # Write beggining of file
  #-----------------------------
  fprintf(fp, 0, "<?xml version='1.0' encoding='ISO-8859-1'?>\n", append=FALSE)
  fprintf(fp, 0, mzXML$header)
  fprintf(fp, 1, "<msRun scanCount='",length(mzXML$scan),"'>\n")
  fprintf(fp, 0, mzXML$parentFile)
  fprintf(fp, 0, mzXML$msInstrument)
  fprintf(fp, 0, mzXML$dataProcessing)
  fprintf(fp, 0, mzXML$separation)
  fprintf(fp, 0, mzXML$spotting)

  #---------------------------------
  # Write scan Section
  #---------------------------------
  indexScan = "  <index name='scan'>\n"
  n   = length(mzXML$scan)
  Num = integer(n)
  for (i in 1:n) Num[i] = mzXML$scan[[i]]$num
  mzXML$scan = mzXML$scan[ order(Num) ]
  for (i in 1:n) Num[i] = mzXML$scan[[i]]$parentNum
  mzXML$scan = mzXML$scan[ order(Num) ]
  for (i in 1:n) Num[i] = mzXML$scan[[i]]$msLevel
  Num = 1-diff(c(Num,1)) # number of </scan> after each scan
  size = (if (precision=="32") 4 else 8)
  for(i in 1:n) if (!is.null(mzXML$scan[[i]])) {
    indexScan = Paste(indexScan, "    <offset id='",mzXML$scan[[i]]$num,"'>",
       seek(fp),"</offset>\n")
    mass  = mzXML$scan[[i]]$mass
    peaks = mzXML$scan[[i]]$peaks
    stopifnot(length(mass)==length(peaks))
    ScanHeader = Paste("<scan num='"   , mzXML$scan[[i]]$num,
                       "' msLevel='"   , mzXML$scan[[i]]$msLevel, 
                       "' peaksCount='", length(peaks),
                       "'	lowMz='"     , min(mass),
                       "'	highMz='"    , max(mass),
                       "'", mzXML$scan[[i]]$scanAttr, ">\n")
    fprintf(fp, 2, ScanHeader)
    fprintf(fp, 3, mzXML$scan[[i]]$scanOrigin)
    fprintf(fp, 0, mzXML$scan[[i]]$precursorMz)
    fprintf(fp, 0, mzXML$scan[[i]]$maldi)
    p = as.vector(rbind(mass,peaks))
    fprintf(fp, 3, Paste("<peaks precision='",precision,
      "' byteOrder='network' pairOrder='m/z-int'>",
      base64encode(p, endian="big", size=size), "</peaks>\n"))
    fprintf(fp, 0, mzXML$scan[[i]]$nameValue)
    if(Num[i]) for (j in 1:Num[i]) fprintf(fp, 2, "</scan>\n")
  }
  indexScan = Paste(indexScan, "  </index>\n")


  #---------------------------------
  # Write end of file
  #---------------------------------
  fprintf(fp, 1, "</msRun>\n")
  n = seek(fp)
  fprintf(fp, 0, indexScan)
  fprintf(fp, 1, "<indexOffset>",n,"</indexOffset>\n")
  cat("  <sha1>", file=fp, sep="")
  n = seek(fp)
  close(fp)
  sha1 = digest(filename, algo="sha1", file=TRUE, length=n)
  cat(sha1, "</sha1>\n</mzXML>\n", file=filename, append=TRUE, sep="")
  invisible(NULL)
}

#==============================================================================

