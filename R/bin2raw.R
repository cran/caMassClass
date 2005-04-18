bin2raw = function(x, size=NA)
{
  if (typeof(x)=="raw") return(x)
  TypeList = c("logical", "integer", "double", "complex", "character", "raw")
  if (is.character(x) && length(x)==1) x = strsplit(x, NULL)[[1]]
  if (!is.vector(x) || mode(x) == "list") 
     stop("bin2raw: can only write vector objects")
  if (!is.na(size)) nBits=size
  else nBits = switch(match(typeof(x), TypeList), 4, 4, 8, 16, 2, 1)
  n = length(x)
  fname = tempfile() 
  writeBin(x, fname, size=size)
  r = readBin(fname, "raw", n = n*nBits)
  file.remove(fname)
  return (r)
}

raw2bin = function(r, what, size=NA, signed = TRUE)
{
  TypeList = c("logical", "integer", "double", "complex", "character", "raw", 
               "numeric", "int")
  if (!is.character(what) || length(what) != 1 || !(what %in% TypeList)) 
    what <- typeof(what)
  if (what=="raw") return(r)
  if (!is.na(size)) nBits=size 
  else nBits = switch(match(typeof(x), TypeList), 4, 4, 8, 16, 2, 1, 8, 4) 
  n = length(r)
  fname = tempfile()      
  writeBin(as.raw(r), fname)
  x = readBin(fname, what, n = n%/%nBits, size=size, signed=signed)
  if (what=="character")  x = paste(x, collapse = "")
  file.remove(fname)
  return (x)
}


