#===========================================================================#
# Written by Jarek Tuszynski. Copyright 2001-2003 SAIC.                     #
# Software developed in conjunction with the National Cancer Institute, and #
# distributed under "caBIO Software License" included in "COPYING" file.    #
#===========================================================================#

msc.preprocess.run = function( X, mzXML=NULL,
    baseline.removal = 0,
      breaks=200, qntl=0, bw=0.005,                    # bslnoff
    min.mass = 3000,                                   # msc.mass.cut
    mass.drift.adjustment = 1,
      shiftPar=0.0005,                                 # msc.mass.adjust
    peak.extraction = 0, 
     PeakFile=0, SNR=2, span=c(81,11), zerothresh=0.9, # msc.peaks.find
     BmrkFile=0, BinSize=c(0.002, 0.008), tol=0.97,    # msc.peaks.align 
     FlBmFile=0, FillType=0.9,                         # msc.biomarkers.fill
    merge.copies = 1,                                  # msc.copies.merge
    verbose = TRUE) 
# Input Parameters:
#  X:  the data in one ( when X is a matrix [nSamp x nFeat]) or multiple copies  (X is an array [nSamp x nFeat x nCopy])
#  SoN, span, sm.span : getPeaks parameters used to customize peak extraction
#  baseline.removal: 0 - no Baseline removal, 1 - baseline removal using PROcess bslnoff call
#  mass.drift.adjustment: 0 - skip all
#                       1 - mass drift adjustment and normalization by matching means
#                       2 - mass drift adjustment and normalization by matching means and mediums
# merge.copies : if copies are to be concatinated (as separate samples) set merge.copies to 1 otherwise set it to 0
#               if copies are to be averaged add 2
#               if worst copy is to be deleted add 8 or 
#               if in case of a single bad copy of a sample one wants to replace it with the best copy add 4
{  
  Record = function(mzXML, verbose, Type, Name, ...) {
    if (verbose) cat(paste(Type,"\n"))
    End = "</dataProcessing>"
    Str = mzXML$dataProcessing
    lp  = attr(regexpr(paste(".*",End,sep=""),Str),'match.length')-nchar(End)
    mzXML$dataProcessing = paste(substring(Str, 1, lp-1),
      "   <processingOperation type='",Type,"' name='",Name,
      "' value='",...,"'/>\n    </dataProcessing>\n", sep="") # insert string
    return(mzXML)
  }
  
  if (is.null(mzXML)) mzXML = new.mzXML()
  mzXML$scan=NULL;
  Version = packageDescription("caMassClass")$Version
  Time    = format(Sys.time(), "%Y-%m-%dT%H:%M:%S")
  mzXML$dataProcessing = paste(mzXML$dataProcessing, 
      "    <dataProcessing>\n      <software type='processing'",
      " name='cran.r-project.org/caMassClass' version='", Version,
      "' completionTime='",Time,"'/>\n    </dataProcessing>\n", sep="")
  
  # X is expected to be a 3D array if it is a matrix than create extra dimention
  if (!is.array(X) && is.matrix(X)) {
    dNames = dimnames(X)
    dim(X) = c(nrow(X), ncol(X),1) 
    rownames(X) = dNames[[1]]
    colnames(X) = dNames[[2]]
  }

  # ================
  # BaseLine Removal
  # ================
  if (baseline.removal) {
    mzXML = Record(mzXML, verbose,  "Baseline Removal", "msc.baseline.subtract",
           "(breaks=",breaks,", qntl=",qntl,", bw=",bw,", method=loess)")
    X = msc.baseline.subtract(X, breaks=breaks, qntl=qntl, bw=bw, method="loess")
  } else if(verbose) print("Baseline Removal - skipped")
  
  # =======================================================
  # remove features/columns with no data and Cut Low Masses
  # =======================================================
  if(min.mass>0) 
    mzXML = Record(mzXML, verbose, "Cut low masses", "msc.mass.cut", 
           "(min.mass=",min.mass,")")
  X = msc.mass.cut(X, min.mass)

  # =====================
  # mass drift adjustment
  # =====================
  if (mass.drift.adjustment) {
    mzXML = Record(mzXML, verbose, "Mass Drift Adjustment", "msc.mass.adjust",
    "(scalePar=",mass.drift.adjustment,", AvrSamp=0, shiftPar=",shiftPar,")")
    X = msc.mass.adjust(X, scalePar=mass.drift.adjustment, AvrSamp=0, 
        shiftPar=shiftPar)
  } else if(verbose) print("Mass Drift Adjustment - skipped")
 
  # ===============
  # peak extraction
  # ===============
  if (peak.extraction) {  # extract and align peaks
    
    extension = function( pattern, text) 
      return(length(grep(pattern, text, ignore.case=TRUE ))>0)
    
    mzXML = Record(mzXML, verbose, "Peak Extraction", "msc.peaks.find",
       "(SNR=",SNR,", span=(",paste(span, collapse=","),"), zerothresh=",zerothresh,")")
    Peaks = msc.peaks.find(X, SNR=SNR, span=span, zerothresh=zerothresh)
    if (is.character(PeakFile)) {
     if (extension(".csv", PeakFile))
	     msc.peaks.write.csv(Peaks, PeakFile)  
     if (extension("xml", PeakFile))
	     msc.peaks.write.mzXML(Peaks, PeakFile, mzXML)  
    }
    
    mzXML = Record(mzXML, verbose, "Peak Aligment", "msc.peaks.align", 
      "(BmrkFile=",BmrkFile,", BinSize=(",paste(BinSize, collapse=","),"), tol=",tol,")")
    out = msc.peaks.align(Peaks, BinSize=BinSize, tol=tol)
    if (is.character(BmrkFile)) {
      if (extension(".csv", BmrkFile))
        msc.biomarkers.write.csv(out$Bmrks, BmrkFile)
      if (extension("xml", BmrkFile)) 
	      msc.rawMS.write.mzXML(out$Bmrks, PeakFile, mzXML)  
    }
    
    mzXML = Record(mzXML, verbose, "Fill Biomarkers", "msc.biomarkers.fill",
     "(FillType=",FillType,")")
    X = msc.biomarkers.fill( X, out$Bmrks, out$BinBounds, FillType=FillType)
    if (is.character(BmrkFile)) {
      if (extension(".csv", FlBmFile))
        msc.biomarkers.write.csv(X, FlBmFile)
      if (extension("xml", FlBmFile))
	      msc.rawMS.write.mzXML(out$Bmrks, FlBmFile, mzXML) 
    }
  } else if(verbose) print("Peak Extraction/Aligment - skipped")
  
  # ==============
  # combine copies
  # ==============
  d = dim(X)
  if (length(d)==3 && d[3]>1 && merge.copies) {
    mzXML = Record(mzXML, verbose, "Merge Samples", "msc.copies.merge",
    "(mergeType=",merge.copies,", PeaksOnly=",!peak.extraction,")")
    X=msc.copies.merge(X, merge.copies, PeaksOnly=!peak.extraction)
  } else if(verbose) print("Merge Samples - skipped")
  if(verbose) print("Done with preprocessing")
  
  attr(X,"mzXML") <- mzXML
  return( X )
}