#===========================================================================#
# Written by Jarek Tuszynski. Copyright 2001-2003 SAIC.                     #
# Software developed in conjunction with the National Cancer Institute, and #
# distributed under "caBIO Software License" included in "COPYING" file.    #
#===========================================================================#

msc.preprocess.run = function( X,
    baseline.removal = 0,
      breaks=200, qntl=0, bw=0.005,                    # bslnoff
    min.mass = 3000,                                   # msc.mass.cut
    mass.drift.adjustment = 1,
      shiftPar=0.0005,                                 # msc.mass.adjust
    peak.extraction = 0, 
     PeakFile=0, SNR=2, span=c(81,11), zerothresh=0.9, # msc.peaks.find
     BmrkFile=0, BinSize=c(0.002, 0.008), tol=0.97,    # msc.peaks.align 
     FlBmFile=0, FillType=0.9,                         # msc.biomarkers.fill
    merge.copies = 4+2+1,                              # msc.copies.merge
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
  #print(match.call(expand.dots=TRUE))
  
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
    if(verbose) cat("Baseline Removal\n")
    X = msc.baseline.subtract(X, breaks=breaks, qntl=qntl, bw=bw, method="loess")
  } else if(verbose) cat("Baseline Removal - skipped\n")
  
  # =======================================================
  # remove features/columns with no data and Cut Low Masses
  # =======================================================
  if(verbose) cat("Cut low masses\n")
  X = msc.mass.cut(X, min.mass)

  # =====================
  # mass drift adjustment
  # =====================
  if (mass.drift.adjustment) {
    if(verbose) cat("Mass Drift Adjustment\n")
    X = msc.mass.adjust(X, scalePar=mass.drift.adjustment, AvrSamp=0, 
        shiftPar=shiftPar)
  } else if(verbose) print("Mass Drift Adjustment - skipped\n")
 
  # ===============
  # peak extraction
  # ===============
  if (peak.extraction) {  # extract and align peaks
    if(verbose) cat("Peak Extraction\n")
    Peaks = msc.peaks.find(X, PeakFile=PeakFile, SNR=SNR, span=span, 
            zerothresh=zerothresh)
    if(verbose) cat("Peak Aligment\n")
    out = msc.peaks.align(Peaks, BmrkFile=BmrkFile, BinSize=BinSize, tol=tol) 
    if(verbose) cat("Fill Biomarkers\n")
    X   = msc.biomarkers.fill( X, out$Bmrks, out$BinBounds, BmrkFile=FlBmFile, 
          FillType=FillType)
  } else if(verbose) print("Peak Extraction/Aligment - skipped\n")
  
  # ==============
  # combine copies
  # ==============
  d = dim(X)
  if (length(d)==3 && d[3]>1 && merge.copies) {
    if(verbose) cat("Merge Samples\n")
    X=msc.copies.merge(X, merge.copies, PeaksOnly=!peak.extraction)
  } else if(verbose) print("Merge Samples - skipped\n")
  if(verbose) cat("Done with preprocessing\n")

  return( X )
}