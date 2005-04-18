#===========================================================================#
# Written by Jarek Tuszynski. Copyright 2001-2003 SAIC.                     #
# Software developed in conjunction with the National Cancer Institute, and #
# distributed under "caBIO Software License" included in "COPYING" file.    #
#===========================================================================#

msc.peaks.clust = function(dM, S, BinSize=c(0,sum(dM)), tol=0.97, verbose=FALSE) 
 # bin = DivCon(dM,S). Aligne peaks into bins in such a way as to make sure 
 # that no two peaks from the same sample are present in the same bin. Bins are split 
 # at the largest gap between bins, or if there are multiple similar size gaps than in 
 # such a way as to minimize number of repeats on each smaller bin.
 # Takes 2 arrays of equal size: 
 # dM - stores distances between sorted peaks
 # S  - stores peak sample number. 
 # tol - gaps bigger than tol*max(gap) are assumed to be the same size as the largest gap
 # BinSize - upper and lower bound of bin size measured as (R-L)/mean(R,L)
 # verbose - print out degguging info
 # The output is binary array where left boundaries of each clusters-bin (biomarker) are marked
 {
  #repeats = function(x) return (sum(diff(sort(x))==0)) # how many numbers repeat themselves in vector x?
   repeats = function(x) return (sum( duplicated(x) ) ) # how many numbers repeat themselves in vector x?
   
   nS = length(S)
   if (nS!=length(dM) & (nS-1)!=length(dM)) 
     stop("Error in msc.peaks.clust: vectors S and dM have to have the same length")
   if (length(BinSize)!=2) BinSize=c(0,1)
   if (BinSize[1]>BinSize[2]) { a=BinSize[1]; BinSize[1]=BinSize[2]; BinSize[2]=a } 
   if (verbose) 
     print("Stack #Gaps [parent_bin ](bin_size) -> [Left_child ](#reps bin_size) + [right_child](#reps bin_size)  gap_chosen") 
   mStack = max(20, as.integer(3*log(nS))) # initial stack size   
   nStack = 1                   # stack will be seeded with a single record
   Stack = matrix(0, 2, mStack) # "Stack" used to store bins that still have to be divided
   Stack[1,nStack] = 1          # field 1 stores left boundary of the bin (index number in mass array)
   Stack[2,nStack] = nS         # field 2 stores right boundary of the bin (index number in mass array)
   bin = numeric(nS)            # bin is binary array where bin (cluster) left boundaries are stored
   while (nStack>0) {           # if there are any repeats within cluster than try to remove them
     from= Stack[1,nStack]      # get first item from the queue
     to  = Stack[2,nStack]      # [from:to] is a range of peaks that needs splitting
     nStack = nStack - 1        # decrease stack size
     
     #==================================
     # Find biggest gap to divide the bin
     #==================================
     gap = dM[from:(to-1)]             # gap sizes within range in question
     size= sum(gap)                    # bin size
     mx  = tol*max(gap)                # find largest gap and gaps of the "similar" size
     idx = which(gap>=mx)              # find biggest gap between peaks first
     len = length(idx)                 # number of equal size gaps
     if (len>1) {                      # if there are multiple equal size gaps
       score = numeric(len)            # then find one that reduces number of repeats the most
       for (j in 1:len) {              # test each gap in decending order until we find a split that will reduce number of repeats
         Cut = idx[j]+from             # proposed split point
         nL = repeats(S[from:(Cut-1)]) # how many repeats in the group to the left
         nR = repeats(S[Cut:to])       # how many repeats in the group to the right
         score[j] = nL+nR              # number of repeats will decrease if we split
       }
       i = which(score==min(score))    # find split that reduces number of repeats the most
       idx = idx[i] 
       m   = length(idx) 
       if (m>1) {                      # if there are multiple ones that perform the same way
         k = to-from-1 
         if (idx[1]<k-idx[m]) idx=idx[1] # try to keep one of the bins the maximum size
         else                 idx=idx[m]
       }
     }
     Cut = idx[1]+from
     
     #==================================
     # Left subspace processing
     #==================================
     if (from<=Cut-2) sL=sum(dM[from:(Cut-2)]) else sL=0 # size of the bin to the left
     nL=-1                           # special number for verbose mode means "bin too small" rule was used
     if (sL>BinSize[2]) {            # left subspace is too big - add it to the stack
       nStack = nStack + 1
       Stack[,nStack] = c(from, Cut-1)
       nL=-2                         # special number for verbose mode means "bin too big" rule was invoked
     } else if (sL>BinSize[1]) {     # left subspace could be too small and will not be splited 
       nL = repeats(S[from:(Cut-1)]) # how many repeats in the group to the left
       if (nL>0) {                   # left subspace still have repeats - add it to the stack
         nStack = nStack + 1
         Stack[,nStack] = c(from, Cut-1)
       }
     } else if (sL==0) nL=0
     
     #==================================
     # Right subspace processing
     #==================================
     if(Cut<=to-1) sR=sum(dM[Cut:(to-1)]) else sR=0 # size of the bin to the right
     nR=-1                           # special number for verbose mode means "bin too small" rule was used
     if (sR>BinSize[2]) {            # right subspace is too big - add it to the stack to be splited
       nStack = nStack + 1
       Stack[,nStack] = c(Cut, to)
       nR=-2                         # special number for verbose mode means "bin too big" rule was invoked
     } else if (sR>BinSize[1]) {     # right subspace could be too small and will not be splited 
       nR = repeats(S[Cut:to])       # how many repeats in the group to the right
       if (nR>0) {                   # right subspace still have repeats - add it to the stack
         nStack = nStack + 1
         Stack[,nStack] = c(Cut, to)
       }
     } else if (sR==0) nR=0
     
     if (nStack==mStack) {           # incease stack size if needed
       mStack = trunc(1.5*mStack)    # hopefully rarely used
       Stack  = rbind(Stack, matrix(0,2, mStack-nStack))
     }
     bin[Cut]=1
     if (verbose) 
       print(sprintf("%5i %5i [%5i %5i](%8.5f) -> [%5i %5i](%5i, %7.4f) + [%5i %5i](%5i, %7.4f)  gap=%6.4f", 
        as.integer(nStack), as.integer(len), 
        as.integer(from), as.integer( to  ), size,
        as.integer(from), as.integer(Cut-1), as.integer(nL),  sL, 
        as.integer(Cut ), as.integer( to  ), as.integer(nR),  sR, dM[Cut-1]))
   } 
   bin[1 ] = 1
   bin[nS] = 0
   return(bin)
}

