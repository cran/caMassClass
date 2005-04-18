\name{runfunc}
\alias{runmean}
\alias{runmin}
\alias{runmax}
\alias{runmad}
\alias{runquantile}
\alias{EndRule}
\title{Moving Window Analysis of a Vector}
\description{A collection of functions to perform fast moving window (running, 
  rolling window) analysis of vectors. }
\usage{
  runmean(x, k, endrule=c("NA", "trim", "keep", "constant", "func"))
  runmin (x, k, endrule=c("NA", "trim", "keep", "constant", "func"))
  runmax (x, k, endrule=c("NA", "trim", "keep", "constant", "func"))
  runmad (x, k, center=runmed(x,k,endrule="keep"), constant=1.4826,  
          endrule=c("NA", "trim", "keep", "constant", "func"))
  runquantile(x, k, probs, type=7, 
          endrule=c("NA", "trim", "keep", "constant", "func"))
  EndRule(x, y, k, 
          endrule=c("NA", "trim", "keep", "constant", "func"), Func, \dots)
}

\arguments{
  \item{x}{numeric vector of length n}
  \item{k}{width of moving window; must be an odd integer bigger than one. }
  \item{endrule}{character string indicating how the values at the beginning 
    and the end, of the data, should be treated. Only first and last \code{k2} 
    values at both ends are affected, where \code{k2} is the half-bandwidth 
    \code{k2 = k \%/\% 2}.
     \itemize{
       \item \code{"trim"} - trim the ends output array length is equal to 
         \code{length(x)-2*k2 (out = out[(k2+1):(n-k2)])}. This option mimics 
         output of \code{\link{apply}} \code{(\link{embed}(x,k),1,FUN)} and other 
         related functions.
       \item \code{"keep"} - fill the ends with numbers from \code{x} vector 
         \code{(out[1:k2] = x[1:k2])}
       \item \code{"constant"} - fill the ends with first and last calculated 
         value in output array \code{(out[1:k2] = out[k2+1])}
       \item \code{"NA"} - fill the ends with NA's \code{(out[1:k2] = NA)}
       \item \code{"func"} - applies the underlying function to smaller and 
       smaller sections of the array. For example in case of mean: 
       \code{for(i in 1:k2) out[i]=mean(x[1:i])}. This option is not optimized 
       and could be very slow for large windows.
     }
     Similar to \code{endrule} in \code{\link{runmed}} function which has the 
     following options: \dQuote{\code{c("median", "keep", "constant")}} .
  }
  \item{center}{moving window center used by \code{runmad} function defaults 
    to running median (\code{\link{runmed}} function). Similar to \code{center}  
    in \code{\link{mad}} function. }
  \item{constant}{scale factor used by \code{runmad}, such that for gaussian 
    distribution X, \code{\link{mad}}(X) is the same as \code{\link{sd}}(X). 
    Same as \code{constant} in \code{\link{mad}} function.}
  \item{probs}{numeric vector of probabilities with values in [0,1] range 
    used by \code{runquantile}. For example \code{Probs=c(0,0.5,1)} would be 
    equivalent to running \code{runmin}, \code{\link{runmed}} and \code{runmax}.
    Same as \code{probs} in \code{\link{quantile}} function. }
  \item{type}{an integer between 1 and 9 selecting one of the nine quantile 
    algorithms, same as \code{type} in \code{\link{quantile}} function. 
    Another even more readable description of nine ways to calculate quantiles 
    can be found at \url{http://mathworld.wolfram.com/Quantile.html}. }
  \item{y}{numeric vector of length n, which is partially filled output of 
    one of the \code{run} functions. Function \code{EndRule} will fill the 
    remaining begining and end sections using method chosen by \code{endrule} 
    argument.}
  \item{Func}{Function name that \code{EndRule} will use in case of 
    \code{endrule="func"}.}
  \item{\dots}{Additionam parameters that \code{EndRule} will use in case of 
    \code{endrule="func"}.}
}

\details{
  Apart from the end values, the result of y = runFUN(x, k) is the same as 
  \dQuote{\code{for(j=(1+k2):(n-k2)) y[j]=FUN(x[(j-k2):(j+k2)])}}, where FUN 
  stands for min, max, mean, mad or quantile functions.

  The main incentive to write this set of functions was relative slowness of 
  majority of moving window functions available in R and its packages.  With 
  exception of \code{\link{runmed}}, a running window median function, all 
  functions listed in "see also" section are slower than very inefficient 
  \dQuote{\code{\link{apply}(\link{embed}(x,k),1,FUN)}} approach. Speeds of 
  above functions are as follow:
     \itemize{
       \item \code{runmin}, \code{runmax}, \code{runmean} run at O(n)
       \item \code{runquantile} and \code{runmad} run at O(n*k)
       \item \code{\link{runmed}} - related R function run at O(n*log(k)) 
      }
  
  Functions \code{runquantile} and \code{runmad} are using insertion sort to 
  sort the moving window, but gain speed by remembering results of the previous 
  sort. Since each time the window is moved, only one point changes, all but one 
  points in the window are already sorted. Insertion sort can fix that in O(k) 
  time.
  
  Function \code{runquantile} when run in single probability mode automatically
  recognizes probabilities: 0, 1/2, and 1 as special cases and return output 
  from functions: \code{runmin}, \code{\link{runmed}} and \code{runmax} 
  respectivly. 

  All \code{run*} functions are written in C, but \code{runmin}, \code{runmax} 
  and \code{runmean} also have R code versions that can be used if DLL is not 
  loaded. 
  
  Function \code{EndRule} applies one of the five methods (see \code{endrule} 
  argument) to process end-points of the input array \code{x}. 
}


\value{
  Functions \code{runmin}, \code{runmax}, \code{runmean} and \code{runmad} 
  return a numeric vector of the same length as \code{x}. 
  Function \code{runquantile} returns a matrix of size [n \eqn{\times}{x} 
  length(probs)].
} 

\references{
  Hyndman, R. J. and Fan, Y. (1996) Sample quantiles in statistical packages, 
  American Statistician, 50, 361. 
} 

\author{Jarek Tuszynski (SAIC) \email{jaroslaw.w.tuszynski@saic.com}} 

\seealso{
  Links related to each function:
  \itemize{       
   \item \code{runmin} - \code{\link{min}}, \code{\link[fSeries]{rollMin}} from 
     \pkg{fSeries} library
   \item \code{runmax} - \code{\link{max}}, \code{\link[fSeries]{rollMax}} from 
     \pkg{fSeries} library
   \item \code{runmean} - \code{\link{mean}}, \code{\link{kernapply}}, 
     \code{\link{filter}}, 
     \code{\link[fSeries]{rollMean}} from \pkg{fSeries} library, 
     \code{\link[magic]{subsums}} from \pkg{magic} library
   \item \code{runquantile} - \code{\link{quantile}}, \code{\link{runmed}}, 
     \code{\link{smooth}}
   \item \code{runmad} - \code{\link{mad}}, \code{\link[fSeries]{rollVar}} from 
     \pkg{fSeries} library
   \item generic running window functions: \code{\link{apply}
     (\link{embed}(x,k), 1, FUN)} (fastest), \code{\link[fSeries]{rollFun}} 
     from \pkg{fSeries} (slow), \code{\link[gtools]{running}} from \pkg{gtools} 
     package (extrimly slow for this purpose), \code{\link[magic]{subsums}} from 
     \pkg{magic} library can perform running window operations on data with any 
     dimensions.
   \item \code{EndRule} - \code{\link{smoothEnds}(y,k)} function is similar to 
     \code{EndRule(x,y,k,endrule="func", median)}
  }
}

\examples{
  # test runmin, runmax and runmed
  k=15; n=200;
  x = rnorm(n,sd=30) + abs(seq(n)-n/4)
  col = c("black", "red", "green", "blue", "magenta", "cyan")
  plot(x, col=col[1], main = "Moving Window Analysis Functions")
  lines(runmin(x,k), col=col[2])
  lines(runmed(x,k), col=col[3])
  lines(runmax(x,k), col=col[4])
  legend(0,.9*n, c("data", "runmin", "runmed", "runmax"), col=col, lty=1 )

  #test runmean and runquantile
  y=runquantile(x, k, probs=c(0, 0.5, 1, 0.25, 0.75), endrule="constant")
  plot(x, col=col[1], main = "Moving Window Quantile")
  lines(runmean(y[,1],k), col=col[2])
  lines(y[,2], col=col[3])
  lines(runmean(y[,3],k), col=col[4])
  lines(y[,4], col=col[5])
  lines(y[,5], col=col[6])
  lab = c("data", "runmean(runquantile(0))", "runquantile(0.5)", 
  "runmean(runquantile(1))", "runquantile(.25)", "runquantile(.75)")
  legend(0,0.9*n, lab, col=col, lty=1 )

  #test runmean and runquantile
  k =25
  m=runmed(x, k)
  y=runmad(x, k, center=m)
  plot(x, col=col[1], main = "Moving Window Analysis Functions")
  lines(m    , col=col[2])
  lines(m-y/2, col=col[3])
  lines(m+y/2, col=col[4])
  lab = c("data", "runmed", "runmed-runmad/2", "runmed+runmad/2")
  legend(0,1.8*n, lab, col=col, lty=1 )

  # speed comparison
  x=runif(100000); k=991;
  system.time(runmean(x,k))
  system.time(filter(x, rep(1/k,k), sides=2)) #the fastest alternative
  k=91;
  system.time(runmad(x,k))
  system.time(apply(embed(x,k), 1, mad)) #the fastest alternative
}

\keyword{ts}
\keyword{smooth}
\keyword{array}
\keyword{utilities}
