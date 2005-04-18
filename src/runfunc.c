#include <stdlib.h>
#include <math.h>

/*============================================================================================*/

void insertion_sort(const double *V, int *idx, const int nIdx) 
{ 
  /* Each iteration of an insertion sort removes an element from the input data, inserting  */
  /*  it at the correct position in the already sorted list, until no elements are left     */
  /* in the input. For unsorted data is much less efficient than the more advanced          */
  /* algorithms such as Quicksort, Heapsort, or Merge sort, but it is very efficient on     */
  /* data sets which are already substantially sorted (almost O(n))                         */
  /* Referances:                                                                            */
  /*   R. Sedgewick: Algorithms. Addison-Wesley (1988) (page 99)                            */
  /*   http://en.wikipedia.org/wiki/Insertion_sort                                          */
  /*   http://www.cs.ubc.ca/spider/harrison/Java/sorting-demo.html                          */
  /* Input:                                                                                 */
  /*   V    - data array we operate on remains unchanged (we assume that no number in idx   */
  /*          array is longer than V                                                        */
  /*   idx  - index numbers of elements in V to be partially sorted                         */
  /*   nIdx - length of idx array                                                           */
  /* Output:                                                                                */
  /*   idx  - index numbers of sorted elements in V w                                       */
  int i, j, id;
  double v;
  for (i=1; i<nIdx; i++) {
    id = idx[i];
    v  = V[id];
    for(j=i; j>0; j--) {
      if (V[idx[j-1]]<v) break;
      idx[j] = idx[j-1];
    }
    idx[j] = id;
  }
}

/*============================================================================================*/

void runmean(double *In, double *Out, const int *nIn, const int *nWin)
/* Mean function applied to (running) window                      */  
/* Input :                                                        */
/*   In   - array to run moving window over will remain umchanged */
/*   Out  - empty space for array to store the results            */
/*   nIn  - size of arrays In and Out                             */
/*   nWin - size of the moving window                             */
/* Output :                                                       */
/*   Out  - results of runing moving window over array In and colecting window mean */
{
  int i, k, n=*nIn, m=*nWin;
  double *in, *out;
  k = m/2;                              /* half of moving window */                           
  for(i=0; i<=k; i++) Out[i]=Out[n-i-1]=0;
  if (m<n) {
    in=In; out=Out+k; 
    for(i=0; i<m; i++) *out += *(in++); /* sum of initial window out[k]=sum(in[0:m-1]) */
    for(out++; i<n; i++, out++, in++)   /* "in" already set: in=In+m  */
      *out = *(out-1) + *in - *(in-m); 
    out=Out+k;    
    for(i=m; i<=n; i++) *(out++) /= m;  /* out = out/m       */
  }
}

/*============================================================================================*/

void runquantile(double *In, double *Out, const int *nIn, const int *nWin, const double *Prob, const int *nProb)
/* quantile function applied to (running) window                  */ 
/* Input :                                                        */
/*   In   - array to run moving window over will remain umchanged */
/*   Out  - empty space for array to store the results            */
/*          Out is assumed to have reserved memory for nIn*nProbs elements */
/*   nIn  - size of arrays In and Out                             */
/*   nWin - size of the moving window (odd)                       */
/*   Prob - in the moving window with sorted elements X[1]:X[nWin] which element to return? */
/*          1 is minimum, nWin is maximum, (nWin+1)/2 is median, non-integers cause interpolation */
/*          multiple probs are allowed.                           */
/*   nProb - How many elements in Probs array?                    */
/* Output :                                                       */
/*   Out  - results of runing moving window over array In and colecting window mean */
{
  int i, j, k, *idx, d, n=*nIn, m=*nWin, nPrb=*nProb;
  double *Win, *in, *out, pointOut, ext, r, ip;
  k   = m/2;                          /* half of window size */
  in  = In;
  out = Out+k;

  if (nPrb==1 && (*Prob==1 || *Prob==m)) { /* trivial case shortcut - if prob is 0 or 1 than wind windows min or max */
    d = (*Prob==1 ? -1 : 1);          /* runmin d=-1; runmax d=1*/
    pointOut=ext=0;
    for(i=m-1; i<n; i++) {
      if(pointOut==ext) {             /* if point comining out of the window was window's extreme than ... */
        ext=in[0];                    /* in the 2 lines below we do things diffferent when extreme is a max vs. min */
        if (d==1) { 
          for(j=1; j<m; j++) if (ext<in[j]) ext=in[j];  /* find maximum over a window of length m */
        } else {   
          for(j=1; j<m; j++) if (ext>in[j]) ext=in[j];  /* find minimum over a window of length m */
        } 
      } else                          /* if point comining out of the window was NOT window extreme than we know ... */
        if (ext*d<in[m-1]*d) ext=in[m-1]; /* ... extreme of window's first m-1 points, so we have to add a single point */
      pointOut = *(in++);             /* store point comming out of the window for future use and move window */
      *(out++) = ext;                 /* and fill the space with window extreme and move window */
    }
  } else {                            /* non-trivial case */
    idx = (int   *) malloc(m*sizeof(int   )); /* index will hold partially sorted index numbers of Save array */
    Win = (double*) malloc(m*sizeof(double)); /* stores all points of the current running window */
    for(i=0; i<m; i++) {
      Win[i] = *(in++);               /* initialize running window */
      idx[i] = i;                     /* and its index */
    }
    in--;                             /* last point of the window will be placed again */
    for(j=i=m-1; i<n; i++) {
      Win[j] = *(in++);               /* Move Win to the right: replace a[i-m] with a[m] point  */
      insertion_sort(Win,idx,m);      /* sort current Win */
      for(d=0; d<nPrb; d++) {
        r = modf( Prob[d], &ip );     /* Divide p into its fractional and integer parts            */
        k = (int) ip-1;               /* k-1 instead of k because in C arrays are 0 based and in R they are 1 based */
        if (r) r = Win[idx[k]]*(1-r) + Win[idx[k+1]]*r; /* interpolate */
        else   r = Win[idx[k]];  
        *(out+d*n) = r;
      }
      out++;
      j = (j+1)%m;                    /* index goes from 0 to m-1, and back to 0 again  */
    }
    free(Win);
    free(idx);
  }
}

/*============================================================================================*/

void runmad(double *In, double *Ctr, double *Out, const int *nIn, const int *nWin)
/* MAD function applied to moving (running) window                */ 
/* Input :                                                        */
/*   In   - array to run moving window over will remain umchanged */
/*   Ctr  - array storing results of runmed                       */
/*   Out  - empty space for array to store the results            */
/*   nIn  - size of arrays In and Out                             */
/*   nWin - size of the moving window                             */
/* Output :                                                       */
/*   Out  - results of runing moving window over array In and colecting window mean */
{ 
  int i, k, j, l, *idx2, n=*nIn, m=*nWin;
  double *Win1, *Win2, *in, *out, *ctr, med0, med;

  idx2 = (int   *) malloc(m*sizeof(int   )); /* index will hold partially sorted index numbers of Save array */
  Win1 = (double*) malloc(m*sizeof(double)); /* stores all points of the current running window */
  Win2 = (double*) malloc(m*sizeof(double)); /* stores all points of the current running window */
  k    =  m/2;                   /* half of moving window size */  
  in   = In;                     /* initialize pointer to input In vector */
  out  = Out+k;                  /* initialize pointer to output Mad vector */
  ctr  = Ctr+k;                  /* initialize pointer to output Mad vector */
  med0 = 0;                      /* med0 - will save previous center (median) sowe know it changed */
  for(i=0; i<m; i++) {
    Win1[i] = Win2[i] = *(in++); /* initialize running windows */
    idx2[i] = i;                 /* and its indexes */
  }
  in--;                          /* last point of the window will be placed again */
  for(j=i=m-1; i<n; i++) {       /* the main loop */
    Win1[j] = *(in++);           /* Move Win to the right: replace a[i-m] with a[m] point  */    
    med     = *(ctr++);          /* find median of current Win1 */
    if (med==med0)               /* median did not changed */
      Win2[j]=fabs(Win1[j]-med); /* recalculate one point */
    else for(l=0; l<m; l++)      /* median did not changed */
      Win2[l]=fabs(Win1[l]-med); /* recalculate whole window Win2 */
    insertion_sort(Win2,idx2,m); /* sort Win2 - if medians did not change than  data should be sorted*/
    *(out++) = Win2[idx2[k]];    /* find mad of current Win1 and store it */
    med0 = med;                  /* save previous median */
    j = (j+1)%m;                 /* index goes from 0 to m-1, and back to 0 again  */
  }
  free(Win2);
  free(Win1);
  free(idx2);
}

