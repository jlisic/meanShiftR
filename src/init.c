/* Copyright (c) 2015-2016  Jonathan Lisic 
 * Last edit: 16/08/24 - 09:20:44
 * License: GPL (>=2) 
 */  

#include <stdio.h> 
#include <stdlib.h>

#include "R.h"
#include "Rinternals.h"
#include "Rmath.h"
#include <R_ext/Rdynload.h>
#include "kdtree.h"
#include "meanShift.h"
#if defined _OPENMP
  #include <omp.h>
#endif

/***********************************/
/* Register SO's                   */

static R_NativePrimitiveArgType R_meanShift_t[] = {
      REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, 
      REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP
};

static R_NativePrimitiveArgType R_knn_sparse_t[] = {
  REALSXP, // 1 query points
  INTSXP,  // 2 query index 
  INTSXP,  // 3 query pointers 
  INTSXP,  // 4 query index n 
  INTSXP,  // 5 query pointers n 
  REALSXP, // 6 x
  INTSXP,  // 7 x index 
  INTSXP,  // 8 x pointers 
  INTSXP,  // 9 x index n 
  INTSXP,  // 10 x pointers n 
  INTSXP,  // 11 xnrowPtr
  INTSXP,  // 12 nrowPtr
  INTSXP,  // 13 ncolPtr
  REALSXP, // 14 kDist
  INTSXP,  // 15 indexInt
  INTSXP,  // 16 kPtr
  REALSXP, // 17 weight
  INTSXP,  // 18 leafSizePtr
  REALSXP, // 19 maxDist
  INTSXP   // 20 sparse
};

static R_NativePrimitiveArgType R_knn_t[] = {
  REALSXP, // 1 query points
  REALSXP, // 2 x
  INTSXP,  // 3 xnrowPtr
  INTSXP,  // 4 nrowPtr
  INTSXP,  // 5 ncolPtr
  REALSXP, // 6 kDist
  INTSXP,  // 7 indexInt
  INTSXP,  // 8 kPtr
  REALSXP, // 9 weight
  INTSXP,  // 10 leafSizePtr
  REALSXP  // 11 maxDist
};

static const R_CMethodDef cMethods[] = {
     {"R_meanShift", (DL_FUNC) &R_meanShift, 14, R_meanShift_t},
     {"R_knn", (DL_FUNC) &R_knn, 11, R_knn_t},
     {"R_knn_sparse", (DL_FUNC) &R_knn_sparse, 20, R_knn_sparse_t},
     {NULL, NULL, 0, NULL}
};

void R_init_myLib(DllInfo *info)
{
     R_registerRoutines(info, cMethods, NULL, NULL, NULL);
     R_useDynamicSymbols(info, TRUE); 
}


