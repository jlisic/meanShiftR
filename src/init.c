/* Copyright (c) 2015-2016  Jonathan Lisic 
 * Last edit: 16/08/24 - 09:20:44
 * License: GPL (>=2) 
 */  

#include <stdio.h> 
#include <stdlib.h>
//#include <time.h>

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

static R_NativePrimitiveArgType R_knn_t[] = {
      REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, REALSXP, INTSXP, INTSXP, 
      REALSXP, INTSXP, REALSXP
};

static const R_CMethodDef cMethods[] = {
     {"R_meanShift", (DL_FUNC) &R_meanShift, 14, R_meanShift_t},
     {"R_knn", (DL_FUNC) &R_knn, 14, R_knn_t},
        {NULL, NULL, 0, NULL}
};

void R_init_myLib(DllInfo *info)
{
     R_registerRoutines(info, cMethods, NULL, NULL, NULL);
     R_useDynamicSymbols(info, TRUE); 
}


