/* 
* License: GPL (>=2) 
 */  

#ifndef DFTOSPARSE_HEADER

#define DFTOSPARSE_HEADER

#include <stdio.h> 
#include <stdlib.h>
//#include <time.h>

#include "R.h"
#include "Rinternals.h"
#include "Rmath.h"

#if defined _OPENMP
  #include <omp.h>
#endif



/***********************************************************************************************/
/* This is a function to turn an R data frame into a sparse matrix                             */
/***********************************************************************************************/

void R_dfToSparse(
 int * x ,     // 1 Not used yet
 int *n,          // 2 number of rows
 int *p,          // 3 number of columns
 int *new_i,      // 4 column indexes
 int *new_p,      // 5 pointers for sparse matrix
 int *max_cat,    // 6 number of categories per column
 int *max_cat_cumsum , // 7 cumulative sume of the number of categories per column
 int *type              // 8 Not used yet
 );


#endif

