/* Copyright (c) 2015-2016  Jonathan Lisic 
 * Last edit: 16/08/24 - 09:20:44
 * License: GPL (>=2) 
 */  

#ifndef MEANSHIFT_HEADER

#define MEANSHIFT_HEADER

#include <stdio.h> 
#include <stdlib.h>
//#include <time.h>

#include "R.h"
#include "Rinternals.h"
#include "Rmath.h"
#include "kdtree.h"

#if defined _OPENMP
  #include <omp.h>
#endif

/***********************************************************************************************/
/* print matrix                                                                                */
/* NOTE: printing is done in column major form                                                 */
/***********************************************************************************************/

/* function that prints out a matrix of type integer */
void printMatrixFullInt(int * X , size_t row, size_t col ); 


/* function that prints out a matrix of type double */
void printMatrixFullDbl(double * X , size_t row, size_t col ) ;


/* function that prints out a matrix of type integer */
void printMatrixFullSize_t(size_t * X , size_t row, size_t col );


/***********************************************************************************************/
/* this is an alternative kernel, with pre computed distances                                  */
/***********************************************************************************************/
/*
 * query      - query matrix    queryRow by queryCol
 * train      - train matrix    trainRow by queryCol
 * neighbors  - neighbors       nNeighbors by queryRow
 * queryRow   - rows for query
 * trainRow   - rows for train
 * queryCol   - columns for query and train
 * nNeighbors - number of neighbors
 * bandwidth  - bandwidth parameter
 * alpha      - alpha parameter
 */
void normalKernelNewton_preDist(
  double * query,        /* query matrix */
  double * train,        /* train matrix */
  size_t * neighbors,      /* neighbors matrix */
  size_t queryRow,
  size_t trainRow,
  size_t queryCol,
  size_t nNeighbors,
  double * bandwidth2,    /* the squared reciprical of the bandwidth */
  double * distance,
  double alpha
  );



/***********************************************************************************************/
/* this is a normal kernel                                                               */
/***********************************************************************************************/
/*
 * query      - query matrix    queryRow by queryCol
 * train      - train matrix    trainRow by queryCol
 * neighbors  - neighbors       nNeighbors by queryRow
 * queryRow   - rows for query
 * trainRow   - rows for train
 * queryCol   - columns for query and train
 * nNeighbors - number of neighbors
 * bandwidth  - bandwidth parameter
 * alpha      - alpha parameter
 */
void normalKernelNewton(
  double * query,        /* query matrix */
  double * train,        /* train matrix */
  size_t * neighbors,      /* neighbors matrix */
  size_t queryRow,
  size_t trainRow,
  size_t queryCol,
  size_t nNeighbors,
  double * bandwidth2,    /* the squared reciprical of the bandwidth */
  double alpha
  );


/***********************************************************************************************/
/* this is a Epanechnikov product kernel                                                       */
/***********************************************************************************************/
/*
 * query      - query matrix    queryRow by queryCol
 * train      - train matrix    trainRow by queryCol
 * neighbors  - neighbors       nNeighbors by queryRow
 * queryRow   - rows for query
 * trainRow   - rows for train
 * queryCol   - columns for query and train
 * nNeighbors - number of neighbors
 * bandwidth  - bandwidth parameter
 * alpha      - alpha parameter
 */
void epanechnikov_Kernel(
  double * query,        /* query matrix */
  double * train,        /* train matrix */
  size_t * neighbors,      /* neighbors matrix */
  size_t queryRow,
  size_t trainRow,
  size_t queryCol,
  size_t nNeighbors,
  double * bandwidth    /* the squared reciprical of the bandwidth */
  );



/***********************************************************************************************/
/* this is a biweight product kernel                                                            */
/***********************************************************************************************/
/*
 * query      - query matrix    queryRow by queryCol
 * train      - train matrix    trainRow by queryCol
 * neighbors  - neighbors       nNeighbors by queryRow
 * queryRow   - rows for query
 * trainRow   - rows for train
 * queryCol   - columns for query and train
 * nNeighbors - number of neighbors
 * bandwidth  - bandwidth parameter
 * alpha      - alpha parameter
 */
void biweight_Kernel(
  double * query,        /* query matrix */
  double * train,        /* train matrix */
  size_t * neighbors,      /* neighbors matrix */
  size_t queryRow,
  size_t trainRow,
  size_t queryCol,
  size_t nNeighbors,
  double * bandwidth    /* the squared reciprical of the bandwidth */
  );


/* mean shfit nearest neighbors */
void R_meanShift(
  double * query,                /* data to query for, in row major form */
  double * train,                /* reference data, in row major form */
  int * assignment,              /* assignment */
  int * queryRowPtr,             /* number of rows of query data */
  int * trainRowPtr,             /* number of rows of reference data */
  int * queryColPtr,             /* number of columns of query and reference*/
  int * nNeighborsPtr,           /* number of Neighbors */
  int * iterationsPtr,           /* number of iterations */
  double * bandwidth,    
  double * alphaPtr,             /* ms to newton ajustment parameter */
  double * epsilonPtr,           /* not used */
  double * epsilonClusterPtr,           /* not used */
  int * kernelEnumPtr,           /* not used */
  int * algorithmEnumPtr,        /* not used */
  int * intParameters,
  int * dblParameters
  );


#endif
