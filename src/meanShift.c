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
#include "kdtree.h"
#include "meanShift.h"

#if defined _OPENMP
  #include <omp.h>
#endif



/***********************************************************************************************/
/* print matrix                                                                                */
/* NOTE: printing is done in column major form                                                 */
/***********************************************************************************************/

/* function that prints out a matrix of type integer */
void printMatrixFullInt(int * X , size_t row, size_t col ) {
  size_t i,j;

  for(i = 0; i < row; i++) {
    Rprintf("%d:\t",(int) i);
    for(j = 0; j < col; j++) {
      Rprintf("%d\t",X[i*col + j]);
    }
    Rprintf("\n");
  }

}


/* function that prints out a matrix of type double */
void printMatrixFullDbl(double * X , size_t row, size_t col ) {
  size_t i,j;

  for(i = 0; i < row; i++) {
    Rprintf("%d:\t",(int) i);
    for(j = 0; j < col; j++) {
      Rprintf("%0.4f\t",X[i*col + j]);
    }
    Rprintf("\n");
  }

}


/* function that prints out a matrix of type integer */
void printMatrixFullSize_t(size_t * X , size_t row, size_t col ) {
  size_t i,j;

  for(i = 0; i < row; i++) {
    Rprintf("%d:\t",(int) i);
    for(j = 0; j < col; j++) {
      Rprintf("%d\t",(int) X[i*col + j]);
    }
    Rprintf("\n");
  }

}


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
  ) {

  size_t i; /* query data index */
  size_t k; /* neighbor index */
  size_t d; /* column index */


  double * numerator;
  double w;
  double * denominator;
  double * queryValue;
  double * neighborValue;
  double tmp;


//  omp_set_num_threads(1);

  /* iterate over each query point */ 
  #pragma omp parallel private(i,d,k,neighborValue, queryValue,w,numerator,denominator)
  {
 
  // alloc from heap 
  numerator   = (double * ) calloc( queryCol, sizeof(double)); // used to store numerator
  denominator = (double * ) calloc( queryCol, sizeof(double)); // used to store denominator
  

  #pragma omp for
  for( i = 0; i < queryRow; i++)  {

    // check if we are already at the terminal condition 
    if( neighbors[i*nNeighbors] == nNeighbors) continue;


    /* focal weight to divide by for the denominator */
    for( d = 0; d < queryCol; d++) { 
      numerator[d] = 0; 
      denominator[d] = 0; 
    }
      
    queryValue = query + i* queryCol;

    /* iterate over each neighbor */
    for( k = 0; k < nNeighbors; k++)  { 
      
      if( distance[ i*nNeighbors + k] == INFINITY ) continue;
     
      neighborValue = train + neighbors[i*nNeighbors + k] * queryCol;

      // for each point query over each dimension 
      // double dnorm( double x, double mu, double sigma, int give_log)   
      //w *=  dnorm( z[d] , 0, 1 , 0 );
      w =  exp( - 0.5 * distance[ i*nNeighbors + k] );
      
      /* aggregate over all variables */
      for( d = 0; d < queryCol; d++) {
        tmp =  queryValue[d] - neighborValue[d] ;
        numerator[d] += tmp * w;
        if( alpha > 0 ) {
          denominator[d] += ( 1 - alpha * tmp * tmp) * w ;
        } else { 
          denominator[d] += w;
        }
      }

    }

    for( d = 0; d < queryCol; d++)  
      if( denominator[d] != 0 )  queryValue[d] = queryValue[d] - numerator[d] / denominator[d];
     
  }
  
  free(denominator); 
  free(numerator);

  } // for omp

  return;
}



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
  ) {

  size_t i; /* query data index */
  size_t k; /* neighbor index */
  size_t d; /* column index */


  double * numerator;
  double w;
  double * denominator;
  double * z;
  size_t currentNeighbor;


//  omp_set_num_threads(1);

  /* iterate over each query point */ 
  #pragma omp parallel private(i,d,k,z,w,numerator,denominator,currentNeighbor)
  {
 
  // alloc from heap 
  z           = (double * ) calloc( queryCol, sizeof(double)); // used to store difference
  numerator   = (double * ) calloc( queryCol, sizeof(double)); // used to store numerator
  denominator = (double * ) calloc( queryCol, sizeof(double)); // used to store denominator
  

  #pragma omp for
  for( i = 0; i < queryRow; i++)  {

    // check if we are already at the terminal condition 
    if( neighbors[i*nNeighbors] == nNeighbors) continue;


    /* focal weight to divide by for the denominator */
    for( d = 0; d < queryCol; d++) { 
      numerator[d] = 0; 
      denominator[d] = 0; 
    }

    /* iterate over each neighbor */
    for( k = 0; k < nNeighbors; k++)  { 
      
      currentNeighbor = neighbors[i*nNeighbors + k];

      /* set weight to 1 */ 
      w = 1;

      // for each point query over each dimension 
      for( d = 0; d < queryCol; d++)  { 
        z[d] = ( query[i * queryCol + d] - train[currentNeighbor * queryCol + d] );
        // double dnorm( double x, double mu, double sigma, int give_log)   
        //w *=  dnorm( z[d] , 0, 1 , 0 );
        w *= exp( - 0.5 * z[d]*z[d] * bandwidth2[d] );
      }

      /* aggregate over all variables */
      for( d = 0; d < queryCol; d++) {
        numerator[d] +=  z[d] * w ;
        if( alpha > 0 ) {
          denominator[d] += ( 1 - alpha * z[d] * z[d]) * w ;
        } else { 
          denominator[d] += w;
        }
      }

    }

    for( d = 0; d < queryCol; d++)  
      if( denominator[d] != 0 )  query[i* queryCol + d] = query[i*queryCol +d] - numerator[d] / denominator[d];
     
  }
  
  free(z);
  free(denominator); 
  free(numerator);

  } // for omp

  return;
}


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
  ) {

  size_t i; /* query data index */
  size_t k; /* neighbor index */
  size_t d; /* column index */

  double * numerator;
  double * denominator;
  double w;
  size_t currentNeighbor;
  double * u2; 
  double l2_distance;

  //omp_set_num_threads(1);

  /* iterate over each query point */ 
  #pragma omp parallel private(i,d,k,w,currentNeighbor,numerator,denominator,u2,l2_distance)
  {

  // alloc from heap 
  numerator = (double * ) calloc( queryCol, sizeof(double)); // used to store difference
  denominator = (double * ) calloc( queryCol, sizeof(double)); // used to store difference
  u2  = (double * ) calloc( queryCol, sizeof(double)); // used to store difference
  
  #pragma omp for
  for( i = 0; i < queryRow; i++)  {
   
    // check if we are already at the terminal condition 
    if( neighbors[i*nNeighbors] == nNeighbors) continue;

    for( d = 0; d < queryCol; d++) {
      numerator[d] = 0; 
      denominator[d] = 0; 
    }
   
    // iterate over each neighbor 
    for( k = 0; k < nNeighbors; k++)  { 
      
      currentNeighbor = neighbors[i*nNeighbors + k];
      w = 1;
      l2_distance = 0;

      // calculate total and l2 distance components 
      for( d = 0; d < queryCol; d++)  { 
        u2[d]=( query[i * queryCol + d] - train[currentNeighbor * queryCol + d] ) / bandwidth[d]; 
        u2[d] *= u2[d];
        w *= (1-u2[d]); 
        l2_distance += u2[d];
      }

      if( l2_distance > 1 ) w = 0;
      
      // aggregate the result to the numerator and denominator 
      if( w > 0 ) { 
        for( d = 0; d < queryCol; d++)  { 
          if( u2[d] != 1) {
            numerator[d] += train[currentNeighbor * queryCol + d] * w / (1-u2[d]); 
            denominator[d] += w / (1-u2[d]); 
          }
        }
      }  
    }

    // divide by bandwidth 
    for( d = 0; d < queryCol; d++) {
      query[i*queryCol + d] = numerator[d]/denominator[d] ;  
    }
  }
  free(u2); 
  free(numerator);
  free(denominator);

  } // for omp

  return;
}



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
  ) {

  size_t i; /* query data index */
  size_t k; /* neighbor index */
  size_t d; /* column index */

  double * numerator;
  double * denominator;
  double w;
  size_t currentNeighbor;
  double * u2;
  double  l2_distance;

  //omp_set_num_threads(1);

  /* iterate over each query point */ 
  #pragma omp parallel private(i,d,k,w,currentNeighbor,numerator,denominator,u2,l2_distance)
  {

  // alloc from heap 
  numerator = (double * ) calloc( queryCol, sizeof(double)); // used to store difference
  denominator = (double * ) calloc( queryCol, sizeof(double)); // used to store difference
  u2  = (double * ) calloc( queryCol, sizeof(double)); // used to store difference
  
  #pragma omp for
  for( i = 0; i < queryRow; i++)  {
   
    // check if we are already at the terminal condition 
    if( neighbors[i*nNeighbors] == nNeighbors) continue;

    for( d = 0; d < queryCol; d++) {
      numerator[d] = 0; 
      denominator[d] = 0;
    }
   
    // iterate over each neighbor 
    for( k = 0; k < nNeighbors; k++)  { 
      currentNeighbor = neighbors[i*nNeighbors + k];
      w = 1;
      l2_distance = 0;

      // calculate total and l2 distance components 
      for( d = 0; d < queryCol; d++)  { 
        u2[d]=( query[i * queryCol + d] - train[currentNeighbor * queryCol + d] ) / bandwidth[d]; 
        u2[d] *= u2[d];
        w *= (1-u2[d]); 
        w *= (1-u2[d]); 
        l2_distance += u2[d];
      }

      if( l2_distance > 1 ) w = 0;
      
      // aggregate the result to the numerator and denominator 
      if( w > 0 ) { 
        for( d = 0; d < queryCol; d++)  { 
          if( u2[d] != 1) {
            numerator[d] += train[currentNeighbor * queryCol + d] * w / (1-u2[d]); 
            denominator[d] += w / (1-u2[d]); 
          }
        }
      }  
      
    }

    // divide by bandwidth 
    for( d = 0; d < queryCol; d++) {
      query[i*queryCol + d] = numerator[d]/denominator[d] ;  
    }
  }
  free(u2); 
  free(numerator);
  free(denominator);

  } // for omp

  return;
}





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
  ) {

  // simplify the interface by dereferencing the pointers 
  size_t queryRow =       (size_t) * queryRowPtr; 
  size_t trainRow =       (size_t) * trainRowPtr; 
  size_t queryCol =       (size_t) * queryColPtr;
  size_t nNeighbors =     (size_t) * nNeighborsPtr;
  size_t iterations =     (size_t) * iterationsPtr;
  double alpha = *alphaPtr; 
  double epsilon = *epsilonPtr; 
  double epsilonCluster = *epsilonClusterPtr; 

  double dist_den, dist_num, dist_tmp; // variables for terminal distance calculation

  size_t i,j,k; // index 

  int m;  // max number of clusters
  double min_dist; // min dist for cluster assignment
  size_t min_j = 0;  // current min for cluster assignment

  // used to simplify bandwidth calculation
  double * bandwidth2 = (double *) calloc( queryCol, sizeof(double)); 

  // ******************************************** 
  // data structure for k-d tree 
  // ******************************************** 
  rootNodePtr kdTree = NULL;
  size_t * kdTreeIndex = NULL; // index of points for training points for kdtree
  double tieBreak = -1; // value to handle tie breaks (not used) 

  // create a data structure to store prior iterations, used to check
  // terminal conditions
  double * queryNew = (double *) calloc(queryRow*queryCol, sizeof(double));

  // create a data set for identification of nearest neighbors 
  size_t * neighbors = (size_t *) calloc(queryRow*nNeighbors, sizeof(size_t));
  double * distances = NULL;
   
  
  // ************************************************************************* 
  // Generic Setup
  // ************************************************************************* 
  for(i = 0; i < queryCol; i++) bandwidth2[i] = 1/(bandwidth[i] * bandwidth[i]);


  // ************************************************************************* 
  // Algorithm specific Setup
  // ************************************************************************* 
  // Linear 
  // *********************************** 
  if( *algorithmEnumPtr == 0 ) {
  
    #pragma omp for private(i)
    for( j=0; j < queryRow; j++) 
      for( i=0; i < nNeighbors; i++) neighbors[j*nNeighbors + i]=i;  
   
  // *********************************** 
  // K-d Tree
  // *********************************** 
  } else if( *algorithmEnumPtr == 1 ) {
      // create a place to store distances
      distances = (double *) calloc(queryRow*nNeighbors, sizeof(double));
      
      kdTree = createTree(queryCol, intParameters[0], trainRow, train);
      // this will get freed
      kdTreeIndex = calloc(trainRow, sizeof(size_t));
      #pragma omp for 
      for( i = 0; i < trainRow; i++) kdTreeIndex[i] = i;
  
      // build the tree 
      kdTree->root = buildIndex( 
        kdTree,      // root pointer 
        0,           // current dim
        trainRow,    // current length of obs
        kdTreeIndex  // pointer to obs indexes 
      ); 
      
  }  
    
  /* begin iteration */
  for( i =0; i < iterations; i++) {
    
    #pragma omp for 
    for( j = 0; j < queryRow*queryCol; j++) queryNew[j] = query[j];


    // ************************************************************************* 
    // Algorithm specific neighbor searching 
    // ************************************************************************* 
    // K-d Tree
    // *********************************** 
    if( *algorithmEnumPtr == 1 ) {

      // Setup Initial Data Sets 
      #pragma omp for 
      for(j = 0; j < queryRow*nNeighbors; j++) distances[j] = INFINITY; 

      // find neighbors
      #pragma omp for 
      for(j = 0; j < queryRow; j++)  
        if( neighbors[j*nNeighbors] !=  nNeighbors )
        find_knn( kdTree, kdTree->root, 
            nNeighbors,                // number of neighbors
            &(query[j*queryCol]),      // point to query
            &(neighbors[j*nNeighbors]), // indexes returned
            &(distances[j*nNeighbors]), // distances returned
            dblParameters[0],           // initial min distance
            bandwidth2,
            &tieBreak                  // not used
            );

    }  
 
  
  
    /**********************************************/
    /* 1. calculate kernel and update query data*/
    /**********************************************/
    /*
     * double * query      - query matrix    queryRow by queryCol  updates are done to query
     * double * train      - train matrix    trainRow by queryCol
     * size_t * neighbors  - neighbors       nNeighbors by queryRow
     * size_t queryRow   - rows for query
     * size_t trainRow   - rows for train
     * size_t  queryCol   - columns for query and train
     * size_t  nNeighbors - number of neighbors
     * double * bandwidth  - bandwidth parameter    queryCol
     * double   alpha      - alpha parameter
     */

    // ************************************************************************* 
    // Algorithm specific kernels 
    // ************************************************************************* 
    // linear 
    // *********************************** 
    
    if( *algorithmEnumPtr == 0 ) {
      if( *kernelEnumPtr == 0 ) {
        normalKernelNewton(
          query,      
          train,      
          neighbors,    
          queryRow,
          trainRow,
          queryCol,
          nNeighbors,
          bandwidth2,    
          alpha
        ); 
      }
      else if( *kernelEnumPtr == 1 ) {
        epanechnikov_Kernel(
          query,        
          train,        
          neighbors,     
          queryRow,
          trainRow,
          queryCol,
          nNeighbors,
          bandwidth    
        );
      } 
      else if( *kernelEnumPtr == 2 ) {
        biweight_Kernel(
          query,        
          train,        
          neighbors,     
          queryRow,
          trainRow,
          queryCol,
          nNeighbors,
          bandwidth    
        );
      } 
    }  
    // *********************************** 
    // K-d Tree
    // *********************************** 
    else if( *algorithmEnumPtr == 1 ) {
      normalKernelNewton_preDist(
        query,      
        train,      
        neighbors,    
        queryRow,
        trainRow,
        queryCol,
        nNeighbors,
        bandwidth2,    
        distances,
        alpha
      ); 
    }
    /**********************************************/
    /* 3. check terminal conditions */
    /**********************************************/
  
    for( j = 0; j < queryRow; j++ ) {
      dist_num = 0;
      dist_den = 0;
      for(k = 0; k < queryCol; k++ ) {
        // calculate numerator l2 dist
        dist_tmp = queryNew[j*queryCol+k] - query[j*queryCol+k]; 
        dist_tmp *= dist_tmp;
        dist_num += dist_tmp;
  
        // calculate denominator l2 dist
        dist_tmp = queryNew[j*queryCol+k];
        dist_tmp *= dist_tmp;
        dist_den += dist_tmp;
      }
  
      // set terminal condition
      // this needs to be updated with a mapping at some point in the future
      if( dist_num / dist_den < epsilon ) 
        neighbors[j*nNeighbors] = nNeighbors;
    }
    
  

  /* end iteration */
  }


  /**********************************************/
  /* 4. merge terminal points                   */ 
  /**********************************************/
  // NOTE:  there are sequential dependencies in this section


  // copy over first point
  for(k = 0; k < queryCol;k++) queryNew[k] = query[k];
  //set the first assignment
  assignment[0] = 0; 
  m = 1; 
  
  // this needs to be replaced with a tree
  for(i = 1; i < queryRow; i++ ) {
    min_dist = INFINITY;
    for( j = 0; j < m; j++ ) {
      dist_num = 0;
      dist_den = 0;
      for(k = 0; k < queryCol; k++ ) {
        // calculate numerator l2 dist
        dist_tmp = queryNew[j*queryCol+k] - query[i*queryCol+k]; 
        dist_tmp *= dist_tmp;
        dist_num += dist_tmp;
  
        // calculate denominator l2 dist
        dist_tmp = queryNew[j*queryCol+k];
        dist_tmp *= dist_tmp;
        dist_den += dist_tmp;
      }
        
  
      if( dist_num / dist_den < min_dist ) {
        min_dist = dist_num / dist_den;
        min_j = j; 
      }

    }

    // assign cluster 
    if (min_dist < epsilonCluster ) {
      assignment[i] = min_j;
      for(k = 0; k < queryCol;k++) 
        query[i*queryCol +k] = queryNew[min_j*queryCol + k];

    // create new cluster
    } else {
      assignment[i] = m; // create assignment
      // save new cluster point
      for(k = 0; k < queryCol;k++) 
        queryNew[m*queryCol +k] = query[i*queryCol + k];
      m++; // update m
    }

  } 

  // ************************************************************************* 
  // Algorithm specific clean up 
  // ************************************************************************* 
  // K-d Tree
  // *********************************** 
  if( *algorithmEnumPtr == 1 ) {
      deleteTree( kdTree );
  }  


  free(bandwidth2);
  free(distances);
  free(neighbors);
  free(queryNew);

  return; 
}

