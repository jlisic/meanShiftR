#include <flann/flann.hpp>
#include <flann/io/hdf5.h>
#include <stdio.h> 
#include <stdlib.h>
#include <vector>

#include "R.h"
#include "Rinternals.h"
#define CSTACK_DEFNS 7
#include "Rinterface.h"

#include "Rmath.h"
#include "omp.h"





/**************** print matrix ******************************************/

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



/****************************** ADT stuff to fix issues with the compiler e.g. headers not working on all platforms *****/

/* create a queue ADT */
struct queue { 
  size_t item;
  struct queue * head;
  struct queue * tail;
};

typedef struct queue queue;
typedef struct queue * queuePtr;


/* function to create a new Q */
queuePtr createQueue( size_t item ) {

  queuePtr Q; /* declare our new pointer */

  Q = (queuePtr) malloc(sizeof(queue));
  //assert( Q != NULL );

  Q->head = NULL;
  Q->tail = NULL;
  Q->item = item;
  return(Q);
} 



/* function to remove a node in a queue */
/* here we remove qi from q */
queuePtr removeQueue( queuePtr q, queuePtr qi) { 

  //assert( q != NULL);
  //assert( qi != NULL);

  /* here we go towards the head, then look back at it's tail and replace it with qi's tail */ 
  if( qi->head != NULL)  {
    (qi->head)->tail = qi->tail;
  }
  else { 
    q = qi->tail;
  }

  if( qi->tail != NULL)  {
    (qi->tail)->head = qi->head;
  } 

  free(qi);
  qi = NULL;

  return(q);
}

/* function to delete a Q and set to NULL */
void deleteQueue( queuePtr q) {

while( q != NULL )
 q = removeQueue(q,q);  

return;
}



/* function to add a node before qi in a queue with a particular item */
/* if qi is null we assign to tail */
queuePtr addQueue(queuePtr q, queuePtr qi, size_t item) {
  queuePtr qNew,qs;
 
  /* create new queue */
  qNew = createQueue( item );

  /* in case there isn't anything there */
  if( q == NULL) {
    q = qNew;
    return(q);
  } 
  /* mess around with the tail node */ 
  else if( qi == NULL ) {
    /* we need to find the last pointer */
    qs = q;
    while(qs->tail != NULL) qs = qs->tail; 
    qNew->head = qs;
    qs->tail = qNew;
  }   
  /* mess around with the head node */ 
  else if( qi == q ) {
    qNew->tail = q;  /* qNew -> q */
    q->head = qNew;   
    q = qNew;
  }
  /* for everyone else */ 
  else {
    qNew->tail = qi;
    qNew->head = qi->head;
    (qi->head)->tail = qNew;
    qi->head   = qNew;
  }

  return(q); 
}



/* function to search a Q for an item, this function returns the first occurence */
queuePtr searchQueue( queuePtr q, size_t item) {
  queuePtr qi;

  //assert(q != NULL );

  qi = q;

  while( qi != NULL ){
    if( qi->item == item ) break; 
    qi = (qi->tail);
  } 

  return(qi);
}


/* function to create a queueArray */
queuePtr * queueArrayCreate( size_t * x, size_t n) {
 
  size_t i;

  queuePtr * queueArray = (queuePtr * ) calloc(n, sizeof(queuePtr));
  for( i = 0; i < n ; i++) queueArray[i] = createQueue(x[i]);    

  return( queueArray );
}


/* function to copy an element of a queueArray to another element */
void queueArrayCopy( queuePtr * queueArray, size_t from, size_t to ) {

  queuePtr updatePtr = queueArray[to]; 
  
  if( queueArray[from] == NULL) return; //Nothing to do
  
  if( queueArray[to] != NULL ) { 
    while( updatePtr->tail != NULL ) updatePtr=updatePtr->tail;
    updatePtr->tail = queueArray[from];  // add to chain
    (updatePtr->tail)->head = updatePtr;  // update head
  } else {
    queueArray[to] = queueArray[from];  // since to is NULL we just copy
  }

  queueArray[from] = NULL;   // set original to null

  return;
}


/* function to copy an element of a queueArray to another element */
/* but this has a check to track down where things have moved to */
void queueArrayCopyTogether( queuePtr * queueArray, size_t from, size_t to, size_t * check ) {

  queuePtr updatePtr = NULL; 
  
  if( queueArray[from] == NULL) {
    Rprintf("Invalid Move, %d has already moved to %d, therefore it can't mvoe to %d.\n", (int) from, (int) check[from], (int) to); 
    return; //Nothing to do
  }

  // this is the big difference from the regular copy, if we are moving to a null we figure out where it went
  if( queueArray[to] == NULL) to = check[to];  
 
  // something went wrong 
  if( queueArray[to] == NULL) { 
    Rprintf("Poor house keeping, %d is NULL.\n", (int) to); 
    return; //Nothing to do
  }
  
  //first we update check about where we are moving to 
  updatePtr = queueArray[from]; 
  check[ updatePtr->item ] = to;
  while( updatePtr->tail != NULL ) {
    updatePtr=updatePtr->tail;
    check[ updatePtr->item ] = to;
  }
  
  updatePtr = queueArray[to]; 
  while( updatePtr->tail != NULL ) updatePtr=updatePtr->tail;
  updatePtr->tail = queueArray[from];  // add to chain
  queueArray[from]->head = updatePtr;  // update head

  queueArray[from] = NULL;   // set original to null

  return;
}



/* function to create a queueArray */
void queueArrayPrint( queuePtr * x, size_t n) {
 
  size_t i;
  queuePtr updatePtr; 
  
  Rprintf("\n");

  for( i = 0; i < n ; i++) {
    updatePtr =  x[i];
    Rprintf("%d\t:", (int) i);
    while( updatePtr != NULL) {
      Rprintf("%d ", (int) updatePtr->item); 
      updatePtr=updatePtr->tail; 
    }
    Rprintf("\n");
  }

  return;
}


/* function to create a queueArray */
void queueArrayDelete( queuePtr * x, size_t n) {
 
  size_t i;
 
  /* clean up */ 
  for( i = 0; i < n ; i++) {
    if( x[i] != NULL) deleteQueue( x[i] ), x[i] = NULL; 
  }

  free(x);

  return;
}




/***********************************************************************************************/
/* Variable size array handling                                                                */
/***********************************************************************************************/
size_t *  updateVariableArray( size_t * x, size_t xRow , size_t  * xLength) {

  size_t * y;
  size_t j;

  // check if we need to resize
  if( xRow >= *xLength) {
    //Rprintf("resizing\n");
    y = (size_t *) calloc( xRow + *xLength/2 + 1, sizeof( size_t) );
    for(j = 0; j < *xLength; j++) y[j] = x[j]; 

    *xLength = xRow + *xLength/2 + 1; // update size
  
    // clean up  
    free(x);
    x = y;
    y = NULL;
  }
 
  // update 
  return(x);
}




/***********************************************************************************************/
/* Variable size matrix handling                                                               */
/***********************************************************************************************/
double *  updateVariableMatrix( double * x, size_t xRow, size_t xCol, size_t  * xLength ) {

  double * y;
  size_t j;

  // check if we need to resize
  if( xRow  >= *xLength) {
    //Rprintf("resizing\n");
    y = (double *) calloc( (xRow + *xLength/2 + 1)*xCol , sizeof(double) );
    for(j = 0; j < *xLength * xCol; j++) y[j] = x[j]; 

    *xLength = xRow + *xLength/2 + 1; // update size
  
    // clean up  
    free(x);
    x = y;
    y = NULL;
  }
 
  // update 
  return(x);
}



/***********************************************************************************************/
/* function for updating assignment, y, and x2y                                                */
/***********************************************************************************************/
    /* let's define the input and output here */
    /* INPUT */
    /* 
     * assignment     -- current group assignment        length xLength    size_t * 
     * queueArray     -- an array of queues              length xLength    queuePtr
     * xRow           -- length of x array               1                 size_t
     * xRow           -- length of y array               1                 size_t
     * x2y            -- mapping from x 2 y              length xLength    size_t *
     * x2yNew         -- mapping from x 2 updated y      length xLength    size_t *
     * neighborsQuery -- set of neighbors from z         length yLength    int *
     * distancesQuery -- set of distances from 1nn query length yLength    double *
     * prob           -- current kde value for a cluster length yLength    double *
     * mergeTreeIndex -- index for merge tree            length zLength?   int *
     * mergeTreeLength
     * mergeTreeSize
     * y              -- matrix of group locations       length yRow by xCol double *
     * z              -- matrix of historic group locations       length zRow by xCol double *
     */ 
     /* OUTPUT, destructive changes */
    /*
     * assignment     -- current group assignment        length xLength    size_t * 
     * x2yNew         -- mapping from x 2 updated y      length xLength    size_t *
     * y              -- matrix of group locations       length yRow by xCol double *
     * z              -- matrix of historic group locations       length zRow by xCol double *
     */

void mergePath( 
     size_t * assignment, 
     double ** zPtr,
     double * y,
     queuePtr * queueArray,
     size_t ** x2yPtr,
     size_t ** x2yNewPtr,
     size_t ** mergeTreeIndexPtr, 
     size_t * mergeTreeLengthPtr,
     size_t * mergeTreeSizePtr,
     size_t * zRowPtr,
     size_t * yRowPtr,
     size_t * zSizePtr,

     const double * epsilon,
     const size_t xRow,
     const size_t xCol,
     const int * neighborsQuery,
     const double * distancesQuery,
     const double * prob
    ) {

  size_t * x2y = *x2yPtr;
  size_t * x2yNew = *x2yNewPtr;
  size_t j, k;
  size_t j_y;        //variable to get the row number in yMatrix
  size_t j_neighbor;        // neighbor value in relation to x from queryTree
  size_t j_neighbor_y;     //j_neighbor in relation to y 

  double * z = *zPtr;
  size_t * mergeTreeIndex = *mergeTreeIndexPtr;

    *yRowPtr = 0;
    // find neighbors with higher prob within epsilon 
    for( j =0; j < xRow; j++) {
      if( assignment[j] == j ) {  // only run for points that haven't been merged

        if( queueArray[j] == NULL ) {
          Rprintf("\n\n############ j=%d Inconsistent State! ###############\n\n", (int) j);
          queueArrayPrint( queueArray, xRow); 
          return;
        }

        
        /* perform adjustments based on j */
        j_y = x2y[j];            // get the row number for j in Y
        j_neighbor = assignment[ mergeTreeIndex[ neighborsQuery[j_y] ] ];    // value for j's neighbor in x domain
        j_neighbor_y = x2y[j_neighbor]; // value for j's neighbor in Y
        
        // check if both sufficiently close and higher probability 
        if( 
            (assignment[j_neighbor] != j ) &         // check to see if the nearest neighbor's assignment isn't already j
            ( distancesQuery[j_y] < *epsilon ) & // check to see if distance threshold is met 
            ( prob[j_y] < prob[ j_neighbor_y ] )   // check to see if there is an increas in probability
          ) {

          // merge together  
          queueArrayCopyTogether( queueArray, j, j_neighbor, assignment);
          
        } else {
          /* check if we need to resize the vector */

          mergeTreeIndex = updateVariableArray( mergeTreeIndex, *mergeTreeLengthPtr + *yRowPtr, mergeTreeSizePtr );
          mergeTreeIndex[ *mergeTreeLengthPtr + *yRowPtr ] = j;
          

          z = updateVariableMatrix( z, *zRowPtr + *yRowPtr, xCol, zSizePtr );

          /* to ensure we only add unmerged points into the tree */ 
          /* note it's safe to write over y at this point */
          //Rprintf("  %d:\t", (int) (*zRowPtr + *yRowPtr) );
          for( k=0; k < xCol; k++) {
          //  Rprintf("%f ", y[ x2y[j]*xCol + k]);
            y[(*yRowPtr)*xCol +k] = y[ x2y[j] * xCol+k];
            z[(*zRowPtr + *yRowPtr)*xCol +k] = y[ x2y[j] * xCol+k];
          //  Rprintf("(%f) ", z[(*zRowPtr + *yRowPtr)*xCol +k] ); 
          }
          //Rprintf("\n");

          /* we can't write over x2y now since we will use it for later requets */
          /* so x2yNew is kept to refer to the updated y */ 
          x2yNew[j] = *yRowPtr;

          if( *yRowPtr > j ) {
            Rprintf("yRow is greater than j\n");
            return; 
          }
          if( *yRowPtr > x2y[j] ) {
            Rprintf("x2y[j] is greater than j\n");
            return; 
          }

          // increment 
          *yRowPtr = *yRowPtr + 1;

        }
      }
    }

    //update number of objects in query Index
    *mergeTreeLengthPtr = *mergeTreeLengthPtr + *yRowPtr;
    *zRowPtr = *zRowPtr + *yRowPtr;

    //update x2y
    *x2yPtr = x2yNew;
    *x2yNewPtr = x2y;
 
    /* pass pointers back */ 
    *zPtr = z;
    *mergeTreeIndexPtr = mergeTreeIndex;
    
    //Rprintf("z - postMerge1:\n"); 
    //printMatrixFullDbl(z,*zRowPtr,xCol);

    return;
}



/***********************************************************************************************/
/* this is an alternative kernel                                                               */
/***********************************************************************************************/
void normalKernelNewton(
  double * query,        /* query matrix */
  double * train,        /* train matrix */
  int * neighbors,      /* neighbors matrix */
  double * prob,
  size_t queryRow,
  size_t queryCol,
  size_t nNeighbors,
  double * bandwidth,    /* bandwidth parameters */
  double alpha, 
  double * weights // for debugging
  ) {

  size_t i; /* query data index */
  size_t j; /* build data index */ 
  size_t k; /* neighbor index */
  size_t d; /* column index */


  double * numerator;
  double w, w0;
  double * denominator;
  double denominatorMS;
  double * z;
  int currentNeighbor;

/*
6 # calculate the derivative
7 dnorm.deriv1 <- function(x,X,h) {
8   n <- nrow(X)
9   -1/(n * h^3 ) * t(x - X) %*% dnorm( (x - X)/h )
10 }
11
12 # calculate the original function
13 dnorm.deriv2 <- function(x,X,h) {
14   n <- nrow(X)
15   1/(n * h^5 ) * t( ( x - X)^2 - h^2)  %*% dnorm( (x - X)/h )
16 }
*/ 

//  omp_set_num_threads(1);

  /* iterate over each query point */ 
  #pragma omp parallel private(i,d,k,z,w,numerator,denominator,denominatorMS,currentNeighbor,w0)
  {
  
  z         = (double * ) calloc( queryCol, sizeof(double));
  numerator = (double * ) calloc( queryCol, sizeof(double));
  denominator = (double * ) calloc( queryCol, sizeof(double));
  

  #pragma omp for
  for( i = 0; i < queryRow; i++)  {
    prob[i] = 0;

    /* fotal weight to divide by for the denominator */
    for( d = 0; d < queryCol; d++) { 
      numerator[d] = 0; 
      denominator[d] = 0; 
    }
    denominatorMS = 0;

    /* iterate over each neighbor */
    for( k = 0; k < nNeighbors; k++)  { 
      
      currentNeighbor = neighbors[i*nNeighbors + k];

      /* set weight to 1 */ 
      w = 1;

      // for each point query over each dimension 
      for( d = 0; d < queryCol; d++)  { 
        
        // save a copy of the point 
        z[d] = ( query[i * queryCol + d] - train[currentNeighbor * queryCol + d] )/bandwidth[d];
        
        // double dnorm( double x, double mu, double sigma, int give_log)   
        //w *=  dnorm( z[d] , 0, 1 , 0 );
        w *=  0.5 * exp( - 0.5 * z[d]*z[d] );

      }

      // for debug
      if( weights != NULL) weights[i * nNeighbors + k] = w;

      /* aggregate over all variables */
      for( d = 0; d < queryCol; d++) {
        prob[i] += w;
        numerator[d] +=  z[d] * bandwidth[d] * w ;
        if( alpha > 0 ) {
          denominator[d] += ( 1 - alpha * z[d] * z[d]) * w ;
        } else { 
          denominator[d] += w;
        }
      }

/*
    if( k > 0 ) {
        if( (w > 0) & (w / w0 < 0.0001 ) ) break;  
      } else { 
        w0 = w;
      }
*/
     
      
      denominatorMS += w; 
    }

    for( d = 0; d < queryCol; d++) { 
      if( denominator[d] != 0 )  query[i* queryCol + d] = query[i*queryCol +d] - numerator[d] / denominator[d];
    }
     
  }
  
  free(z);
  free(denominator); 
  free(numerator);

  } // for omp

  return;
}


/* 
 * y points to find neighbors for 
 * z pool of neighbors 
 */
void getOneNeighbor( 
    double * y, 
    double * z, 
    size_t yRow, 
    size_t zRow, 
    size_t yCol, 
    int * neighbors, 
    double * distance
  ) {
  
  flann::SearchParams exactSearchParameters;
  //exactSearchParameters.eps = 0.0001;
  exactSearchParameters.checks = 128;
  exactSearchParameters.cores = 0;
  
  flann::Matrix<double> yMatrix( y, yRow, yCol ); 
  flann::Matrix<double> zMatrix( z, zRow, yCol ); 
  flann::Matrix<double> distancesQueryMatrix; 
  flann::Matrix<int> neighborsQueryMatrix; 
  neighborsQueryMatrix = flann::Matrix<int>(neighbors, yRow, 1);
  distancesQueryMatrix = flann::Matrix<double>(distance, yRow, 1);
 
  /* create Index */ 
  flann::Index< flann::L2<double> > queryIndex( zMatrix, flann::KDTreeIndexParams( 4 ) ); 
  queryIndex.buildIndex();

  queryIndex.knnSearch(
      yMatrix, 
      neighborsQueryMatrix, 
      distancesQueryMatrix, 
      1, 
      exactSearchParameters
      ); 

  return;
} 




/*
R side, x:
 1  5  9  x1
 2  6 10  x2
 3  7 11  x3
 4  8 12  x4

x:
x1 x2 x3 x4
 1  2  3  4
 5  6  7  8
 9 10 11 12

*/

// extern start for .C interface
extern "C" {


/* mean shfit nearest neighbors */
void R_meanShiftNN(
  double * x,                /* data to query for */
  double * r,                /* reference data */
  int * neighbors,
  double * distances,
  double * prob,
  int * neighborsQuery,              /* assignment */
  double * distancesQuery,
  int * xRowPtr,                     /* number of rows of data to query */
  int * rRowPtr,                     /* number of rows of data for kde */
  int * xColPtr,                     /* number of columns of data to form r */
  int * nNeighborsPtr,               /* number of Neighbors */
  int * intSearchParameters,         /* integer parameters for search*/
  double * doubleSearchParameters,   /* double parameters for search */ 
  int * intAlgorithmParameters,      /* integer parameters for algorithm*/
  double * doubleAlgorithmParameters,/* double parameters for algorithm*/
  int * kernelMethodPtr,             /* kernel methods */ 
  double * bandwidth,                /* bandwidth for kernel methods */
  double * alpha,                    /* alpha - for newton approximation */
  int * iterations,                  /* number of iterations */
  double * epsilon,                  /* cut off for classification */
  double * interpolate,
  int * interpolateIndex,
  int * interpolateRowPtr,
  int * assignmentsDebug,
  double * weightsDebug,
  int * neighborsDebug,              /* record all neighbors for each r record */
  double * valuesDebug,
  double * clusterEpsilon            /* cut off for cluster distances */  

) {

  /* simplify the interface by dereferencing the pointers */
  size_t xRow =     (size_t) * xRowPtr; 
  size_t rRow = (size_t) * rRowPtr; 
  size_t xCol =     (size_t) * xColPtr;
  size_t kernelMethod =   (size_t) * kernelMethodPtr;
  size_t nNeighbors = (size_t) * nNeighborsPtr;
  size_t interpolateRow = (size_t) * interpolateRowPtr;

  size_t debug = (size_t) assignmentsDebug[1];   // status of recording observations via debug

  /* indexes */
  size_t yRow = xRow; // number of rows for the update matrix for the QueryTree 

  /* calculate initial size for dynamic arrays */
  size_t mergeTreeSize = 10 * xRow + 1;  // current max number of elements in the variable size array, treeIndex
  size_t * mergeTreeIndex = (size_t *) calloc( mergeTreeSize, sizeof( size_t) ); // variable size array used to map elements in queryTree to x
  size_t mergeTreeLength = xRow;   // current number of elements of queryTree
 
  size_t zSize = 10 * xRow + 1;  // current max number of elements in the variable size matrix, treeIndex
  double * z = (double *) calloc( zSize *  xCol, sizeof( double) ); // variable size matrix used to map elements in queryTree to x
  size_t zRow = xRow;   // current number of elements of queryTree
  
  /* create Matrix ADT for adding to second tree */
  double * y = (double *) calloc( xRow * xCol, sizeof(double) ); 

  size_t i,j,k; /* index */
  size_t j_neighbor;
  size_t j_neighbor_y;

  /* for interpolation */
  double * interpolateDistance;

  /* fixes for R with OpenMP */
  R_CStackLimit=(uintptr_t)-1;

  /* create flann matricies */
  flann::Matrix<double> yMatrix; 
  flann::Matrix<double> distancesMatrix; 
  flann::Matrix<int> neighborsMatrix; 
  
  flann::Matrix<double> interpolateMatrix; 
  flann::Matrix<double> distancesInterpolateMatrix; 
  flann::Matrix<int> neighborsInterpolateMatrix; 

  flann::Matrix<double> rMatrix( r, rRow, xCol ); 
  
  /* index to determine if something has been assigned */
  size_t * assignment = (size_t * ) calloc( xRow, sizeof(size_t));
  size_t * x2y = (size_t * ) calloc( xRow, sizeof(size_t));
  size_t * x2yNew = (size_t * ) calloc( xRow, sizeof(size_t));

  /* check min cluster distance */
  queuePtr tmpPtr;
  
  /* copy data over for y and x */
  for( i = 0; i < xRow * xCol; i++) {
    y[i] = x[i];
    z[i] = x[i];
  }
  /* for everything that needs an increasing initial value */
  for( i = 0; i < xRow; i++) assignment[i] = i, x2y[i] = i, x2yNew[i]=i;
  for( i = 0; i < rRow; i++) mergeTreeIndex[i]=i;
  queuePtr * queueArray = queueArrayCreate( assignment, xRow); 
  queuePtr * queueArrayY = (queuePtr *) calloc( xRow, sizeof( queuePtr) ); 
     
  // ----- debugging variables  
  // offsets for debugging 
  size_t offset ;  // obs * iteration
  size_t offsetNeighbors;
  size_t offsetCol;
  double * weights = NULL;
  if( debug == 1) weights = (double *) calloc( nNeighbors * xRow, sizeof( double ) );  // create weight matrix for debugging 

  double dist;

  // double epsilon since flann uses the unsquared difference
  *epsilon = (*epsilon) * (*epsilon);
  *clusterEpsilon = (*clusterEpsilon) * (*clusterEpsilon);

  /******************* CREATE INDEX *********************/ 
  /* if KDTreeIndex */ 
  flann::Index< flann::L2<double> > index( rMatrix, flann::KDTreeIndexParams( intAlgorithmParameters[0] ) ); 

  /* build the index */ 
  index.buildIndex();

  /******************* SEARCH *********************/ 
  flann::SearchParams mpSearchParameters;
  mpSearchParameters.cores = 0;
  mpSearchParameters.checks = intSearchParameters[0];
  
     
  // create initial xMatrix 
  yMatrix = flann::Matrix<double>(y, yRow, xCol);
  neighborsMatrix = flann::Matrix<int>(neighbors, yRow, nNeighbors);
  distancesMatrix = flann::Matrix<double>(distances, yRow, nNeighbors);
    
//  Rprintf("x:\n"); 
//  printMatrixFullDbl(x,xRow,xCol);

  /* perform searches */ 
  /* Implementation Notes
   * yRow is a monotonically decreaseing sequence based on consecutive merges.
   * Therefore the following items shrink (this may just be with respect to using less rows):
   *   y
   *   yMatrix
   *   neighborsMatrix
   *   distancesMatrix
   *   yRow
   *   neighborsMatrixQuery
   *   distancesQueryMatrix
   *   prob
   * 
   * The following things stay the same
   *   r
   *   xCol
   *   nNeighbors
   *   bandwidth
   *   alpha
   *   x2y
   *   assignment
   *   queueArray
   *   x
   *
   * The following things grow in size
   *   mergeTreeIndex
   *   queryIndex
   *
   */
  for( i = 0; i < *iterations; i++ ) {
//Rprintf("Iteration: %d ----------------------------------------------------------- \n", (int) i);
   

    /* look for points in the ring set */
    //Rprintf("KNN Search: %d ----------------------------------------------------------- \n", (int) i);
    index.knnSearch(yMatrix, neighborsMatrix, distancesMatrix, nNeighbors, mpSearchParameters); 
 
    if(kernelMethod != 3) {
     Rprintf("This kernel is not supported");
    } 
  
    //Rprintf("Starting Kernel: %d ----------------------------------------------------------- \n", (int) i);
    if( kernelMethod == 3) {
      normalKernelNewton(
        y,            /* query matrix */
        r,            /* reference matrix */
        neighbors,    /* neighbors matrix */
        prob,
        yRow,
        xCol,
        nNeighbors,
        bandwidth,    /* bandwidth parameters */
        *alpha,
        weights  // for debug
      );
    }
  
    // query  querytree and find closest point 
    getOneNeighbor(y, z, yRow, zRow, xCol, neighborsQuery, distancesQuery); 

    /***************** DEBUG *****************************/ 
    // record record level changes for each iteration
    
    if( debug == 1) {
     
      offset = xRow * i; 
      offsetNeighbors = offset * nNeighbors;
      offsetCol = offset * xCol;

      // copy assignment back over for debug
      //#pragma omp for private(k,j_neighbor_y)
      for( j = 0; j < xRow; j++) {

        // handle assignment
        assignmentsDebug[ offset + j] = (int) assignment[j];

        // handle weights and neighbors 
        j_neighbor_y = x2y[assignment[j]];   // get mapping of current record to y

        for( k = 0; k < nNeighbors; k ++) {
          neighborsDebug[ offsetNeighbors +  j*nNeighbors + k] = neighbors[ j_neighbor_y * nNeighbors +k];  
          weightsDebug[ offsetNeighbors +  j*nNeighbors + k]   = weights[ j_neighbor_y * nNeighbors +k];  
        }

        for( k = 0; k < xCol; k++ ) valuesDebug[ offsetCol +  j*xCol + k] = y[ j_neighbor_y * xCol +k];  

      }
    }
    /***************** DEBUG *****************************/ 


     /* perform merge path */
     mergePath( 
       assignment, 
       &z,
       y,
       queueArray,
       &x2y,
       &x2yNew,
       &mergeTreeIndex,
       &mergeTreeLength,
       &mergeTreeSize,
       &zRow,
       &yRow,
       &zSize,
       epsilon,
       xRow,
       xCol,
       neighborsQuery,
       distancesQuery,
       prob
    ); 
/* 
    if( i >= 9 ) {
      printMatrixFullSize_t(assignment, xRow , 1 );
      printMatrixFullDbl(y, yRow , xCol );
    }
*/
    
    Rprintf("Iteration: %d Decrease = %04.2f\n", (int) i + 1, 1 - (float) yRow / (float) xRow );
          
    if( i < *iterations - 1 ) { // don't do a final query

      // update matricies
      yMatrix = flann::Matrix<double>(y, yRow, xCol);
      neighborsMatrix = flann::Matrix<int>(neighbors, yRow, nNeighbors);
      distancesMatrix = flann::Matrix<double>(distances, yRow, nNeighbors);
    
    }

  }  // end iteraations

  // figure out what assignments we may need to change
  for( i=0; i < xRow; i++) {
    if( queueArray[i] != NULL) {
      k= x2y[i];
      queueArrayY[k] = queueArray[i]; 
    } 
  }



  // do a final merge
  for( i=0; i < yRow; i++) {
    for( j=i+1; j < yRow; j++) {

      // get distance
      dist = 0;
      for( k=0; k < xCol; k++) dist += (y[i*xCol + k] - y[j*xCol+k])* (y[i*xCol + k] - y[j*xCol+k]);
      
      // check distance
      if (*clusterEpsilon >= dist) {

        k = assignment[ queueArrayY[i]->item ];
        tmpPtr = queueArrayY[j]; 

        // update assignment
        while( tmpPtr != NULL ){        
          assignment[tmpPtr->item] = k;
          tmpPtr = tmpPtr->tail;
        }

      }
    }
  }  
  
 

  // interpolate
  if( interpolateIndex[0] != -1) { 

    interpolateDistance = (double * ) calloc(interpolateRow, sizeof(double)); //create some space
    
    interpolateMatrix = flann::Matrix<double>(interpolate, interpolateRow, xCol);
    neighborsInterpolateMatrix = flann::Matrix<int>(interpolateIndex, interpolateRow, 1);
    distancesInterpolateMatrix = flann::Matrix<double>(interpolateDistance, yRow, 1);

    index.knnSearch(
        interpolateMatrix, 
        neighborsInterpolateMatrix, 
        distancesInterpolateMatrix, 
        1, 
        mpSearchParameters
        ); 
  
    // copy assignment back over for interpolate
    #pragma omp for private(k,j_neighbor,j_neighbor_y)
    for( i = 0; i < interpolateRow; i++) {
      //j_neighbor    = treeIndex[ interpolateIndex[i] ];  // get neighbor
      j_neighbor    = interpolateIndex[i];  // get neighbor
      j_neighbor_y = x2y[ assignment[j_neighbor] ]; // get neighbor in Y
    
      interpolateIndex[i] = assignment[j_neighbor];

      // copy back location
      for( k = 0; k < xCol; k ++) interpolate[i*xCol + k] = y[ j_neighbor_y * xCol +k];  
  
    }

    free(interpolateDistance);
  }


  // copy assignment back over
  #pragma omp for private(k,j)
  for( i = 0; i < xRow; i++) {
    neighborsQuery[i] = (int) assignment[i];

    j = x2y[assignment[i]];

    // copy back location
    for( k = 0; k < xCol; k ++) x[i*xCol + k] = y[ j * xCol +k];  

    // using distance to copy back
    distancesQuery[i] = prob[ j ];     
  }



  free(z);
  free(weights);
  free(assignment);
  free(y);
  free(x2y);
  free(x2yNew);
  free(mergeTreeIndex);

  free(queueArrayY);
  for( i = 0; i < xRow; i++) deleteQueue( queueArray[i] );
  free(queueArray);
  Rprintf("\n");

  return; 
}


// extern end for .C interface
}



