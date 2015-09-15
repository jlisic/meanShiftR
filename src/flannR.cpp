#include <flann/flann.hpp>
#include <flann/io/hdf5.h>
#include <stdio.h> 
#include <stdlib.h>

#include "R.h"
#include "Rinternals.h"
#define CSTACK_DEFNS 7
#include "Rinterface.h"

#include "Rmath.h"
#include "omp.h"


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
    Rprintf("resizing\n");
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
  if( xRow >= *xLength) {
    Rprintf("resizing\n");
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
        z[d] = ( query[i * queryCol + d] - train[currentNeighbor * queryCol + d] );
        
        // double dnorm( double x, double mu, double sigma, int give_log)   
        w *=  dnorm( z[d] , 0, 1 , 0 );

      }

      // for debug
      if( weights != NULL) weights[i * nNeighbors + k] = w;

      /* aggregate over all variables */
      for( d = 0; d < queryCol; d++) {
        prob[i] += w;
        numerator[d] +=  z[d] * w ;
        denominator[d] += ( 1 - alpha * z[d] * z[d]) * w ;
      }

      if( k > 0 ) {
        if( (w > 0) & (w / w0 < 0.0001 ) ) break;  
      } else { 
        w0 = w;
      }
     
      
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
 * x points to find neighbors for 
 * y pool of neighbors 
 */
void getNeighbors( 
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

  queryIndex.knnSearch(yMatrix, neighborsQueryMatrix, distancesQueryMatrix, 1, exactSearchParameters); 

  return;
}; 












// extern start for .C interface
extern "C" {


/* mean shfit nearest neighbors */
void R_meanShiftNN(
  double * x,                /* data to query for */
  double * train, 
  int * neighbors,
  double * distances,
  double * prob,
  int * neighborsQuery,
  double * distancesQuery,
  int * xRowPtr,                     /* number of rows of data to query */
  int * trainRowPtr,                 /* number of rows of data to train on */
  int * xColPtr,                     /* number of columns of data to form train */
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
  int * neighborsDebug,              /* record all neighbors for each train record */
  double * valuesDebug       

) {

  /* simplify the interface by dereferencing the pointers */
  size_t xRow =     (size_t) * xRowPtr; 
  size_t trainRow = (size_t) * trainRowPtr; 
  size_t xCol =     (size_t) * xColPtr;
  size_t kernelMethod =   (size_t) * kernelMethodPtr;
  size_t nNeighbors = (size_t) * nNeighborsPtr;
  size_t interpolateRow = (size_t) * interpolateRowPtr;

  size_t debug = (size_t) assignmentsDebug[1];   // status of recording observations via debug
  size_t inNeighbors;

  /* indexes */
  int jAdjust;        // neighbor value in relation to x from queryTree
  int jAdjustInY;     //jAdjust in relation to y 
  size_t jInY;        //variable to get the row number in yMatrix
  size_t yRow = xRow; // number of rows for the update matrix for the QueryTree 


  size_t currentTreeIndexLength = 10* xRow + 1;  // current max number of elements in the variable size array, treeIndex
  size_t * treeIndex = (size_t *) calloc( currentTreeIndexLength, sizeof( size_t) ); // variable size array used to map elements in queryTree to x
  size_t currentTreeIndexMax = trainRow;   // current number of elements of queryTree
  
  size_t currentZIndexLength = 10* xRow + 1;  // current max number of elements in the variable size matrix, treeIndex
  double * z = (double *) calloc( currentZIndexLength * xCol, sizeof( double) ); // variable size matrix used to map elements in queryTree to x
  size_t zRow = trainRow;   // current number of elements of queryTree

  size_t i,j,k; /* index */

  double * tmpPtr;

  /* for interpolation */
  double * interpolateDistance;

  /* fixes for R with OpenMP */
  R_CStackLimit=(uintptr_t)-1;

  /* create Matrix ADT for adding to second tree */
  double * y = (double *) calloc( xRow * xCol, sizeof(double) ); 
  for( i = 0; i < xRow * xCol; i++) y[i] = x[i];
  for( i = 0; i < trainRow * xCol; i++) z[i] = train[i];

  flann::Matrix<double> yMatrix; 
  flann::Matrix<double> distancesMatrix; 
  flann::Matrix<int> neighborsMatrix; 
  
  flann::Matrix<double> interpolateMatrix; 
  flann::Matrix<double> distancesInterpolateMatrix; 
  flann::Matrix<int> neighborsInterpolateMatrix; 

  flann::Matrix<double> trainMatrix( train, trainRow, xCol ); 
  
  /* index to determine if something has been assigned */
  size_t * assignment = (size_t * ) calloc( xRow, sizeof(size_t));
  size_t * x2y = (size_t * ) calloc( xRow, sizeof(size_t));
  size_t * x2yNew = (size_t * ) calloc( xRow, sizeof(size_t));
  size_t * x2yTmp = NULL;
  
  /* for everything that needs an increasing initial value */
  for( i = 0; i < xRow; i++) assignment[i] = i, x2y[i] = i, x2yNew[i]=i;
  for( i = 0; i < trainRow; i++) treeIndex[i]=i;
  queuePtr * queueArray = queueArrayCreate( assignment, xRow); 
     
  // ----- debugging variables  
  // offsets for debugging 
  size_t offset ;  // obs * iteration
  size_t offsetNeighbors;
  size_t offsetCol;
  double * weights = NULL;
  if( debug == 1) weights = (double *) calloc( nNeighbors * xRow, sizeof( double ) );  // create weight matrix for debugging 

  double dist;
  //FILE * currentFile;
  //FILE * queryFile;
  //FILE * query2File;
  //FILE * query3File;
  //currentFile = fopen("query.txt","w");
  //queryFile = fopen("queryTree.txt","w");
  //query2File = fopen("query2Tree.txt","w");
  //query3File = fopen("query3Tree.txt","w");


  /******************* CREATE INDEX *********************/ 
  /* if KDTreeIndex */ 
  flann::Index< flann::L2<double> > index( trainMatrix, flann::KDTreeIndexParams( intAlgorithmParameters[0] ) ); 

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
   *   train
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
   *   treeIndex
   *   queryIndex
   *
   */
  for( i = 0; i < *iterations; i++ ) {
//Rprintf("Iteration: %d ----------------------------------------------------------- \n", (int) i);
   

    /* look for points in the training set */
    //Rprintf("KNN Search: %d ----------------------------------------------------------- \n", (int) i);
    index.knnSearch(yMatrix, neighborsMatrix, distancesMatrix, nNeighbors, mpSearchParameters); 
 
    if(kernelMethod != 3) {
     Rprintf("This kernel is not supported");
    } 
  
    //Rprintf("Starting Kernel: %d ----------------------------------------------------------- \n", (int) i);
    if( kernelMethod == 3) {
      normalKernelNewton(
        y,            /* query matrix */
        train,        /* train matrix */
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
    
    
    // start off with all elements merged 
    
    /***************** DEBUG *****************************/ 
    /*
    fprintf( currentFile, "PRESEARCH \n");
    for( j =0; j < xRow; j++) {

      jInY = x2y[assignment[j]];

      fprintf(currentFile, "y: %d, %d, %d ", 
           (int) i,
           (int) j,
           (int) jInY
       ); 
      for( k=0; k < xCol; k++) fprintf( currentFile, "%f,", y[jInY*xCol +k] );
      fprintf( currentFile, "\n");
    }
    fprintf( currentFile, "\n");
    */
    /***************** DEBUG *****************************/ 

    // query  querytree and find closest point 
    getNeighbors(y, z, yRow, zRow, xCol, neighborsQuery, distancesQuery); 
   


    /***************** DEBUG *****************************/ 
    // record record level changes for each iteration
    
    if( debug == 1) {
     
      offset = xRow * i; 
      offsetNeighbors = offset * nNeighbors;
      offsetCol = offset * xCol;

      // copy assignment back over for debug
      //#pragma omp for private(k,jAdjustInY)
      for( j = 0; j < xRow; j++) {

        // handle assignment
        assignmentsDebug[ offset + j] = (int) assignment[j];

        // handle weights and neighbors 
        jAdjustInY = x2y[assignment[j]];   // get mapping of current record to y

        for( k = 0; k < nNeighbors; k ++) {
          neighborsDebug[ offsetNeighbors +  j*nNeighbors + k] = neighbors[ jAdjustInY * nNeighbors +k];  
          weightsDebug[ offsetNeighbors +  j*nNeighbors + k]   = weights[ jAdjustInY * nNeighbors +k];  
        }

        for( k = 0; k < xCol; k++ ) valuesDebug[ offsetCol +  j*xCol + k] = y[ jAdjustInY * xCol +k] * bandwidth[k];  

      }
    }
    
    /***************** DEBUG *****************************/ 

    /***************** DEBUG *****************************/ 
    /*
    // log 
    for( j =0; j < currentTreeIndexMax; j++) {

      if( neighborsQuery[jInY] >=  0 ) {
      tmpPtr = queryIndex.getPoint( j );
      fprintf(queryFile, "%d, %d, %d ", 
           (int) i,
           (int) j,
           (int) treeIndex[j] 
       ); 
      for( k=0; k < xCol; k++) fprintf( queryFile, "%f,", tmpPtr[k] );
      fprintf( queryFile, "\n");
      } else {
        fprintf(queryFile, "%d \n", neighborsQuery[jInY]);
      }

    }
    fprintf(queryFile, "\n");
    */
    /***************** DEBUG *****************************/ 

    /***************** DEBUG *****************************/ 
    // log
    /*
    for( j =0; j < xRow; j++) {

      jInY = x2y[assignment[j]];

      fprintf(currentFile, "y: %d, %d, %d ", 
           (int) i,
           (int) j,
           (int) jInY
       ); 
      for( k=0; k < xCol; k++) fprintf( currentFile, "%f,", y[jInY*xCol +k] );
      fprintf( currentFile, "\n");

      if( neighborsQuery[jInY] >=  0 ) {
      tmpPtr = queryIndex.getPoint( neighborsQuery[jInY] );
      fprintf(currentFile, "z: %d, %d, %d ", 
           (int) i,
           (int) j,
           (int) jInY 
       ); 
      for( k=0; k < xCol; k++) fprintf( currentFile, "%f,", tmpPtr[k] );
      fprintf( currentFile, "\n");
      } else {
        fprintf(currentFile, "z: %d \n", neighborsQuery[jInY]);
      }

    }
    fprintf(currentFile, "\n");
    */
    /***************** DEBUG *****************************/ 
    /*
    Rprintf("X2Y:\n");
    for( j =0; j < xRow; j++) Rprintf("%d ",j);
    Rprintf("\n");
    for( j =0; j < xRow; j++) Rprintf("%d ", x2y[assignment[j]]);
    Rprintf("\n");
   
    Rprintf("Assignment:\n");
    for( j =0; j < xRow; j++) Rprintf("%d ",j);
    Rprintf("\n");
    for( j =0; j < xRow; j++) Rprintf("%d ", assignment[j]);
    Rprintf("\n");
    
    Rprintf("Tree Index:\n");
    for( j =0; j < currentTreeIndexMax; j++) Rprintf("%d: %d\n", (int) j, (int) treeIndex[j]);
    Rprintf("\n");
    
    Rprintf("NeighborsQuery:\n");
    for( j =0; j < yRow; j++) Rprintf("%d: %d\n", (int) j, (int) neighborsQuery[j]);
    Rprintf("\n");
    */
    
    yRow = 0;
    // find neighbors with higher prob within epsilon 
    for( j =0; j < xRow; j++) {
      if( assignment[j] == j ) {  // only run for points that haven't been merged

        if( queueArray[j] == NULL ) {
          Rprintf("\n\n############ j=%d Inconsistent State! ###############\n\n", (int) j);
          queueArrayPrint( queueArray, xRow); 
          return;
        }

        /* perform adjustments based on j */
        jInY = x2y[j];            // get the row number for j in Y
        jAdjust = assignment[ treeIndex[ neighborsQuery[jInY] ] ];    // value for j's neighbor
        jAdjustInY = x2y[jAdjust]; // value for j's neighbor in Y


// DEBUG 
//Rprintf("%d: j=%d x2y=%d neighbor=%d migrateTo= %d (%f %f)", (int) i, (int) j, (int) jInY, (int) neighborsQuery[jInY], (int) treeIndex[neighborsQuery[jInY]], distancesQuery[jInY], *epsilon);
//Rprintf("[%f, %f]\n", prob[j], prob[ jAdjustInY ] );
// DEBUG 

        
        // check if both sufficiently close and higher probability 
        if( 
            (assignment[jAdjust] != j ) &         // check to see if the nearest neighbor's assignment isn't already j
            ( distancesQuery[jInY] < *epsilon ) & // check to see if distance threshold is met 
            ( prob[jInY] < prob[ jAdjustInY ] )   // check to see if there is an increas in probability
          ) {

          // merge together  
//Rprintf("\tj = %d, jAdjust = %d \n", (int) j, (int) jAdjust);
          queueArrayCopyTogether( queueArray, j, jAdjust, assignment);
          
        } else {
          /* check if we need to resize the vector */

//Rprintf("currentTreeIndexLength = %d, add = %d, currentTreeIndexMax = %d\n", currentTreeIndexMax, currentTreeIndexMax + yRow, currentTreeIndexLength);
          treeIndex = updateVariableArray( treeIndex, currentTreeIndexMax + yRow, &currentTreeIndexLength );
//Rprintf("currentTreeIndexLength = %d, add = %d, currentTreeIndexMax = %d\n", currentTreeIndexMax, currentTreeIndexMax + yRow, currentTreeIndexLength);
          treeIndex[ currentTreeIndexMax + yRow ] = j;

//Rprintf("zRow = %d, add = %d, currentZIndexMax = %d\n", zRow, zRow + yRow, currentZIndexLength);
          z = updateVariableMatrix( z, zRow + yRow, xCol, &currentZIndexLength );
//Rprintf("zRow = %d, add = %d, currentZIndexMax = %d\n", zRow, zRow + yRow, currentZIndexLength);

          /* to ensure we only add unmerged points into the tree */ 
          for( k=0; k < xCol; k++) {
            y[yRow*xCol +k] = y[ x2y[j] * xCol+k];
            z[(zRow + yRow)*xCol +k] = y[ x2y[j] * xCol+k];
          }
          x2yNew[j] = yRow;

          if( yRow > j ) {
            Rprintf("yRow is greater than j\n");
            return; 
          }
          if( yRow > x2y[j] ) {
            Rprintf("x2y[j] is greater than j\n");
            return; 
          }

          // increment 
          yRow++;

        }
      }
    }
    //Rprintf("Iteration: %d Decrease = %04.2f\r", (int) i + 1, 1 - (float) yRow / (float) xRow );
    Rprintf("Iteration: %d Decrease = %04.2f\n", (int) i + 1, 1 - (float) yRow / (float) xRow );

    //update number of objects in query Index
    currentTreeIndexMax += yRow;
    zRow += yRow;

    //update x2y
    x2yTmp = x2y;
    x2y = x2yNew;
    x2yNew = x2yTmp;
    x2yTmp = NULL; 
          
//queueArrayPrint( queueArray, xRow); 
      
    if( i < *iterations - 1 ) { // don't do a final query

      // update matricies
      yMatrix = flann::Matrix<double>(y, yRow, xCol);
      neighborsMatrix = flann::Matrix<int>(neighbors, yRow, nNeighbors);
      distancesMatrix = flann::Matrix<double>(distances, yRow, nNeighbors);
    
      
    /***************** DEBUG *****************************/ 
    /*
      for( j =0; j < yRow; j++) {
        fprintf(query2File, "y: %d, %d, ", 
             (int) i,
             (int) j
         ); 
        for( k=0; k < xCol; k++) fprintf( query2File, "%f,", y[j*xCol +k] );
        fprintf( query2File, "\n");
      }
      fprintf( query2File, "\n");
    */
    /***************** DEBUG *****************************/ 
 

    /***************** DEBUG *****************************/ 
    // log 
    /*
    for( j =0; j < currentTreeIndexMax; j++) {

      tmpPtr = queryIndex.getPoint( j );
      fprintf(query3File, "%d, %d, %d ", 
           (int) i,
           (int) j,
           (int) treeIndex[j] 
       ); 
      for( k=0; k < xCol; k++) fprintf( query3File, "%f,", tmpPtr[k] );
      fprintf( query3File, "\n");

    }
    fprintf(query3File, "\n");
    */
    /***************** DEBUG *****************************/ 


    }
  

  }  // end iteraations
  
  // interpolate
  if( interpolateIndex[0] != -1) { 

    interpolateDistance = (double * ) calloc(interpolateRow, sizeof(double)); //create some space
    
    interpolateMatrix = flann::Matrix<double>(interpolate, interpolateRow, xCol);
    neighborsInterpolateMatrix = flann::Matrix<int>(interpolateIndex, interpolateRow, 1);
    distancesInterpolateMatrix = flann::Matrix<double>(interpolateDistance, yRow, 1);

    index.knnSearch(interpolateMatrix, neighborsInterpolateMatrix, distancesInterpolateMatrix, 1, mpSearchParameters); 
  
    // copy assignment back over for interpolate
    #pragma omp for private(k,jAdjust,jAdjustInY)
    for( i = 0; i < interpolateRow; i++) {
      //jAdjust    = treeIndex[ interpolateIndex[i] ];  // get neighbor
      jAdjust    = interpolateIndex[i];  // get neighbor
      jAdjustInY = x2y[ assignment[jAdjust] ]; // get neighbor in Y
    
//      if( i < 10) Rprintf("%d: jAdjustInY %d [%d] jAdjustInY %d [%d]\n" , (int) i, (int) jAdjust, (int) xRow, (int) jAdjustInY, (int) yRow);

      interpolateIndex[i] = assignment[jAdjust];

      // copy back location
      for( k = 0; k < xCol; k ++) interpolate[i*xCol + k] = y[ jAdjustInY * xCol +k] * bandwidth[k];  
  
    }

    free(interpolateDistance);
  }


  // copy assignment back over
  #pragma omp for private(k,j)
  for( i = 0; i < xRow; i++) {
    neighborsQuery[i] = (int) assignment[i];

    j = x2y[assignment[i]];

    // copy back location
    for( k = 0; k < xCol; k ++) x[i*xCol + k] = y[ j * xCol +k] * bandwidth[k];  

    // using distance to copy back
    distancesQuery[i] = prob[ j ];     
  }


  // free everything
  //fclose(currentFile);
  //fclose(queryFile);
  //fclose(query2File);
  //fclose(query3File);

  free(z);
  free(weights);
  free(assignment);
  free(y);
  free(x2y);
  free(x2yNew);
  free(treeIndex);

  Rprintf("\n");

  return; 
}


// extern end for .C interface
}



    /*********************************************** DEBUG *************************************/ 
   /* 
    for( j = 0; j < yRow; j++) {
      if( j == 0) {
        Rprintf("\n y: yRow = %d, xCol = %d \n\t", (int) yRow, (int) xCol );
        for( k = 0; k < xCol; k++) Rprintf("%d\t", (int) k); 
        Rprintf("\n");
      }
      Rprintf("%d:\t", (int) j);
      for( k = 0; k < xCol; k++) Rprintf("%4.2f\t", y[ xCol * j + k]); 
      Rprintf("\n"); 
    }
    
    for( j = 0; j < yRow; j++) {
      if( j == 0) {
        Rprintf("\n distances: yRow = %d, nNeighbors = %d \n\t", (int) yRow, (int) nNeighbors );
        for( k = 0; k < nNeighbors; k++) Rprintf("%d\t", (int) k); 
        Rprintf("\n");
      }
      Rprintf("%d:\t", (int) j);
      for( k = 0; k < nNeighbors; k++) Rprintf("%4.2f\t", distances[ nNeighbors * j + k]); 
      Rprintf("\n"); 
    }
    
    for( j = 0; j < yRow; j++) {
      if( j == 0) {
        Rprintf("\n neighbors: yRow = %d, nNeighbors = %d\n\t", (int) yRow, (int) nNeighbors );
        for( k = 0; k < nNeighbors; k++) Rprintf("%d\t", (int) k); 
        Rprintf("\n");
      }
      Rprintf("%d:\t", (int) j);
      for( k = 0; k < nNeighbors; k++) Rprintf("%d\t", (int) neighbors[ nNeighbors * j + k]); 
      Rprintf("\n"); 
    }
    
    for( j = 0; j < yRow; j++) {
      if( j == 0) {
        Rprintf("\n distancesQuery: yRow = %d, 1 = %d \n\t", (int) yRow, (int) 1 );
        for( k = 0; k < 1; k++) Rprintf("%d\t", (int) k); 
        Rprintf("\n");
      }
      Rprintf("%d:\t", (int) j);
      for( k = 0; k < 1; k++) Rprintf("%4.2f\t", distancesQuery[ 1 * j + k]); 
      Rprintf("\n"); 
    }
    
    for( j = 0; j < yRow; j++) {
      if( j == 0) {
        Rprintf("\n neighborsQuery: yRow = %d, 1 = %d\n\t", (int) yRow, 1 );
        for( k = 0; k < 1; k++) Rprintf("%d\t", (int) k); 
        Rprintf("\n");
      }
      Rprintf("%d:\t", (int) j);
      for( k = 0; k < 1; k++) Rprintf("%d (%d)\t", (int) neighborsQuery[ 1 * j + k] , (int) treeIndex[ neighborsQuery[ 1 * j + k] ]); 
      Rprintf("\n"); 
    }
    
   
    Rprintf("\n treeIndex: length = %d\n ", (int) currentTreeIndexMax );
    for( j = 0; j < currentTreeIndexMax; j++) Rprintf("%d ", (int) treeIndex[j]); 
    Rprintf("\n"); 

    for( j = 0; j < yRow; j++) {
      if( j == 0) {
        Rprintf("\n prob: yRow = %d, 1 = %d \n\t", (int) yRow, (int) 1 );
        for( k = 0; k < 1; k++) Rprintf("%d\t", (int) k); 
        Rprintf("\n");
      }
      Rprintf("%d:\t", (int) j);
      for( k = 0; k < 1; k++) Rprintf("%4.2f\t", prob[ 1 * j + k]); 
      Rprintf("\n"); 
    }
    
    for( j = 0; j < xRow; j++) {
      if( j == 0) {
        Rprintf("\n assignment: xRow = %d, 1 = %d\n\t", (int) xRow, 1 );
        for( k = 0; k < 1; k++) Rprintf("%d\t", (int) k); 
        Rprintf("\n");
      }
      Rprintf("%d:\t", (int) j);
      for( k = 0; k < 1; k++) Rprintf("%d \t", (int) assignment[ 1 * j + k]);
      Rprintf("\n"); 
    }

    Rprintf("Print Queue\n");
    queueArrayPrint( queueArray, xRow); 
    */ 
    /********************************************** END DEBUG *******************************/

