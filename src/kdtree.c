/* Copyright (c) 2015-2016  Jonathan Lisic 
 * Last edit: 16/08/24 - 09:20:44
 * License: GPL (>=2) 
 */  

#include "kdtree.h"
#include "median.h"

/* printf fixing */
#ifdef CLI 
 #define PRINTF printf
 #define RUNIF (double) rand() / RAND_MAX  
#else 
 #include "R.h"
 #include "Rmath.h"
 #define PRINTF Rprintf
 #define RUNIF runif(0.0,1.0)
#endif

#ifdef DDEBUG 
  #define KNNLDEBUG PRINTF("To the left %f < %f = dist (%p)\n", distMin, *dist, (void *) c), 
  #define KNNRDEBUG PRINTF("To the right %f < %f = dist (%p)\n", distMin, *dist, (void*) c),
  //#define KNNNODEDEBUG PRINTF("Node: %p\n", (void *) c); 
  //#define KNNLISTDEBUG for( i = 0; i < k; i ++) PRINTF(" %d : %f \n", (int) indexes[i], dist[i]);
  #define KNNNODEDEBUG
  #define KNNLISTDEBUG
#else 
  #define KNNLDEBUG 
  #define KNNRDEBUG 
  #define KNNNODEDEBUG
  #define KNNLISTDEBUG
#endif



// function to create a new Tree 
rootNodePtr createTree( size_t K, size_t leafSize, size_t n, double * data ) {
  
  rootNodePtr y = malloc( sizeof( rootNode ) );
  
  y->pointerIndex = calloc( n, sizeof( size_t * ) );

  // setup root node 
  y->K = K;
  y->leafSize = leafSize;
  y->root = NULL;
  y->data = data;
  y->n = n;
  y->type = NULL;

  return(y);
}




// function to create a new Tree 
void deleteTree( rootNodePtr r ) {
 
  if(r->pointerIndex != NULL) free( r->pointerIndex ); 
  r->pointerIndex = NULL;
  r->data = NULL;
  r->type = NULL;
    
  deleteNode( r, r->root ); 

  free(r);
  return;
}




// add an element to the tree 
nodePtr buildIndex( 
    rootNodePtr r,      // root pointer 
    size_t dim,         // current dim
    size_t m,           // current length of obs
    size_t * indexPtr   // pointer to obs indexes 
  ) {
 
  size_t i,K; 
  size_t * indexLeftPtr = NULL;
  size_t * indexRightPtr = NULL;
  size_t indexLeftSize;
  size_t indexRightSize;

  nodePtr c = createNode(r);
  c->indexUsed = m;
  c->index = indexPtr;
  c->dim = dim;

  K = r->K;
   
  // do we have too many points? 
  if( m <= r->leafSize ) {

    // save the final pointer locations 
    for( i = 0; i < m; i++) 
      r->pointerIndex[ indexPtr[i] ] = &( indexPtr[i] );

    return c;
  } 

  // if we are here we have too many points 
  // create children
  // figure out our new dim
  // split data and give to children 
  c-> split = splitData( 
    r->data,
    c->index, 
    &indexLeftPtr,
    &indexRightPtr,
    &indexLeftSize,
    &indexRightSize,
    m, 
    K,
    dim  
    ); 

  free(c->index);
  c->index = NULL; 

  // move current contents to new children
  c->left  = buildIndex( r, (dim+1) % K, indexLeftSize , indexLeftPtr);
  c->right = buildIndex( r, (dim+1) % K, indexRightSize, indexRightPtr);

  return c;
}




// function to create a simple node 
nodePtr createNode( rootNodePtr r ) {

  nodePtr c = malloc( sizeof( node ) );

  c->index = NULL;
  c->left  = NULL;
  c->right = NULL;
  c->indexUsed = 0;
  c->split = 0;
  c->dim = 0;

  return(c);
}




// function to delete a node 
void deleteNode( rootNodePtr r, nodePtr c ) {

  if( c == NULL ) return;

  if( c->index != NULL) free(c->index);

  // to the left
  if( c->left != NULL ) {
    deleteNode( r, c->left );
    c->left = NULL;
  }

  // to the right
  if( c->right != NULL ) {
    deleteNode( r, c->right );
    c->right = NULL;
  }

  free(c);
  c = NULL;
  return;
}




// split and create children
double splitData( 
    double * y,
    size_t * index, 
    size_t ** indexLeft,
    size_t ** indexRight,
    size_t * indexLeftSize,
    size_t * indexRightSize,
    size_t n, 
    size_t p,
    size_t dim 
    ) {

  double split;
  size_t splitIndex,i;

  // get the median 
  double * x =  calloc( n, sizeof(double) );        //allocate some temporary space for finding the median
  double ** xPtr =  calloc( n, sizeof( double * ) ); //allocate some temporary space for finding the median
  
  // create input for qsort
  for( i = 0; i < n; i++) {
    x[i] = y[ index[i] * p + dim];
    xPtr[i] = &(x[i]);
  }
  // use quick sort to find the median 
  // get split
  splitIndex = n/2;
  
  split = quantile_quickSelectIndex( xPtr, splitIndex, n ); 
    
  *indexLeftSize  = n / 2;
  *indexRightSize = n - *indexLeftSize;
  
  *indexLeft  = calloc(*indexLeftSize , sizeof(size_t) );
  *indexRight = calloc(*indexRightSize , sizeof(size_t) );
    
  // now let's have some fun with pointer math 
  for( i = 0; i < *indexLeftSize ; i++) (*indexLeft)[i]  = index[xPtr[i] - x];
  for( i = 0; i < *indexRightSize; i++) (*indexRight)[i] = index[ xPtr[*indexLeftSize + i] -x ];
  
  free(xPtr);
  free(x); 

  return split;
}




// a function to print the tree 
void printTree( rootNodePtr r, nodePtr c ) {

  size_t i;

  PRINTF("node: %p\n", (void *) c);
  if( c->index != NULL) {
    for( i=0; i < c->indexUsed; i++) PRINTF("%d ", (int) c->index[i]); 
  } 
  PRINTF("\n\tleft: %p right %p (split = %f) \n", (void *) c->left, (void*) c->right, c->split );
  PRINTF("\n");

  if( c->left ) {
    PRINTF("left ");
    printTree( r, c->left);
  }

  if( c->right ) {
    PRINTF("right ");
    printTree( r, c->right);
  }

}

// function to find the minimal Euclidian distance ( with handling for avoiding cat unpack)
void getClosest2( 
    rootNodePtr r, 
    nodePtr c, 
    size_t k,
    double * queryPoint,
    size_t * indexes,
    double * dist, 
    double * weight,
    double * tieBreak 
  ) {

  size_t i,j,l,d;

  size_t K = r->K;
  double * x = r->data;            
  double currentDist;
  double tmp;

  // iterate over all obs on the leaf 
  for( i = 0; i < c->indexUsed; i++) {

    // get first index 
    j = c->index[i]; 

//printf("  j = %d:  ", (int) j);

    // calculate L2 distance
    for( d=0, currentDist=0; d < K; d++) { 
      if( r->type[d] == 1 ) {
        tmp = x[j * K + d] - queryPoint[d];
        tmp *= tmp * weight[d];
      } else {
        if( x[j * K + d] != queryPoint[d] ) {
          tmp = weight[d];
        }
        else {
          tmp = 0;
        }
      }
//printf(" %f, ", tmp);
      currentDist += tmp;
    }

//printf(" tot= %f\n", currentDist);

    
    // if smaller than prior index update 
    // dist is ordered from biggest dist to smallest
    if( currentDist < dist[0] ) {

      for(l=1;l<k;l++) { 
        // if the distance is bigger than current dist we will shift things down.
        if( currentDist < dist[l] ) {
          dist[l-1] = dist[l];
          indexes[l-1] = indexes[l];
        }
        else break;
      }
   
      dist[l-1] = currentDist;
      indexes[l-1] = j;
          
    } 
    /* 
    // if it is a tie
    else if( currentDist == *dist ) {
    
      // generate a deviate on (0,1) and pick the biggest
      // each obs has the same prob of being largest due
      // to exchangability  
      newTieBreak = RUNIF;

      if( *tieBreak < 0 ) *tieBreak = RUNIF;  // if no tie was set

      if( newTieBreak > *tieBreak) *tieBreak = newTieBreak;
      closestIndex = i;
    }
    */
  }

  return;
}



// function to find the minimal Euclidian distance 
void getClosest( 
    rootNodePtr r, 
    nodePtr c, 
    size_t k,
    double * queryPoint,
    size_t * indexes,
    double * dist, 
    double * weight,
    double * tieBreak 
  ) {

  size_t i,j,l,d;

  size_t K = r->K;
  double * x = r->data;            
  double currentDist;
  double tmp;

  // iterate over all obs on the leaf 
  for( i = 0; i < c->indexUsed; i++) {

    // get first index 
    j = c->index[i]; 

//    printf("  j = %d\n", (int) j);

    // calculate L2 distance
    for( d=0, currentDist=0; d < K; d++) { 
      tmp = x[j * K + d] - queryPoint[d];
      tmp *= tmp * weight[d];
      currentDist += tmp;
    }
    
    // if smaller than prior index update 
    // dist is ordered from biggest dist to smallest
    if( currentDist < dist[0] ) {

      for(l=1;l<k;l++) { 
        // if the distance is bigger than current dist we will shift things down.
        if( currentDist < dist[l] ) {
          dist[l-1] = dist[l];
          indexes[l-1] = indexes[l];
        }
        else break;
      }
   
      dist[l-1] = currentDist;
      indexes[l-1] = j;
          
    } 
    /* 
    // if it is a tie
    else if( currentDist == *dist ) {
    
      // generate a deviate on (0,1) and pick the biggest
      // each obs has the same prob of being largest due
      // to exchangability  
      newTieBreak = RUNIF;

      if( *tieBreak < 0 ) *tieBreak = RUNIF;  // if no tie was set

      if( newTieBreak > *tieBreak) *tieBreak = newTieBreak;
      closestIndex = i;
    }
    */
  }

  return;
}



// find the k nearest neighbors with type checks
void find_knn2( 
    rootNodePtr r, 
    nodePtr c, 
    size_t k,            // number of items to find
    double * queryPoint, // location of point in support 
    size_t * indexes,    // k indexes
    double * dist,       // distances
    double medianDist,
    double * weight,
    double * tieBreak   // tie break
  ) {

  double distMin;
  
  // return if c == NULL 
  if( c == NULL ) {
    PRINTF(" not good\n");
    return;
  }

  // is there anything here ? 
  // if there is we get the closest item 
  if( c->index != NULL ) { 
    getClosest2(r,c,k,queryPoint,indexes,dist,weight,tieBreak); 
    // debug macros 
    KNNNODEDEBUG KNNLISTDEBUG
    return;
  }
   
  // as we iterate through the tree we calculate the distance between the current split element and the query point 
  //       x_1
  //        |
  //        |    . y_1
  //        |
  //        |
  //        |
  //
  //       
  //        |   y_2 
  //        |    .
  //   -----+----- x_2
  //        |
  //        |
  //

  // this is a dist calculation
  if( r->type[ c->dim ] == 1 ) {
    distMin = (queryPoint[c->dim] - c->split);
    distMin *= distMin * weight[c->dim];
  } else {
    if( queryPoint[c->dim] != c->split ) {
      distMin = weight[c->dim];
    } else {
      distMin = 0;
    }
  }


  if( distMin < medianDist ) medianDist = distMin;


  // first check if the query point is less than split 
  if( queryPoint[c->dim] <= c->split ) {
    
      if(medianDist < *dist) {
        KNNLDEBUG find_knn2( r, c->left, k, queryPoint, indexes, dist, medianDist, weight, tieBreak );  
      } else 
      {
      //  printf("%f > %f, dim=%d, query =%f, split=%f FAIL!\n", distMin, *dist, (int) c->dim, queryPoint[c->dim], c->split ); 
      } 
    
      // now check if there is a point in the split that can be close 
      if(medianDist < *dist) {
        KNNRDEBUG find_knn2( r, c->right, k, queryPoint, indexes, dist, medianDist, weight, tieBreak ); 
      } 
      else {
      //  printf("%f > %f, FAIL!\n", distMin, *dist); 
      } 
      
  } else { // the query point is greater than the split   

      if(medianDist < *dist) {
        KNNRDEBUG find_knn2( r, c->right, k, queryPoint, indexes, dist, medianDist, weight, tieBreak );  
      } 
      else {
      //  printf("%f > %f, FAIL!\n", distMin, *dist); 
      } 
      
      // now check if there is a point in the split that can be close 
      if(medianDist < *dist) {
        KNNLDEBUG find_knn2( r, c->left, k, queryPoint, indexes, dist, medianDist, weight, tieBreak );  
      } 
      else {
      //  printf("%f > %f, FAIL!\n", distMin, *dist); 
      } 
  }

  return;
}




// find the k nearest neighbors 
void find_knn( 
    rootNodePtr r, 
    nodePtr c, 
    size_t k,            // number of items to find
    double * queryPoint, // location of point in support 
    size_t * indexes,    // k indexes
    double * dist,       // distances
    double medianDist,
    double * weight,
    double * tieBreak   // tie break
  ) {

  double distMin;
  
  // return if c == NULL 
  if( c == NULL ) {
    PRINTF(" not good\n");
    return;
  }

  // is there anything here ? 
  // if there is we get the closest item 
  if( c->index != NULL ) { 
    getClosest(r,c,k,queryPoint,indexes,dist,weight,tieBreak); 
    // debug macros 
    KNNNODEDEBUG KNNLISTDEBUG
    return;
  }
   
  // as we iterate through the tree we calculate the distance between the current split element and the query point 
  //       x_1
  //        |
  //        |    . y_1
  //        |
  //        |
  //        |
  //
  //       
  //        |   y_2 
  //        |    .
  //   -----+----- x_2
  //        |
  //        |
  //

  // this is a dist calculation
  distMin = (queryPoint[c->dim] - c->split);
  distMin *= distMin * weight[c->dim];


  if( distMin < medianDist ) medianDist = distMin;


  // first check if the query point is less than split 
  if( queryPoint[c->dim] <= c->split ) {
    
      if(medianDist < *dist) {
        KNNLDEBUG find_knn( r, c->left, k, queryPoint, indexes, dist, medianDist, weight, tieBreak );  
      } else 
      {
      //  printf("%f > %f, dim=%d, query =%f, split=%f FAIL!\n", distMin, *dist, (int) c->dim, queryPoint[c->dim], c->split ); 
      } 
    
      // now check if there is a point in the split that can be close 
      if(medianDist < *dist) {
        KNNRDEBUG find_knn( r, c->right, k, queryPoint, indexes, dist, medianDist, weight, tieBreak ); 
      } 
      else {
      //  printf("%f > %f, FAIL!\n", distMin, *dist); 
      } 
      
  } else { // the query point is greater than the split   

      if(medianDist < *dist) {
        KNNRDEBUG find_knn( r, c->right, k, queryPoint, indexes, dist, medianDist, weight, tieBreak );  
      } 
      else {
      //  printf("%f > %f, FAIL!\n", distMin, *dist); 
      } 
      
      // now check if there is a point in the split that can be close 
      if(medianDist < *dist) {
        KNNLDEBUG find_knn( r, c->left, k, queryPoint, indexes, dist, medianDist, weight, tieBreak );  
      } 
      else {
      //  printf("%f > %f, FAIL!\n", distMin, *dist); 
      } 
  }

  return;
}



/* test code */
#ifdef CLI 
int main () {

/*
    Distances for x:
               0           1         2          3         4          5          6
    0  0.0000000 0.602768283 0.3622356 0.31157833 0.2268138 0.51773743 0.20586722
   [1] 0.6027683 0.000000000 0.4998974 0.15602813 0.5304343 0.03717622 0.13046443
    2  0.3622356 0.499897351 0.0000000 0.66348365 0.0189031 0.64578733 0.43482101
   [3] 0.3115783 0.156028132 0.6634836 0.00000000 0.5856961 0.05270530 0.02438501
    4  0.2268138 0.530434288 0.0189031 0.58569605 0.0000000 0.62905469 0.37149263
   [5] 0.5177374 0.037176221 0.6457873 0.05270530 0.6290547 0.00000000 0.07089969
   [6] 0.2058672 0.130464432 0.4348210 0.02438501 0.3714926 0.07089969 0.00000000
   [7] 0.5823510 0.003165375 0.4320187 0.17973374 0.4686796 0.05702175 0.13663419
    8  0.0913264 0.261416675 0.4087474 0.06750114 0.3126583 0.17907200 0.02503165
    9  0.1487531 0.402522596 0.7035774 0.07134651 0.5619861 0.24253339 0.08122322

    For finding neighbors in the restricted case of index 0 and the left node
         1,      7,      5,      3,     6 
    0.6028, 0.5824, 0.5178, 0.3116, 0.2059 

    and the right node
         
         2      4      9       8      0
    0.3622 0.2268 0.1488  0.0913 0.0000

    for the top 3 knn

         9     8      0
    0.1488 0.913 0.0000
*/

  double x[20] = {
    0.68512126, 0.3251399, 
    0.05171296, 0.7740967,  
    0.74261974, 0.9242472, 
    0.13036790, 0.3870030,  
    0.77980495, 0.7918827, 
    0.01413735, 0.5849822,  
    0.25770368, 0.4773944, 
    0.09543018, 0.8095111,  
    0.39014922, 0.3908506,  
    0.32050716, 0.1994035   
  };


  size_t ncol=2;
  size_t nrow=10;
  size_t i;
  size_t j;
  double * queryPoint;
  double dist;
  double tieBreak = -1;
  size_t k = 3;

  rootNodePtr myTree = NULL;

  myTree = createTree(ncol, 5, nrow, x);

  // this will get freed
  size_t * index = calloc(nrow, sizeof(size_t));
  for(i = 0; i < nrow; i++) index[i] = i;

  // k index
  size_t * kIndex = calloc(k, sizeof(size_t));
  // k dist 
  double * kDist = calloc(k, sizeof(double));
  // k dist 
  double * weight = calloc(k, sizeof(double));


  myTree->root = buildIndex( 
    myTree,      // root pointer 
    0,           // current dim
    nrow,        // current length of obs
    index        // pointer to obs indexes 
  ); 

  printTree( myTree, myTree->root );
  
  /********************************************************************************/ 
  /* This is a test for k nearest neighbors applied at the leaf node              */ 
  /********************************************************************************/ 
 
 

  // print out data set 
  for( i = 0; i < k; i ++) printf(" %d : %f \n", (int) kIndex[i], kDist[i]);


  // Setup Test Data Set 
  for(i = 0; i < k; i++) kIndex[i] = myTree->n;
  for(i = 0; i < k; i++) kDist[i] = INFINITY; 
  for(i = 0; i < k; i++) weight[i] = 1.0; 


  // test of getting the k closest units on the left node 
  getClosest( 
    myTree, 
    myTree->root->left, 
    k,
    x,
    kIndex,
    kDist, 
    weight,
    &tieBreak 
  ); 

  for( i = 0; i < k; i ++) printf(" %d : %f \n", (int) kIndex[i], kDist[i]);

  
  // Setup Test Data Set on the right node 
  for(i = 0; i < k; i++) kIndex[i] = myTree->n;
  for(i = 0; i < k; i++) kDist[i] = INFINITY; 
  for(i = 0; i < k; i++) weight[i] = 1.0; 

  getClosest( 
    myTree, 
    myTree->root->right, 
    k,
    x,
    kIndex,
    kDist, 
    weight,
    &tieBreak 
  ); 

  for( i = 0; i < k; i ++) printf(" %d : %f \n", (int) kIndex[i], kDist[i]);
  
  
  /********************************************************************************/ 
  /* This is a test for k nearest neighbors applied at the root node              */ 
  /********************************************************************************/ 
 
  // run test for k = 3

  // Setup Test Data Set 
  printf("the full knn k=3:\n"); 

  for(i = 0; i < k; i++) kIndex[i] = myTree->n;
  for(i = 0; i < k; i++) kDist[i] = INFINITY; 
  for(i = 0; i < k; i++) weight[i] = 1.0; 

  printf("final results:\n");
  find_knn( myTree, myTree->root, k, x, kIndex, kDist, INFINITY, weight, &tieBreak);

  for( i = 0; i < k; i ++) printf(" %d : %f \n", (int) kIndex[i], kDist[i]);

  // run test for k = 1 
  k=1;
  
  printf("the full knn k=1:\n"); 

  // Setup Test Data Set 
  for(i = 0; i < k; i++) kIndex[i] = myTree->n;
  for(i = 0; i < k; i++) kDist[i] = INFINITY; 
  for(i = 0; i < k; i++) weight[i] = 1.0; 

  find_knn( myTree, myTree->root, k, x, kIndex, kDist, INFINITY, weight, &tieBreak);
  printf("final results:\n");
  for( i = 0; i < k; i ++) printf(" %d : %f \n", (int) kIndex[i], kDist[i]);



  deleteTree( myTree );
  free(kIndex);
  free(kDist);

  return 0;
}

#else 

  /********************************************************************************/ 
  /* An R interface for KNN                                                       */ 
  /********************************************************************************/ 


void R_knn( 
  double * queryPoints,  // point to query for
  double * x,           // data to reference for the query
  int * type, // new
  int * xnrowPtr,
  int * nrowPtr,        // number of rows
  int * ncolPtr,        // number of columns
  double * kDist,       // distance vector
  int * indexInt,       // length of index
  int * kPtr,           // number of nearest neighbors
  double * weight,
  int * leafSizePtr,     // leaf size
  double * maxDist
 ) {

  size_t i;
  size_t j;
  size_t typeEval=0;
  size_t d;
  double * dist;
  double * queryPoint; 
  double tieBreak = -1;
  size_t * kIndex = NULL;

  size_t k        = (size_t) * kPtr;
  size_t ncol     = (size_t) * ncolPtr;
  size_t nrow     = (size_t) * nrowPtr;
  size_t xnrow     = (size_t) * xnrowPtr;
  size_t leafSize = (size_t) * leafSizePtr;


  rootNodePtr myTree = NULL;

  myTree = createTree(ncol, leafSize, nrow, x);

  // add on type
  myTree->type = type;

  // speed things up by checking if there is a need to evaluate any categorical data

  for( d=0; d < ncol; d++) {
    if( myTree->type[d] == 1 ) typeEval++;
  }

  // index for k-d tree
  // this will get freed
  size_t * index = calloc(nrow, sizeof(size_t));
  for(i = 0; i < nrow; i++) index[i] = i;

  // k index
  // this will get freed
  
  // build the index
  myTree->root = buildIndex( 
    myTree,      // root pointer 
    0,           // current dim
    nrow,        // current length of obs
    index        // pointer to obs indexes 
  ); 


  if( typeEval == ncol ) {

  #pragma omp parallel private(i,j,kIndex,tieBreak,queryPoint,dist)
  {

    kIndex = calloc(k, sizeof(size_t));

    #pragma omp for
    for(i=0; i < xnrow; i++) { 

//      printf("i = %d\n", (int) i );
      for(j = 0; j < k; j++) kIndex[j] = myTree->n;
  
      queryPoint = queryPoints + i*ncol; 
      dist = kDist + i*k; 
      // query tree
      find_knn( myTree, myTree->root, k, queryPoint, kIndex, dist, *maxDist, weight, &tieBreak);


      // copy data over 
      for( j = 0; j < k; j++) indexInt[i*k +j] = 1 + (int) kIndex[j]; 
    }

    free(kIndex);
  }

  } else {

  #pragma omp parallel private(i,j,kIndex,tieBreak,queryPoint,dist)
  {

    kIndex = calloc(k, sizeof(size_t));

    #pragma omp for
    for(i=0; i < xnrow; i++) { 

//      printf("i = %d\n", (int) i );
      for(j = 0; j < k; j++) kIndex[j] = myTree->n;
  
      queryPoint = queryPoints + i*ncol; 
      dist = kDist + i*k; 
      // query tree
      find_knn2( myTree, myTree->root, k, queryPoint, kIndex, dist, *maxDist, weight, &tieBreak);


      // copy data over 
      for( j = 0; j < k; j++) indexInt[i*k +j] = 1 + (int) kIndex[j]; 
    }

    free(kIndex);
  }


  }

  // clean up
  deleteTree( myTree );

  return;
}

#endif



