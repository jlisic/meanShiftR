/* Copyright (c) 2015-2016  Jonathan Lisic 
 * Last edit: 16/08/24 - 09:20:44
 * License: GPL (>=2) 
 */  

#ifndef KDTREE_HEADER

#define KDTREE_HEADER


#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>


/* create a node */
struct node {
  size_t dim;         /* dim split across */
  size_t * index;     /* index rows for x */
  size_t indexUsed;   /* since we are using arrays, this indicates the number of elements of each array */
  double split;       /* point split on for dim */
  struct node * left; /* values less than or equal to split */
  struct node * right; /* value greater than or equal to split */ 
};

/* create typedef */
typedef struct node node;
typedef struct node * nodePtr;


/* tree node */
struct rootNode {
  size_t K;                 /* dims of x */
  size_t leafSize;          /* max number of elements on each leaf */
  size_t n;                 /* number of rows in x */ 
  size_t ** pointerIndex;   /* pointer index to allow easy changing of returned values */
  double * data;            /* pointer to x */
  double * query_data;            /* pointer to x */
  int * sparse_i;
  int * sparse_p;
  int * sparse_i_n;
  int * sparse_p_n;
  int * query_sparse_i;
  int * query_sparse_p;
  int * query_sparse_i_n;
  int * query_sparse_p_n;
  size_t transpose;
  nodePtr root;             /* root node pointer */

};

/* create typedef */
typedef struct rootNode rootNode;
typedef struct rootNode * rootNodePtr;


double sparse_point(
    size_t i, //row
    size_t j, //col
    int * index,  // sparse matrix index
    int * pointer,  // sparse matrix pointer
    double * x, // sparse matrix value
    size_t transpose // transpose indicator
    ); 

double sparse_l2(
    size_t i,       //row
    int * i_index,    // sparse matrix index
    int * i_pointer,  // sparse matrix pointer
    double * i_x,     // sparse matrix value
    size_t j,      //col
    int * j_index,    // sparse matrix index
    int * j_pointer,  // sparse matrix pointer
    double * j_x,     // sparse matrix value
    double * weight,
    size_t p,      // NCOL
    double * i_vector,
    double * j_vector,
    double * w_vector
    ); 

/* function to print a tree */
void printTree( rootNodePtr r, nodePtr c ); 

/* function to create a new Tree */
rootNodePtr createTree( 
    size_t K,        // number of columns of data
    size_t leafSize, // maximum size of leaf nodes
    size_t n,        // number of rows of data
    double * data    // pointer to the data
  );

/* delete tree */
void deleteTree( rootNodePtr r );

/* build index for the tree */
nodePtr buildIndex( 
    rootNodePtr r,    // root node pointer
    size_t dim,       // dim to base next node on
    size_t m,         // size of index
    size_t * indexPtr // index array to add to new node
  ); 

/* build index for the tree */
nodePtr buildIndex_sparse( 
    rootNodePtr r,    // root node pointer
    size_t dim,       // dim to base next node on
    size_t m,         // size of index
    size_t * indexPtr // index array to add to new node
  ); 

// create Node 
nodePtr createNode( rootNodePtr r );

// delete a node and all of it's children 
void deleteNode( rootNodePtr r, nodePtr c ); 

// split and create children
double splitData( 
    double * y,
    size_t * index, 
    size_t ** indexLeft,
    size_t ** indexRight,
    size_t * indexLeftSize,
    size_t * indexRighSize,
    size_t n, 
    size_t p,
    size_t dim 
    ); 

// split and create children
double splitData_sparse( 
    rootNodePtr r, 
    size_t * index, 
    size_t ** indexLeft,
    size_t ** indexRight,
    size_t * indexLeftSize,
    size_t * indexRighSize,
    size_t n, 
    size_t p,
    size_t dim 
    ); 


// function to get the closest neighbor, with tie handling 
void getClosest_sparse( 
    rootNodePtr r, 
    nodePtr c, 
    size_t k,
    size_t queryPoint_index,
    size_t * indexes,
    double * dist, 
    double * weight,
    double * tieBreak
  ); 

// function to get the closest neighbor, with tie handling 
void getClosest( 
    rootNodePtr r, 
    nodePtr c, 
    size_t k,
    double * queryPoint,
    size_t * indexes,
    double * dist, 
    double * weight,
    double * tieBreak 
  );

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
    double * tieBreak    // tie break
  ); 

// find the k nearest neighbors 
void find_knn_sparse( 
    rootNodePtr r, 
    nodePtr c, 
    size_t k,            // number of items to find
    size_t queryPoint_index, // index of point in support 
    size_t * indexes,    // k indexes
    double * dist,       // distances
    double medianDist,
    double * weight,
    double * tieBreak   // tie break
  );



// R Knn interface (sparse)
void R_knn_sparse( 
  double * queryPoints,             // 1 - point to query for
  int * queryPoints_sparse_i,       // 2 - row index in sparse format 
  int * queryPoints_sparse_p,       // 3 - pointer to sparse observations
  int * queryPoints_sparse_i_n,     // 4 - length of i and values in sparse format
  int * queryPoints_sparse_p_n,     // 5 - length of p in sparse format
  double * x,                       // 6 - data to reference for the query
  int * x_sparse_i,                 // 7 - row index in sparse format 
  int * x_sparse_p,                 // 8 - pointer to sparse observations
  int * x_sparse_i_n,               // 9 - length of i and values in sparse format
  int * x_sparse_p_n,               // 10 - length of p in sparse format
  int * xnrowPtr,                   // 11 - 
  int * nrowPtr,                    // 12 - number of rows
  int * ncolPtr,                    // 13 - number of columns
  double * kDist,                   // 14 - distance vector
  int * indexInt,                   // 15 - length of index
  int * kPtr,                       // 16 - number of nearest neighbors
  double * weight,                  // 17 -
  int * leafSizePtr,                // 18 - leaf size
  double * maxDist,                 // 19 - 
  int * transpose                   // 20 transpose indicator
  );


// R Knn interface
void R_knn( 
  double * queryPoints,  // point to query for
  double * x,           // data to reference for the query
  int * xnrowPtr,
  int * nrowPtr,        // number of rows
  int * ncolPtr,        // number of columns
  double * kDist,       // distance vector
  int * indexInt,       // length of index
  int * kPtr,           // number of nearest neighbors
  double * weight,
  int * leafSizePtr,     // leaf size
  double * maxDist
 ); 

#endif
