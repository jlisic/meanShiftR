/* Copyright (c) 2019  Jonathan Lisic 
 * last editied 19/02/25 - 10:08:19 
 * License: GPL (>=2) 
 */  

#include "dfToSparse.h" 


void R_dfToSparse(
 int *x ,         // 1 input from original matrix (catgories recorded as integers)
 int *n,          // 2 number of rows
 int *p,          // 3 number of columns
 int *new_i,      // 4 column indexes
 int *new_p,      // 5 pointers for sparse matrix
 int *max_cat,    // 6 number of categories per column
 int *max_cat_cumsum , // 7 cumulative sume of the number of categories per column
 int *type              // 8 Not used yet
 ) {

// *** original test r code ***
//  for( i in 1:n ) {
//  a_row <- a[i,]
//  
//    for( j in 1:p ) {
//      if( a_row[j] != 0 ) {
//        new_i = c(new_i, max_cat_cumsum[j] + (a_row[j] -1) )
//      }
//    }
//  }

  size_t i,j,k=0;

  //int * x_row = (int *) malloc( sizeof(int) * (*p) ) ;

  for( i=0; i < *n; i++) {

    for(j=0; j<*p; j++) {
      //printf("%d, ", (int) x[i * (*p) + j] ); 
      new_i[k]  = max_cat_cumsum[j] + (x[j * (*n) +i] - 1) ;
      k++; 
    }
    //printf("\n");

  }

  //free( x_row );
  
  return;
}



