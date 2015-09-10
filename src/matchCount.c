
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "R.h"
#include "Rinternals.h"
//#define CSTACK_DEFNS 7
#include "Rinterface.h"

#include "Rmath.h"
//#include "omp.h"


void rMatchCount( 
    int * pixel,           /* this is the raster image of groups */ 
    int * match,           /* this is the raster image of items within the groups */ 
    int * assign,
    int * count,
    int * mPtr,     /* max number of matches */
    int * nPtr
    ) {
  
  size_t m = *mPtr;

  size_t N = *nPtr;

  size_t i,j; /* iterator */
  int group;
  int item;
  size_t assignIndex;

  /* for each observation find another observation to replace it */
//  #pragma omp parallel private(j,group,item,assignIndex) 
  for( i =0; i < N; i ++){ 
    group = pixel[i];                                /* identify index */ 
    item  = match[i];

    if( group >= 0 ) 
      if( item >= 0 ) 
        for( j = 0; j < m; j++) {
         
          assignIndex = group * m +j;               
         
          if( assign[assignIndex] < 0) {
            
//            #pragma omp critical 
            assign[assignIndex] = item; 
            
//            #pragma omp atomic
            count[assignIndex]++; 

            break;
          }
          if( assign[assignIndex] == item) {
//            #pragma omp atomic
            count[assignIndex]++;

            break;
          } 

        }
  }

  return;
}


