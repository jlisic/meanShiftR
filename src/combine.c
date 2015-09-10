#include <stdio.h>
#include <stdlib.h>

// this is a function to combine data


void combineData( 
    double * x, 
    int * rowPtr, 
    int * colPtr, 
    double * epsilon, 
    int * assign
  ) {

  double sum,z;
  double epsilon2 = (*epsilon) * (*epsilon);
  int i,j,k,jIndex;
  int row = *rowPtr;
  int col = *colPtr;
  int * group  = calloc(row,sizeof(int));


  /* we don't start at 0 because it is already assigned via calloc default*/
  for(i=1; i < row; i++) {

    /* the goal here is to see if a new vector is close to a prior vector */
    /* if they are sufficiently close we cluster them together */
    for( j=0; j < row; j++) {
      if( group[j] == 0 ) {
        if( j != 0) {
          group[j] = i;
          assign[i] = i;
          break;
        }
      }

      /* get the squared l2 distance between two vectors */
      sum = 0;
      jIndex = group[j];

      for( k=0; k < col; k++) {
        z = x[ i*col + k ] - x[ jIndex * col + k ];
        sum += z*z; 
      }

      /* check if the squared distance is less than epsilon */
      if (sum <= epsilon2) {
        assign[i] = jIndex;
        break;
      }
    }
  }

  free(group);
  return;
}







