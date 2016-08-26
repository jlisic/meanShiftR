/* Copyright (c) 2015-2016  Jonathan Lisic 
 * Last edit: 16/08/24 - 09:20:44
 * License: GPL (>=2) 
 */  

#ifndef MEDIAN_HEADER

#define MEDIAN_HEADER

/* a quick select median program */
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#include <time.h>


void partitionIndex(double ** A, size_t * x1, size_t  * x2, const size_t n, const double p); 

double quantile_quickSelectIndex( double ** A, const size_t k, const size_t n ); 

#endif
