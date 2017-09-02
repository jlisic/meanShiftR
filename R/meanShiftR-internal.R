.algorithmEnum <- function( param ) {
  
  if( param == "LINEAR" )    return( 0 );
  if( param == "KDTREE" )    return( 1 );
  #if( param == "DUALTREE" )  return( 2 );
  #if( param == "LSH" )       return( 3 );
  #if( param == "FASTGAUSS" ) return( 4 );
  #if( param == "MERGETREE" ) return( 5 );

  stop( sprintf("Invalid Algorithm = %s",param) )
}

.kernelEnum <- function( param ) {
  if( param == "NORMAL" )       return(0);
  if( param == "EPANECHNIKOV")  return(1);
  if( param == "BIWEIGHT")      return(2);
  if( param == "TRIWEIGHT")     return(3);
  if( param == "TRICUBE")       return(4);
  
  stop( sprintf("Invalid Kernel = %s",param) )
}

