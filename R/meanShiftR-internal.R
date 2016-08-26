.algorithmEnum <-
function( param ) {
  
  if( param == "LINEAR" )    return( 0 );
  if( param == "KDTREE" )    return( 1 );
  #if( param == "DUALTREE" )  return( 2 );
  #if( param == "LSH" )       return( 3 );
  #if( param == "FASTGAUSS" ) return( 4 );
  #if( param == "MERGETREE" ) return( 5 );

  stop( sprintf("Invalid Algorithm = %s",param) )
}
.kernelEnum <-
function( param ) {
  #if( param == "UNIFORM" ) return( 0 );
  #if( param == "TRIANGULAR" ) return( 1 );
  if( param == "NORMAL" )   return( 2 );
  #if( param == "EPANECHNIKOV") return(3);
  #if( param == "BIWEIGHT")     return(4);
  #if( param == "TRIWEIGHT")    return(5);
  #if( param == "TRICUBE")      return(6);
  #if( param == "COSINE")       return(7);
  #if( param == "LOGISTIC")     return(8);
  #if( param == "SIGMOID")      return(9);
  
  stop( sprintf("Invalid Kernel = %s",param) )
}


