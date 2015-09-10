.flannAlgorithmEnum <-
function( param ) {
  
  if( param == "FLANN_INDEX_LINEAR" ) return( 0 );
  if( param == "FLANN_INDEX_KDTREE" ) return( 1 );
  if( param == "FLANN_INDEX_KMEANS" ) return( 2 );
  if( param == "FLANN_INDEX_COMPOSITE" ) return( 3 );
  if( param == "FLANN_INDEX_KDTREE_SINGLE" ) return( 3 );  # I think this is an error in the documentation
  if( param == "FLANN_INDEX_SAVED" ) return( 254 );
  if( param == "FLANN_INDEX_AUTOTUNED" ) return( 255 );

  stop( sprintf("Invalid Algorithm = %s",param) )
}
.returnKernelEnum <-
function( param ) {
  if( param == "NONE" )    return( 0 );
  if( param == "UNIFORM" ) return( 1 );
  if( param == "NORMAL" )  return( 2 );
  if( param == "HYBRID" )  return( 3 );
}
