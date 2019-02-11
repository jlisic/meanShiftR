#' K-d tree based k nearest neighbor search 
#' 
#' \code{knn_meanShift} performs a search for the k nearest neighbors of a single 
#' point, where nearest is determined by the Mahalanobis distance.  This search
#' is performed through a k-d tree.
#'
#' @param points n vectors stored in an n by p matrix.  k nearest neighbors are
#'   found for each vector.
#' @param trainData A matrix or vector of potential nearest neighbors.
#' @param k A scalar indicating the number neighbors to find.
#' @param weight A scalar or vector of length equal to the number of columns of 
#'   \code{trainData}.   This value is used as the diagonal elements for the 
#'   inverse covariance matrix of the Mahalanobis distance.
#' @param leafSize A scalar used to specify the number of points to store in the 
#'   leaf nodes. 
#' @param maxDist A vector specifying the maximum value of the Mahalanobis that
#'   will be considered.
#' @param transpose A boolean determining if you should query by columns instead of
#'   the default of rows (only for sparse matricies).
#' @return A list is returned containing two items: \code{neighbors}, an n by k
#'   matrix of k indexes for each of the n vectors in \code{points}, corresponding to 
#'   the nearest neighbors in \code{trainData}.  \code{value}, a matrix of scalars 
#'   containing the k distances between the neighbors found in \code{trainData} 
#'   and \code{points}.
#'   
#' @examples 
#' x <- matrix(runif(20),10,2)
#' neighbors <- knn_meanShift(c(0,0),x)
#' @useDynLib meanShiftR, .registration = TRUE
#' @export
knn_meanShift <-
function(
  evalData,                      
  trainData,
  k=min(5,NROW(trainData)),
  weight,
  leafSize=40,
  maxDist=Inf,
  transpose = FALSE
) {
  

  # get data size
  if( !transpose ) {
    evalDataRow <- NROW(evalData)
    trainRow <- NROW(trainData)
    trainCol <- NCOL(trainData) 
  } else {
    evalDataRow <- NCOL(evalData)
    trainRow <- NCOL(trainData)
    trainCol <- NROW(trainData)
  }

  # check if k is set, note that if 
  k = min(k,trainRow)

  if( missing(weight) ) weight <- rep(1,trainCol)
  if( length(weight) != trainCol ) stop("Error: weight length is not equal to the number of columns.")

  # allocate space for the return vector for number of neighbors
  neighbors = rep(-1,k*evalDataRow)
  distances = rep(Inf,k*evalDataRow)

 print(class(trainData) )
  # send our data to the C program
  if( !("dgCMatrix" %in% class(trainData) )) {
    #print('start dense')
    # dense matrix handling
    tmp <- -1
    r.result <- .C("R_knn",
      as.double(  t(evalData) ),             # 1 data we query
      as.double(  t(trainData) ),          # 2 set to search for nn 
      as.integer( evalDataRow ),             # 3 number of rows of evalData 
      as.integer( trainRow ),              # 4 number of rows of data to search 
      as.integer( trainCol ),              # 5 number of columns of data to search 
      as.double( distances ),              # 6 distances 
      as.integer( neighbors ),             # 7 k neighbors 
      as.integer(  k ),                    # 8 k 
      as.double( weight ),                 # 9 weights for dist
      as.integer( leafSize ),              # 10 number of nodes in the leaf 
      as.double(maxDist),                  # 11 max acceptable distance
      NAOK=T
    )
  return( list( 
               neighbors=matrix(r.result[[7]],ncol=k,byrow=T) , 
               distances=matrix(r.result[[6]],ncol=k,byrow=T) 
               ) )
  } else {
    #print('start sparse')
    # sparse matrix handling
    r.result <- .C("R_knn_sparse",
      as.double(  evalData@x ),             # 1 data we query
      as.integer( evalData@i ),             # 2 sparse index 
      as.integer( evalData@p ),             # 3 sparse pointers 
      as.integer( length(evalData@i) ),     # 4 length of sparse index 
      as.integer( length(evalData@p) ),     # 5 length of sparse pointers 
      as.double(  trainData@x ),          # 6 set to search for nn 
      as.integer( trainData@i ),          # 7 sparse index 
      as.integer( trainData@p ),          # 8 sparse pointers 
      as.integer( length(trainData@i) ),  # 9 length of sparse index 
      as.integer( length(trainData@p) ),  # 10 length of sparse pointers 
      as.integer( evalDataRow ),            # 11 number of rows of evalData 
      as.integer( trainRow ),             # 12 number of rows of data to search 
      as.integer( trainCol ),             # 13 number of columns of data to search 
      as.double( distances ),             # 14 distances 
      as.integer( neighbors ),            # 15 k neighbors 
      as.integer(  k ),                   # 16 k 
      as.double( weight ),                # 17 weights for dist
      as.integer( leafSize ),             # 18 number of nodes in the leaf 
      as.double(maxDist),                 # 19 max acceptable distance
      as.integer( transpose),             # 20 transpose indicator
      NAOK=T
    )

    return( list( 
               neighbors=matrix(r.result[[15]],ncol=k,byrow=T) , 
               distances=matrix(r.result[[14]],ncol=k,byrow=T) 
               ) )
  }
  
}
