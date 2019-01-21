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
#' @param type allows the user to specify if the value provided is to be treated
#'   as categorical without requiring one-hot encoding.
#' @param leafSize A scalar used to specify the number of points to store in the 
#'   leaf nodes. 
#' @param maxDist A vector specifying the maximum value of the Mahalanobis that
#'   will be considered.
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
  points,                      
  trainData,
  k=min(5,NROW(trainData)),
  weight,
  type,
  leafSize=40,
  maxDist=Inf
) {
  

  # get data size
  pointsRow <- NROW(points)
  trainRow <- NROW(trainData)
  trainCol <- NCOL(trainData) 

  if( missing(type) ) type <- rep(1,trainCol)
  if( length(type) != trainCol ) stop("Error: type length is not equal to the number of columns.")

  
  # check if k is set, note that if 
  k = min(k,NROW(trainData))

  if( missing(weight) ) weight <- rep(1,trainCol)
  if( length(weight) != trainCol ) stop("Error: weight length is not equal to the number of columns.")

  # allocate space for the return vector for number of neighbors
  neighbors = rep(-1,k*pointsRow)
  distances = rep(Inf,k*pointsRow)
 
  # send our data to the C program
  r.result <- .C("R_knn",
    as.double(  t(points) ),             # 1 data we query
    as.double(  t(trainData) ),          # 2 set to search for nn 
    as.integer( type ),                  # 3 type for distance calculation
    as.integer( pointsRow ),             # 4 number of rows of points 
    as.integer( trainRow ),              # 5 number of rows of data to search 
    as.integer( trainCol ),              # 6 number of columns of data to search 
    as.double( distances ),              # 7 number of neighbors
    as.integer( neighbors ),             # 8 k neighbors 
    as.integer(  k ),                    # 9 k 
    as.double( weight ),                 # 10 weights for dist
    as.integer( leafSize ),              # 11 number of nodes in the leaf 
    as.double(maxDist),                  # 12 max acceptable distance
    NAOK=T
  )

  return( list( 
               neighbors=matrix(r.result[[8]],ncol=k,byrow=T) , 
               distances=matrix(r.result[[7]],ncol=k,byrow=T) 
               ) )
  
}
