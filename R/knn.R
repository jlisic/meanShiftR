#' K-d tree based k nearest neighbor search 
#' 
#' \code{knn_meanShift} performs a search for the k nearest neighbors of a single 
#' point, where nearest is determined by the Mahalanobis distance.  This search
#' is performed through a k-d tree.
#'
#' @param point A vector to find k nearest neighbors for. 
#' @param trainData A matrix or vector of potential nearest neighbors.
#' @param k A scalar indicating the number neighbors to find.
#' @param weight A scalar or vector of length equal to the number of columns of 
#'   \code{trainData}.   This value is used as the diagonal elements for the 
#'   inverse covariance matrix of the Mahalanobis distance.
#' @param leafSize A scalar used to specify the number of points to store in the 
#'   leaf nodes. 
#' @return A list is returned containing two items: \code{neighbors}, a vector of 
#'   k indexes corresponding to the nearest neighbors in  \code{trainData}.  \code{value}, a 
#'   vector of scalars containing the k  distances between the neighbors found in 
#'   \code{trainData} and \code{point}.
#'   
#' @examples 
#' x <- matrix(runif(20),10,2)
#' neighbors <- knn_meanShift(c(0,0),x)
#' @useDynLib meanShiftR
#' @export
knn_meanShift <-
function(
  point,                      
  trainData,
  k=min(5,NROW(trainData)),
  weight,
  leafSize=40  
) {
  

  # get data size
  trainRow <- NROW(trainData)
  trainCol <- length(point) 

  # check bandwidth and columns
  if( length(trainData)/trainCol != trainRow) stop("Error: incorrection dimensions of point and trainData.")
  
  # check if k is set, note that if 
  k = min(k,NROW(trainData))

  if( missing(weight) ) weight <- rep(1,trainCol)
  if( length(weight) != trainCol ) stop("Error: weight length is not equal to the number of columns.")

  # allocate space for the return vector for number of neighbors
  neighbors = rep(-1,k)
  distances = rep(Inf,k)
 
  # send our data to the C program
  r.result <- .C("R_knn",
    as.double(  point ),                 # 1 data we query
    as.double(  t(trainData) ),          # 2 set to search for nn 
    as.integer( trainRow ),             # 3 number of rows of data to search 
    as.integer( trainCol ),              # 4 number of columns of data to search 
    as.double( distances ),              # 5 number of neighbors
    as.integer( neighbors ),             # 6 k neighbors 
    as.integer(  k ),                    # 7 k 
    as.double( weight ),
    as.integer( leafSize ),              # 8 number of nodes in the leaf 
    NAOK=T
  )

  return( list( neighbors=r.result[[6]], distances=r.result[[5]] ) )
  
}
