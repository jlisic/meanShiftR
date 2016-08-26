#' Mean shift classification
#' 
#' \code{meanShift} performs classification of a set of query points using
#' steepest ascent to local maxima in a kernel density estimate.
#'
#' @param queryData A matrix or vector of points to be classified by the mean
#'   shift algorithm.  Values must be finite and non-missing.
#' @param trainData A matrix or vector of points used to form a kernel density
#'   estimate.  The local maxima from this kernel density estimate will be used
#'   for steepest ascent classification.
#' @param nNeighbors A scalar indicating the number neighbors to consider for
#'   the kernel density estimate.  This is useful to speed up approximation by
#'   approximating the kernel density estimate.  The default is all data.
#' @param algorithm A string indicating the algorithm to use for nearest neighbor
#'   searches.  Currently, only "LINEAR" and "KDTREE" methods are supported.
#' @param bandwidth A vector of length equal to the number of columns in the 
#'   queryData matrix, or length one when queryData is a vector.  This value will be 
#'   used in the kernel density estimate for steepest ascent classification.  The 
#'   default is one for each dimension.
#' @param alpha A scalar tuning parameter for normal kernels.  When this paramter
#'   is set to zero, the mean shift algorithm will operate as usual.  When this
#'   paramter is set to one, the mean shift algorithm will be approximated through
#'   Newton's Method.  When set to a value between zero and one, a generalization
#'   of Newton's Method and mean shift will be used instead providing a means
#'   to balance convergence speed with stability.  The default is zero, mean shift.
#' @param iterations The number of iterations to perform mean shift.
#' @param epsilon A scalar used to determine when to terminate the iteration of a
#'   individual query point.  If the distance between the query point at iteration \code{i}
#'   and \code{i+1} is less than epsilon, then iteration ceases on this point.
#' @param epsilonCluster A scalar used to determine the minimum distance between distinct
#'   clusters.  This distance is applied after all iterations have finished and in 
#'   order of the rows of \code{queryData}.
#' @param parameters A scalar or vector of paramters used by the specific algorithm.
#'   There are no optional parameters for the "LINEAR" method, "KDTREE" supports an 
#'   optional parameter for the maximum number of points to store in a leaf node.
#' @return A list is returned containing two items: \code{assignment}, a vector of 
#'   classifications.  \code{value}, a vector or matrix containing the location of 
#'   the classified local maxima in the support, each row is associated with the 
#'   classified index in \code{assignment}.  
#'   
#' @examples 
#' x <- matrix(runif(20),10,2)
#' classification <- meanShift(x,x)
#' @useDynLib meanShiftR
#' @export


meanShift <-
function(
  queryData,                      
  trainData,                       
  nNeighbors = NROW(trainData),   
  algorithm = "LINEAR",
  bandwidth,
  alpha = 0.0,
  iterations=10,
  epsilon = 0.00000001,
  epsilonCluster = 0.0001,
  parameters = NULL
) {
  
  # handle missing bandwidth
  if( missing(bandwidth) ){
    if(is.null(ncol(trainData))) {
      bandwidth=1
    } else {
      bandwidth = rep(1,ncol(trainData))
    }
  }
  

  # get data size
  trainRow <- NROW(trainData)
  queryRow <- NROW(queryData)
  queryCol <- length(bandwidth) 

  # check bandwidth and columns
  if( length(trainData)/queryCol != trainRow) 
    stop(sprintf("Error: train rows do not match (%d != %d)", length(trainData)/queryCol, trainRow))
  if( length(queryData)/queryCol != queryRow) 
    stop(sprintf("Error: query rows do not match (%d != %d)", length(queryData)/queryCol, queryRow))
  
  # check if nNeighbors is set, note that if 
  nNeighbors = min(nNeighbors,NROW(trainData))

  # allocate space for the return vector for number of neighbors
  neighbors = rep(0,queryRow)
  distances = rep(0,queryRow*nNeighbors)

  # iterations check
  if( iterations < 1 ) stop("Error: iterations are not at least one.") 

  # convert algorithm ane kernel to enumerated classes
  kernelType = "NORMAL"
  kernelEnum <- .kernelEnum( kernelType ) 
  algorithmEnum <- .algorithmEnum( algorithm ) 


  # algorithm specific parameters
  # kdtree
  if( algorithmEnum == 1 ) {
    if( !is.null(parameters)) {
      intParameters = parameters 
    } else {
      intParameters = min( 40, trainRow ) 
    }
    dblParameters <- 0
  } else {
    intParameters <- 0
    dblParameters <- 0
  }

 
  # send our data to the C program
  r.result <- .C("R_meanShift",
    as.double(  t(queryData) ),                 # 1 data we query
    as.double(  t(trainData) ),                 # 2 sample we use to build index
    as.integer( neighbors ),                    # 3 
    as.integer( queryRow ),                     # 4 number of rows of data to query
    as.integer( trainRow ),                     # 5 number of rows of training data 
    as.integer( queryCol ),                     # 6 number of columns of data to query
    as.integer( nNeighbors ),                   # 7 number of neighbors
    as.integer( iterations ),                   # 8 number of iterations 
    as.double(  bandwidth ),                    # 9 method to form tree
    as.double(  alpha),                         # 10 alpha sequence
    as.double(  epsilon ),                      # 11 difference value to terminate 
    as.double(  epsilonCluster),                # 12 min cluster distance
    as.integer( kernelEnum ),                   # 13 kernel enumeration
    as.integer( algorithmEnum ),                # 14 algorithm enumeration
    as.integer( intParameters ),                # 15 integer parameters
    as.integer( dblParameters )                # 16 double parameters
  )
  
  #prob <- matrix(r.result[[3]],ncol=1)
  assignment <- matrix(r.result[[3]]+1,byrow=TRUE,ncol=1)
  value <- matrix(r.result[[1]],byrow=TRUE,ncol=queryCol) 

  return( list( assignment=assignment, value=value) )
  
}
