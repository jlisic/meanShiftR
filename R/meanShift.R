meanShift <-
function(
  queryData,                      # data to mode seek on
  trainData,                      # data used for training 
  nNeighbors = nrow(trainData),   # number of neighbors
  indexParameters = NULL,
  searchParameters = NULL,
  kernelMethod = "NONE",
  bandwidth = NULL,
  alpha = 0.1,
  iterations=1,
  epsilon = 0.0001,
  clusterEpsilon = 0.0001,
  interpolate,
  debugTrain=F
) {
  

  # get data size
  trainRow <- nrow(trainData)
  queryRow <- nrow(queryData)
  queryCol <- ncol(queryData)

  if( is.null(queryRow) ) {
    # check if the data is sane  
    if( length(queryData) != queryCol ) { stop("Error: queryData columns do not equal trainData columns") } 
    queryData <- matrix(queryData,nrow=1)
  }

  nNeighbors = min(nNeighbors,nrow(trainData))

  # allocate space for the return vector for number of neighbors
  neighbors = rep(0,queryRow*nNeighbors)
  distances = rep(0,queryRow*nNeighbors)

  # iterations check
  if( iterations < 1 ) { stop("Error: iterations are not at least 1") }

  # create prob
  prob <- rep(-1,queryRow)
  neighborsQuery = rep(-1,queryRow)
  distancesQuery = rep(-1,queryRow)
  

  # default search parameters
  integerSearchParameters <- 0
  doubleSearchParameters <- 0

  # default index parameters defined by algorithm
  #if ( algorithm == "FLANN_INDEX_KDTREE" ) {
  integerAlgorithmParameters <- 4 
  doubleAlgorithmParameters <- 0
  #}
  paramKernelMethod <- .returnKernelEnum( kernelMethod ) 
 
  # handle bandwidth for newton and hybrid 
  if( kernelMethod %in% c("NEWTON", "HYBRID")  ) {
    if( is.null(bandwidth) ) stop("bandwidth NULL for kernel method")
    if( length(bandwidth) != queryCol ) stop("bandwidth not equal to columns x")
#    if( length(bandwidth) > 1 ) {
#      queryData <- queryData %*% diag(1/bandwidth) # apply bandwidth for l2 dist
#      trainData <- trainData %*% diag(1/bandwidth) # apply bandwidth for l2 dist
#    } else {
#      queryData <- queryData / bandwidth # apply bandwidth for l2 dist
#      trainData <- trainData / bandwidth # apply bandwidth for l2 dist
#    }
  }


  # handle interpolation
  if( !missing(interpolate) ) {
    print('interpolateFound')
    if( kernelMethod %in% c("NEWTON", "HYBRID")  ) {
      if( is.null(bandwidth) ) stop("bandwidth NULL for kernel method")
      if( length(bandwidth) != queryCol ) stop("bandwidth not equal to columns x")
      if( length(bandwidth) > 1 ) {
        #interpolate <- interpolate %*% diag(1/bandwidth) # apply bandwidth for l2 dist
        interpolateIndex <- 1:nrow(interpolate)
  
      } else {
        #interpolate <- interpolate / bandwidth # apply bandwidth for l2 dist
        interpolateIndex <- 1:nrow(interpolate)
      }
    }
  } else {
    # if there is nothing to interpolate set index to -1
    interpolate <- -1
    interpolateIndex <- -1
  }
 

  # record all 
  #   assignments             assignmentsDebug
  #   weights from neighbors  weightsDebug
  #   neighbors               neighborsDebug
  #   values                  valuesDebug
  if( debugTrain ) {

    print("Warning Debugging is being performed, this may exhaust system memory for large data sets")
    assignmentsDebug <- as.integer(rep(1, NROW(trainData)*iterations ) )   # set as 1 to indicate debug in C program
    weightsDebug     <- as.double( rep(0, NROW(trainData)*iterations*nNeighbors ) )
    neighborsDebug   <- as.integer(rep(0, NROW(trainData)*iterations*nNeighbors ) ) 
    valuesDebug      <- as.double( rep(0, NROW(trainData)*iterations*ncol(trainData) ) ) 
  } else {
    assignmentsDebug <- as.integer(rep(0, 1))
    weightsDebug     <- as.double( rep(0, 1))
    neighborsDebug   <- as.integer(rep(0, 1))
    valuesDebug      <- as.double( rep(0, 1))
  }

  Cprog <- proc.time()

  # send our data to the C program
  r.result <- .C("R_meanShiftNN",
    as.double(  t(queryData) ),                 # 1 data we query
    as.double(  t(trainData) ),                 # 2 sample we use to build index
    as.integer( neighbors ),                   # 3
    as.double(  distances ),                   # 4
    as.double(  prob ),                         # 5
    as.integer( neighborsQuery ),              # 6 
    as.double(  distancesQuery ),              # 7 
    as.integer( queryRow ),                    # 8 number of rows of data to query
    as.integer( trainRow ),                    # 9 number of rows of training data 
    as.integer( queryCol ),                    # 10 number of columns of data to query
    as.integer( nNeighbors ),                  # 11 number of neighbors
    as.integer( integerSearchParameters ),     # 12 integer Parameters for searching
    as.double(  doubleSearchParameters ),       # 13 double Parameters  for searching
    as.integer( integerAlgorithmParameters ),  # 14 integer Parameters for algorithm, first parameter is for the number of trees 
    as.double(  doubleAlgorithmParameters ),    # 15 double Parameters for algorithm
    as.integer( paramKernelMethod ),           # 16 method to form tree
    as.double(  bandwidth ),                    # 17 method to form tree
    as.double(  alpha),                         # 18 alpha sequence
    as.integer( iterations ),                  # 19 number of iterations
    as.double(  epsilon ),                      # 20 number of iterations
    as.double(  t(interpolate) ),               # 21 values to interpolate for
    as.integer( interpolateIndex ),             # 22 indexes from the interpolation
    as.integer( nrow(interpolate) ),             # 23 rows in interpolation matrix
    assignmentsDebug,                           # 24 assignment debug
    weightsDebug,                               # 25 weights debug
    neighborsDebug,                             # 26 neighbors debug
    valuesDebug,                                # 27 values debug
    as.double(  clusterEpsilon )                # 28 cluster cut off
  )
  
  print("C running time")
  print(proc.time() - Cprog) 

  prob <- matrix(r.result[[7]],ncol=1)
  assignment <- matrix(r.result[[6]],byrow=TRUE,ncol=1)
  value <- matrix(r.result[[1]],byrow=TRUE,ncol=queryCol) 

  # interpolate check
  if( interpolateIndex[1] >= 0 ) {
    interpolate      <- matrix(r.result[[21]],byrow=TRUE,ncol=queryCol) 
    interpolateIndex <- matrix(r.result[[22]],byrow=TRUE,ncol=1) 
  
    if( !debugTrain) {
      return( list( 
                   assignment=assignment, 
                   value=value, 
                   probability=prob,
                   interpolateAssignment=interpolateIndex, 
                   interpolateValue=interpolate
                ) 
      )
    }
  }

  # if debug
  if( debugTrain ) {
    assignmentsDebug  <- cbind(rep( 1:iterations, each=trainRow ),  matrix(r.result[[24]],byrow=TRUE,ncol=1) ) 
    weightsDebug  <- cbind(rep( 1:iterations, each=trainRow ),  matrix(r.result[[25]],byrow=TRUE,ncol=nNeighbors)  )
    neighborsDebug  <- cbind(rep( 1:iterations, each=trainRow ),  matrix(r.result[[26]],byrow=TRUE,ncol=nNeighbors) )
    valuesDebug  <- cbind(rep( 1:iterations, each=trainRow ),  matrix(r.result[[27]],byrow=TRUE,ncol=queryCol) )
 
   # debug and interpolate    
   if( interpolateIndex[1] >= 0 ) {
     return( list( 
                   assignment=assignment, 
                   value=value, 
                   probability=prob,
                   interpolateAssignment=interpolateIndex, 
                   interpolateValue=interpolate,
                   assignmentsDebug= assignmentsDebug,
                   weightsDebug= weightsDebug,
                   neighborsDebug= neighborsDebug,
                   valuesDebug= valuesDebug
                ) 
      )
    } else {
     return( list( 
                   assignment=assignment, 
                   value=value, 
                   probability=prob,
                   assignmentsDebug= assignmentsDebug,
                   weightsDebug= weightsDebug,
                   neighborsDebug= neighborsDebug,
                   valuesDebug= valuesDebug
                ) 
      )


    }
  }
    
  return( list( assignment=assignment, value=value, probability=prob ) )
  
}
