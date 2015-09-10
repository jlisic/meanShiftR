sampleClassify <-
function( 
  x, 
  maxIter=1, 
  sampled=FALSE, 
  nNeighbors=500, 
  bandwidth = NULL, 
  kernelMethod = "HYBRID",
  epsilon= 0.01, 
  spatial=F, 
  alpha=0, 
  samplingRate=11 ) {

 
  # grid sampling
  if( is.character(sampled)) {
    if( toupper(sampled) == 'GRID' ) {
      if( samplingRate > max(nrow(x),ncol(x))) stop( "can't sample rate exceeds observations")
      r.row <- nrow(x)
      r.col <- ncol(x)
  
      # random starts
      k.row <- sample(1:r.row,size=1) 
      k.col <- sample(1:r.col,size=1) 
      
      row.values <- sort(unique(c( seq( from=k.row, to=1, by=-1*samplingRate), seq( from=k.row, to=r.row, by=samplingRate) )))
      col.values <- sort(unique(c( seq( from=k.col, to=1, by=-1*samplingRate), seq( from=k.col, to=r.col, by=samplingRate) )))
  
      # create sample
      sampled <- sapply( col.values, function(x) x + r.col * (row.values -1) )
    
    }
  }
 
  x.values <- values(x)
  if( spatial == T ) x.values <- cbind(xyFromCell(x,1:NROW(x.values)), x.values)
 
  # create data set 
  x.nona <- !is.na(rowSums(x.values))
  x.values.nona <- x.values[x.nona,]
 
  # create training set 
  if( spatial == F ) {
    train <- matrix(x.values,ncol=nlayers(x))[sampled,,drop=FALSE]
  } else {
    train <- matrix(x.values,ncol=nlayers(x)+2)[sampled,,drop=FALSE]
  }
  train <- train[rowSums(is.na(train))==0,,drop=FALSE]

  if( length(bandwidth) != ncol(train) ) stop("layers of image and bandwidth do not match")


  nNeighbors = min(nNeighbors,nrow(train))

  nn <- meanShift(
    train, 
    train,
    nNeighbors,
    kernelMethod = "HYBRID", 
    bandwidth=bandwidth,
    alpha=alpha,
    iterations = maxIter, 
    epsilon = 0.5,
    interpolate = x.values.nona,
    debugTrain = F
  ) 

  # copy over values
  value <- x.values
  value[x.nona,] <- nn$interpolateValue

  # copy over the assignment
  x <- subset(x,1)
  x.values <- values(x)
  x.values[!x.nona] <- NA
  x.values[x.nona] <- nn$interpolateAssignment 
  values(x) <- x.values

  return(list(
         assignment=x, 
         value=value,
         train = train, 
         sampleTrain = sampled,
         assignmentTrain=nn$assignment,
         valueTrain=nn$value,
         probabilityTrain=nn$probability
       ) )
}
