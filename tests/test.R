
# example data
  x1 <- matrix( rnorm( 10 ),ncol=2)
  x2 <- matrix( rnorm( 10 ),ncol=2) + 2 
  x <- rbind( x1, x2 ) 

  nn <- getNN(
    x, 
    x, 
    nNeighbors=3,
    kernelMethod = "HYBRID", 
    bandwidth=c(3,3),
    alpha=0,
    iterations = 100,
    debugTrain = T 
  ) 
  print(nn$assignment)


 

