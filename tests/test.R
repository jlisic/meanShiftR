library(meanShiftR)

set.seed(100)
n <- 100 
# example data
  x1 <- matrix( rnorm( n ),ncol=2)
  x2 <- matrix( rnorm( n ),ncol=2) + 2 
  x <- rbind( x1, x2 ) 

  nn <- meanShift(
    x, 
    x, 
    nNeighbors=nrow(x),
    kernelMethod = "HYBRID", 
    bandwidth=c(3,3),
    alpha=0,
    iterations = 100,
    debugTrain = T 
  ) 
  print(nn$assignment)


 

