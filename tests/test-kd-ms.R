### todo 


library(meanShiftR)

# set a seed to make this reproducible 
set.seed(100)

# set the number of iterations to test
iter <- 10 
# set the number of points to simulate
n <- 10 

# set the bandwidth
h <- c(.5,.5)

# create example data
x1 <- matrix( rnorm( n ),ncol=2)
x2 <- matrix( rnorm( n ),ncol=2) + 2 
x <- rbind( x1, x2 ) 

########### meanShiftR kd ###################
run.time <- proc.time()
result <- meanShift(
  x, 
  x, 
  algorithm="KDTREE",
  bandwidth=h,
  alpha=0,
  nNeighbors=5,
  iterations = iter, 
  parameters=c(3,1000)
) 
meanShiftR_kd_runtime <- (proc.time()-run.time)[3]

# assignment
meanShiftR_kd_assignment <- result$assignment

# value
meanShiftR_kd_value <- result$value



## run an unaccelerated check
result <- x 

for( i in 1:iter ) {
  for( j in 1:nrow(x) ) {
  y <- result[j,]
    nNeighbors <- 5 
    
    # get nearest neighbors
    xx <- ((t(x) - y)*(1/h))^2
    d <- colSums(xx) 
    d.sort <- sort(d, index.return=T)
    
    d.exp <- exp( d.sort$x[1:nNeighbors] * -1/2 )   
    w = sum(d.exp)
    
    result[j,] <- colSums(x[d.sort$ix[1:nNeighbors],] * d.exp / w )
  }
}

# compare results  
print( max(abs(result - meanShiftR_kd_value )))


