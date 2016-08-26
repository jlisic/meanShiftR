library(meanShiftR)
library(LPCM)

# set a seed to make this reproducible 
set.seed(100)

# set the number of iterations to test
iter <- 30 

# set the number of points to simulate
n <- 10 

# set the bandwidth
h <- rep(3,2)
h <- c(2,2)

# create example data
  x1 <- matrix( rnorm( n ),ncol=2)
  x2 <- matrix( rnorm( n ),ncol=2) + 2 
  x <- rbind( x1, x2 ) 

# run the linear meanShiftR mean shift algorithm
run.time <- proc.time()
nn.meanShift_kd <- meanShift(
  x, 
  x, 
  bandwidth=h,
  algorithm="KDTREE",
  nNeighbors=5,
  alpha=0,
  iterations = iter,
  parameters=3
) 
run.time.meanShift_kd <- (proc.time()-run.time)[3]

# get the meanShiftR classification and values
nn.meanShift.assignment_kd <- nn.meanShift_kd$assignment
nn.meanShift.value_kd <- nn.meanShift_kd$value

# run the linear meanShiftR mean shift algorithm
run.time <- proc.time()
nn.meanShift <- meanShift(
  x, 
  x, 
  bandwidth=h,
  alpha=0,
  iterations = iter 
) 
run.time.meanShift <- (proc.time()-run.time)[3]

# get the meanShiftR classification and values
nn.meanShift.assignment <- nn.meanShift$assignment
nn.meanShift.value <- nn.meanShift$value


run.time <- proc.time()
nn.ms <- ms(
            x,
            h=h, 
            scaled=FALSE, 
            iter=iter, 
            plotms=-1)
run.time.ms <- (proc.time()-run.time)[3]

print(sprintf("Elapsed time meanShiftR kdtree = %f", run.time.meanShift_kd))
print(sprintf("Elapsed time meanShiftR = %f", run.time.meanShift))
print(sprintf("Elapsed time LPCM ms = %f", run.time.ms))

# get the assignment and values from LPCM's mean shift implementation
nn.ms.assignment <- nn.ms$cluster.label
nn.ms <- nn.ms$cluster.center[nn.ms$cluster.label,]

print("max diff")
print( max(abs(nn.meanShift.value - nn.ms) ))
print( max(abs(nn.meanShift.assignment - nn.ms.assignment)) )
print( max(abs(nn.meanShift.value_kd - nn.meanShift.value) ))
print("min diff")
print( min(abs(nn.meanShift.value - nn.ms) ))
print( min(abs(nn.meanShift.assignment - nn.ms.assignment)) )
print( min(abs(nn.meanShift.value_kd - nn.meanShift.value) ))
print("mean diff")
print( mean(abs(nn.meanShift.value - nn.ms) ))
print( mean(abs(nn.meanShift.assignment - nn.ms.assignment)) )
print( mean(abs(nn.meanShift.value_kd - nn.meanShift.value) ))




result <- x 

for( i in 1:iter ) {
  for( j in 1:nrow(x) ) {
  y <- result[j,]
    nNeighbors <-5
    
    # get nearest neighbors
    xx <- ((t(x) - y)*(1/h))^2
    d <- colSums(xx) 
    d.sort <- sort(d, index.return=T)
    
    d.exp <- exp( d.sort$x[1:nNeighbors] * -1/2 )   
    w = sum(d.exp)
    
    result[j,] <- colSums(x[d.sort$ix[1:nNeighbors],] * d.exp / w )
  }
}

print( max(abs(result - nn.ms) ))
print( max(abs(result - nn.meanShift.value_kd) ))


