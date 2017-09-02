library(meanShiftR)
library(LPCM)

# set a seed to make this reproducible 
set.seed(100)

# set the number of iterations to test
iter <- 20

# set the number of points to simulate
n <- 20 
p <- 10 

# set the bandwidth
h <- rep(0.5, p)


# create example data
x1 <- matrix( rnorm( n/2 * p ),ncol=p)
x2 <- matrix( rnorm( n/2 * p ),ncol=p) + 2 
x <- rbind( x1, x2 ) 

########### meanShiftR ###################
run.time <- proc.time()
result <- meanShift(
  x, 
  x, 
  algorithm="KDTREE",
  bandwidth=h,
  alpha=0,
  nNeighbors=nrow(x),
  iterations = iter, 
  parameters=c(30,17)
) 
meanShiftR_kd_runtime <- (proc.time()-run.time)[3]

# assignment
meanShiftR_kd_assignment <- result$assignment

# value
meanShiftR_kd_value <- result$value


########### meanShiftR ###################
run.time <- proc.time()
result <- meanShift(
  x, 
  x, 
  bandwidth=h,
  alpha=0,
  iterations = iter 
) 
meanShiftR_runtime <- (proc.time()-run.time)[3]

# assignment
meanShiftR_assignment <- result$assignment

# value
meanShiftR_value <- result$value


########### LPCM ###################
runtime <- proc.time()
result <- ms(
            x,
            h=h, 
            scaled=FALSE, 
            iter=iter, 
            plotms=-1)
LPCM_runtime <- (proc.time()-runtime)[3]

# assignment
LPCM_assignment <- result$cluster.label

# value
LPCM_value <- result$cluster.center[LPCM_assignment,]

print(sprintf("Elapsed time meanShiftR kdtree = %f", meanShiftR_kd_runtime))
print(sprintf("Elapsed time meanShiftR = %f", meanShiftR_runtime))
print(sprintf("Elapsed time LPCM ms = %f", LPCM_runtime))


print( max(abs(meanShiftR_value - LPCM_value) ))
print( max(abs(meanShiftR_kd_value - LPCM_value) ))


print( sprintf("Number of differences between meanShiftR and LPCM = %d",
               sum(meanShiftR_assignment != LPCM_assignment)))
print( sprintf("Number of differences between meanShiftR kdtree and LPCM = %d",
               sum(meanShiftR_kd_assignment != LPCM_assignment))) 

