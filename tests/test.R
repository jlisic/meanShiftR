# a comparison of timing between meanShiftR and LPCM's mean shift

library(meanShiftR)
library(LPCM)

# set a seed to make this reproducible 
set.seed(100)

# set the number of iterations to test
iter <- 100

# set the number of points to simulate
n <- 200 

# set the bandwidth
h <- c(0.5,0.5)

# create example data
x1 <- matrix( rnorm( n ),ncol=2)
x2 <- matrix( rnorm( n ),ncol=2) + 2 
x <- rbind( x1, x2 ) 

########### meanShiftR ###################
run.time <- proc.time()
result <- meanShift(
  x, 
  x, 
  algorithm="KDTREE",
  bandwidth=h,
  alpha=0,
  nNeighbors=200,
  iterations = iter, 
  parameters=c(30,0)
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


cat(sprintf("Elapsed time meanShiftR kdtree = %f\n", meanShiftR_kd_runtime))
cat(sprintf("Elapsed time meanShiftR = %f\n", meanShiftR_runtime))
cat(sprintf("Elapsed time LPCM ms = %f\n\n", LPCM_runtime))


cat("max diff:\n")

cat( sprintf( " mean shift, LPCM diff = %f\n", max(abs(meanShiftR_value - LPCM_value) )))
cat( sprintf( " mean shif_kdt, LPCM diff = %f\n", max(abs(meanShiftR_kd_value - LPCM_value) )))


cat( sprintf("Number of differences between meanShiftR and LPCM = %f\n",
               sum(meanShiftR_assignment != LPCM_assignment)/n))
cat( sprintf("Number of differences between meanShiftR kdtree and LPCM = %f\n",
               sum(meanShiftR_kd_assignment != LPCM_assignment)/n)) 
#
