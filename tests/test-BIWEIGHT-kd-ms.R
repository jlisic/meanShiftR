library(meanShiftR)

# function to calculate the shift
mean_shift_test <- function( x,y  ) {
  result_num <- matrix( rep(0,length(y)),nrow=NROW(y),ncol=NCOL(y) )
  result_den <- result_num

  # convert to matrix
  if(is.null(ncol(x))) x <- matrix(x,ncol=1)  
  if(is.null(ncol(y))) y <- matrix(y,ncol=1)  

  for( j in 1:NROW(y) ) {
    x_target <- y[j,] 

    u <- (t(x) - x_target)/bandwidth
   
    u2 <- matrix(colSums(u^2),ncol=1)
   
    w <- (1-t(u)^2)^2

    for(p in 1:NCOL(y) ) { 
      ww <- w
      ww[,p] <- 1-t(u)[,p]^2
      if( j == 1 & p == 1) {
        ww0 <<- ww
        w0 <<- w
        wut <<- u
      }
      ww <-  apply(ww,1,prod)
      result_den[j,p] <- sum(ww*(u2<=1)) 
      result_num[j,p] <- sum(x[,p] * ww *(u2<=1)) 
    }
  
  }

  result <- result_num/result_den 

  # if the square difference is below 1e-08, then return original value
  no_change <- which(rowSums((result_num / result_den - x)^2) < 1e-08)
  result[no_change,] <- x[no_change,]
  result
}


# function to luster the results
mean_shift_cluster_test <- function( result ) {
  
  # convert to matrix if needed
  if(is.null(ncol(result))) result <- matrix(result,ncol=1)  
  
  # handle clustering
  cluster <- matrix(result[1,],ncol=NCOL(result))
  # initial cluster assignment
  assignment <- rep(1,NROW(result))

  print(cluster)

  for(j in 2:NROW(result)) { 

    # if within cluster distance (in order) make cluster assignment 
    min_dist <- colSums(( t(cluster) - result[j,] )^2 )/rowSums(cluster^2)
    if( j == 2) min_dist <<- min_dist 
       
    if(min(min_dist) < 1e-04) {
      assignment[j] <- which.min(min_dist) 
    } else {
      cluster <- rbind(cluster, result[j,])
      assignment[j] <- NROW(cluster)
    }
  }
  
  cluster[assignment,]
}



########################################
# Test 3a.  BIWEIGHT - EXACT 
# Univariate
########################################

## parameters
x <- c(1,2,3)
set.seed(10)
n <- 1000
p <-3 
x <- matrix(rnorm(n*p,mean=rep(c(0,3))), nrow=n, ncol=p)

bandwidth <- rep(3,p) 
iterations <-10


## run mean shift
runTime <- proc.time()
ms_result1 <- meanShift( 
  x,
  algorithm="KDTREE",
  kernelType="BIWEIGHT", 
  bandwidth=bandwidth, 
  nNeighbors=30,
  iterations=iterations
)
print((proc.time() - runTime)[3])



########################################
# Test 2b.  BIWEIGHT - EXACT 
# Multivariate 
########################################



## run mean shift
runTime <- proc.time()
ms_result2 <- meanShift( x,kernelType="BIWEIGHT", bandwidth=bandwidth, iterations=iterations)
print((proc.time() - runTime)[3])


## compare results
dif <- abs(ms_result1$value - ms_result2$value)


## print results
print(sprintf("max dif = %f",max(dif)))
print(sprintf("max dif index = %d",which.max(dif)))






