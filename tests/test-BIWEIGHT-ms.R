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

## paramters
x <- c(1,2,3)
set.seed(10)
n <- 100 
p <- 1 
x <- matrix(rnorm(n*p,mean=rep(c(0,3))), nrow=n, ncol=p)

bandwidth <- rep(3,p) 
iterations <-1


## run mean shift
ms_result <- meanShift( x,x,kernelType="BIWEIGHT", bandwidth=bandwidth, iterations=iterations)


## run test

#execute multiple shifts
result <- x
for( i in 1:iterations) {
  result <- mean_shift_test(x,result)
}

# cluster the results
result2 <- mean_shift_cluster_test(result)


## compare results
dif <- abs(result2 - ms_result$value)


## print results
print(sprintf("max dif = %f",max(dif)))
print(sprintf("max dif index = %d",which.max(dif)))
print( sprintf("check val %f",result2[which.max(dif)]))
print( sprintf("ms value %f",ms_result$value[which.max(dif)]))




########################################
# Test 2b.  BIWEIGHT - EXACT 
# Multivariate 
########################################


## parameters
x <- c(1,2,3)
set.seed(10)
n <- 10 
p <-3 
x <- matrix(rnorm(n*p,mean=rep(c(0,3))), nrow=n, ncol=p)

bandwidth <- rep(3,p) 
iterations <-1


## run mean shift
ms_result <- meanShift( x,x,kernelType="BIWEIGHT", bandwidth=bandwidth, iterations=iterations)


## run test

#execute multiple shifts
result <- x
for( i in 1:iterations) {
  result <- mean_shift_test(x,result)
}

# cluster the results
result2 <- mean_shift_cluster_test(result)


## compare results
dif <- abs(result2 - ms_result$value)


## print results
print(sprintf("max dif = %f",max(dif)))
print(sprintf("max dif index = %d",which.max(dif)))
print( sprintf("check val %f",result2[which.max(dif)]))
print( sprintf("ms value %f",ms_result$value[which.max(dif)]))






