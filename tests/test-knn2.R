# compare FNN to meanShiftR using exact nearest neighbors

library(meanShiftR)

i <- 1698

set.seed(i)

n <- 10 
p <- 3 
k <- 3 

type <- c(1,1,1,2)
x <- matrix( c( rnorm(p*n) , sample(1:3,size=n,replace=TRUE)), nrow=n)

a <- proc.time()
result_ms <- knn_meanShift(x,x,k=k)
print( (proc.time() - a)[3])


#calc dist

mixed_dist <- function( x, y, type ) { 
  sum((x[type == 1] - y[type == 1])^2) +
  sum((x[type != 1] != y[type != 1])) 
}

Dist <- matrix(0,n,n)
for( i in 1:n ) {
  for( j in 1:n ) {
    Dist[i,j] <- mixed_dist(x[i,], x[j,])
  }
}

Dist_rank <- t(apply(Dist,1,order))[,3:1]

print("Sum of rank differences:")
print(sum(abs(Dist_rank - result_ms$neighbors)))

Dist_k <- t(apply(Dist,1,sort))[,3:1]

print("Max abs difference in distance: ")
print(max(abs(Dist_k - result_ms$distance)))



