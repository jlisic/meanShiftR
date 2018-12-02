# compare FNN to meanShiftR using exact nearest neighbors

library(FNN)
library(meanShiftR)


i <- 1698

set.seed(i)

n <- 2
m <- 2000 
p <- 3 
k <- 3 

x <- matrix(rnorm(p*n),ncol=p)
y <- matrix(rnorm(p*m),ncol=p)

labels <- rep(1,n)

a <- proc.time()
result_class <- FNN::knn(x,y,labels,k=k, algorithm="kd_tree")
print( (proc.time() - a)[3])



index_class <- attr(result_class,"nn.index")
dist_class  <- attr(result_class,"nn.dist")^2

a <- proc.time()
result_ms <- knn_meanShift(y,x,k=k )
print( (proc.time() - a)[3])



print(max(abs( result_ms$neighbors - index_class[,k:1])))
