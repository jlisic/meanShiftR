# compare FNN to meanShiftR using exact nearest neighbors

library(FNN)
library(meanShiftR)


i <- 1698

set.seed(i)

n <- 10000 
m <- 20000
p <- 100 
k <- 5 

x <- matrix(exp(rnorm(p*n)),ncol=p)
y <- matrix(exp(rnorm(p*m)),ncol=p)


zero_rate <- 0.9 

# add some zeros
x[sample( 1:length(x), size=round(p*n * zero_rate))] <- 0
y[sample( 1:length(y), size=round(p*m * zero_rate))] <- 0



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


library(Matrix)
x_sparse <- Matrix(x,sparse=TRUE)
y_sparse <- Matrix(y,sparse=TRUE)

#print( "y_sparse@i")
#print( as.integer(y_sparse@i) )
#print( "y_sparse@p")
#print( as.integer(y_sparse@p) )
#print( "y_sparse@i")
#print( as.integer(x_sparse@i) )
#print( "x_sparse@p")
#print( as.integer(x_sparse@p) )

#a <- proc.time()
#print('dense check')
#result_sparse_ms <- knn_meanShift(y_sparse,x_sparse,k=k,leafSize=4 )
#print( (proc.time() - a)[3])


tx_sparse <- t(x_sparse)
ty_sparse <- t(y_sparse)
#print(x_sparse)
#print(y_sparse)


a <- proc.time()
print('dense check transpose')
result_sparse_ms <- knn_meanShift(ty_sparse,tx_sparse,k=k,transpose=1 )
print( (proc.time() - a)[3])





print("distances")
print(max(abs(result_sparse_ms$distances - result_ms$distances)))
print("neighbors")
print(max(abs(result_sparse_ms$neighbors - result_ms$neighbors)))


#print(result_sparse_ms$distances - result_ms$distances)
#
#
#check_dist <- function(i,j,y,x) {
#  sum((y[i,] - x[j,])^2)
#}
#
#print( check_dist(17,7,y,x)) 
#print( check_dist(17,10,y,x)) 
#print( check_dist(17,8,y,x)) 
#
#print( check_dist(17,1,y,x)) 
#print( check_dist(17,7,y,x)) 
#print( check_dist(17,6,y,x)) 
#
#
#print( check_dist(1,4,y,x))
#print( check_dist(1,8,y,x))
#print( check_dist(1,2,y,x))


