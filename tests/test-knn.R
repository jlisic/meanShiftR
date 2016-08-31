library(FNN)
library(meanShiftR)


i <- 1698

#for( i in 1:1 ) {
set.seed(i)

n <- 2000
m <- 2000 
p <- 3 
k <- 200 

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
##if(max(abs( result_ms$neighbors - rev(index_class)))) break
##}
#
#if ( F ) {
#d <-  colSums((t(x) - c(y))^2)
#d.sort <- sort(d, index.return=T)
#d.sort.index <- sort(d, index.return=T)$ix
#d.sort <- sort(d)
#}
#

 #(c( 0.040882,  0.967041, 0.367679) - c(-1.444927, -1.578433, -1.554843))^2
 #sum((c( 0.040882,  0.967041, 0.367679) - c(-1.444927, -1.578433, -1.554843))^2)

