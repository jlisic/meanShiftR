
  x.sd = .50

  set.seed(100)
  n <- 40
  n <- 5 
  x1 <- matrix( c(
        rnorm( n,mean=-1, sd=x.sd) ,
        rnorm( n,mean=1, sd=x.sd)
        ),ncol=2)  
  x2 <- matrix( c(
        rnorm( n,mean=1, sd=x.sd) ,
        rnorm( n,mean=1, sd=x.sd)
        ),ncol=2)  
  x3 <- matrix( c(
        rnorm( n,mean=1, sd=x.sd) ,
        rnorm( n,mean=-1, sd=x.sd)
        ),ncol=2)  
  x4 <- matrix( c(
        rnorm( n,mean=-1, sd=x.sd) ,
        rnorm( n,mean=-1, sd=x.sd)
        ),ncol=2)  

  x <- rbind( x1, x2, x3, x4)


  #png( filename = "scatterPlotNNS.png")
  #plot(x, col=rep(1:4,each=n))
  #dev.off()

  

  # create a 2d density plot
#  library(ggplot2)
#  library(MASS)
#
#  x.df <- as.data.frame(x)
#  colnames(x.df) <- c('x','y')
#  x.df$Class <- as.factor(rep(1:4,each=n)) 
# 
#  m <- ggplot( x.df, aes(x=x, y=y, color=x.df$Class)) + geom_point()  + stat_density2d(aes(fill=..level..), geom="polygon") 
#  plot(m)
  
  





  library(meanShiftR)

bandwidth <- c(.75,.75)


################################################################
################################################################
iter <- 100
################################################################
################################################################
for( i in 1:20) {
  alpha <- 0

  nn <- meanShiftR::meanShift(
    x, 
    x, 
    nNeighbors=nrow(x),
    kernelMethod = "HYBRID", 
    bandwidth=bandwidth,
    alpha=alpha,
    epsilon = -1,
    clusterEpsilon = -1,
    iterations = iter,
    debugTrain=F
  )

  library(LPCM)
  nn.ms <- ms(
    x,
    h=bandwidth, 
    plotms=0,
    thr = -1,
    iter = iter,
    scaled=F
  )

a <- nn$value - nn.ms$cluster.center[nn.ms$cluster.label,]  
print( max(abs(a)))
}
