testShiftNN  <- FALSE 
testShiftNN12 <- FALSE 
testShiftNN2 <- FALSE 
testShiftNN3 <- FALSE 
testShiftNN4 <- FALSE 

if( "testShiftNN" %in% ls() ) {
  if( testShiftNN ) {

    # example data
    x <- matrix( c( 
               1, 1,
               1, 2,
               2, 1, 
               5, 3,
               5, 4,
               4, 5 ), byrow=T, ncol=2)


  #x1 <- matrix( rnorm( 10 ),ncol=2)
  #x2 <- matrix( rnorm( 10 ),ncol=2) + 2 
  #x <- rbind( x1, x2 ) 

  nn <- getNN(
    x, 
    x, 
    #nrow(x),
    nNeighbors=3,
    kernelMethod = "HYBRID", 
    bandwidth=c(3,3),
    #bandwidth=c(3,3),
    alpha=0,
    iterations = 100,
    debugTrain = T 
  ) 
  print(nn$assignment)


  # neighbor restriction check (only 3 neighbors)
  w <- apply( cbind( dnorm( (x[1,1] - x[1:3,1])/3 ), dnorm( (x[1,2] - x[1:3,2])/3 )), 1, prod)
  x1 <- sum(x[1:3,1] * w)/sum(w)
  x2 <- sum(x[1:3,2] * w)/sum(w)
  print(x1)
  print(x2)
  
  # no neighbor restriction check 
  w <- apply( cbind( dnorm( (x[1,1] - x[,1])/3 ), dnorm( (x[1,2] - x[,2])/3 )), 1, prod)
  x1 <- sum(x[,1] * w)/sum(w)
  x2 <- sum(x[,2] * w)/sum(w)
  print(x1)
  print(x2)
 

  }
}

if( "testShiftNN12" %in% ls() ) {
  if( testShiftNN12 ) {

    set.seed(400)

    # load raster
    library(raster)
#    r <- raster("~/src/flannR/data/test.gri")
#
#    Ext <- c( 
#    755774,
#    755960,
#    2064572,
#    2064731
#    )
#
#  r <- crop(r, extent(Ext) )  
 

  years <- c(1,2,3,4,5,6)
  nLayer = length(years)

  load('data/checkRaster')
  r <- subset( checkRaster, c(years,years+6))
  #r <- subset( checkRaster, c(1))
  x <- cbind(xyFromCell(r,1:ncell(r)),values(r))
  samplingRate <- 11
      
  r.row <- nrow(r)
  r.col <- ncol(r)
  
  # random starts
  k.row <- sample(1:r.row,size=1) 
  k.col <- sample(1:r.col,size=1) 
      
  row.values <- sort(unique(c( seq( from=k.row, to=1, by=-1*samplingRate), seq( from=k.row, to=r.row, by=samplingRate) )))
  col.values <- sort(unique(c( seq( from=k.col, to=1, by=-1*samplingRate), seq( from=k.col, to=r.col, by=samplingRate) )))
  
  # create sample
  sampled <- sapply( col.values, function(r) r + r.col * (row.values -1) )

  x.sampled <- x[sampled,]

  bandwidth=c(1000,1000,rep(0.04,nLayer), rep(2,nLayer))


#  sink('test.txt')
  nn <- getNN(
    x.sampled, 
    x.sampled,
    500,
    #nrow(x.sampled),
    kernelMethod = "HYBRID", 
    bandwidth=bandwidth,
    alpha=0.3,
    iterations = 100, 
    interpolate = x, 
    epsilon= 0.5,
    debugTrain = F 
  ) 
#  sink()
  unique.assignment <- unique(nn$assignment)
  reduced.assignment <- sapply( nn$assignment, function(x) ( which(x == unique.assignment) ) )
  colors.assignment <- sample(rainbow( length(unique.assignment)))[reduced.assignment]

  #plot initial picture
  plot( subset(r,1), col=two.colors(256, end='white',start='black', middle='grey')) 
  points( x.sampled[,1:2], col=colors.assignment, pch=16)
  
  png('interpolation.png')
  a <- subset(r,1)
  values(a) <- as.factor(nn$interpolateAssignment)
  plot(a, col=larry.colors() )
  dev.off()

  if( FALSE ) {


#

function( nn, x.sampled ) {

  n <- NROW(x.sampled)

  library(fields)
  if( !(2 %in% dev.list())) x11() ; dev.set(2)

  # handle plot colors
  unique.assignment <- unique(nn$assignment)
  reduced.assignment <- sapply( nn$assignment, function(x) ( which(x == unique.assignment) ) )
  colors.assignment <- sample(rainbow( length(unique.assignment)))[reduced.assignment]

  #plot initial picture
  #plot( r, col=two.colors(256, end='white',start='black', middle='grey')) 
  plot( subset(r,1), col=two.colors(256, end='white',start='black', middle='grey')) 
  points( x.sampled[,1:2], col=colors.assignment, pch=16)
  
#  plot( subset(r,1), col=two.colors(256, start='white',end='black', middle='grey')) 
#  points( nn$value[,1:2], col=as.factor(nn$value[,3]), pch=16)

  # use click to sleect a point
  closestPoint <- click(r,n=1,xy=TRUE)
  closestPoint <- which.min(sqrt(colSums(((t(x.sampled[,1:2]) - as.numeric(closestPoint[1:2])))^2)))

  # mark selected point
  points( x.sampled[closestPoint,1], x.sampled[closestPoint,2] , pch=4, cex=2, col='purple')

  # iterate throug a number or steps plotting important weights, with their magnituded between colors
  for( i in 1:50) {

    closestPointIndex <- n * (i-1) + closestPoint
  
    if( !((2+i) %in% dev.list())) x11() 
    dev.set(2+i)
     
    nNeighbors <- length(nn$weightsDebug[closestPointIndex,-1])
    weight <-  nn$weightsDebug[closestPointIndex,-1] /max(nn$weightsDebug[closestPointIndex,-1]) 
  
    neighbors <- nn$neighborsDebug[closestPointIndex,-1] 
    neighbors <- neighbors[weight > 10^-5]
    weight <- weight[weight > 10^-5]
    
    weightInt <- floor( nNeighbors * weight) +1
    weightColor = two.colors(n=max(weightInt), start='red', end='yellow', middle='orange')
    
    plot( subset(r,1) , col=two.colors(256, start='white',end='black', middle='grey')) 
    title( sprintf("iteration = %d",i))
    points( x.sampled[ neighbors + 1 ,1:2] , col= weightColor[weightInt] , pch=16) 
   
    # get point moved to 
    if( i == 1) {
      points( x.sampled[closestPointIndex,1], x.sampled[closestPointIndex,2] , pch=4, cex=2, col='purple')
    } else {
      points( value[1], value[2] , pch=4, cex=2, col='red')
    }
    value <- nn$valuesDebug[closestPointIndex,-1] 
    points( value[1], value[2] , pch=4, cex=2, col='green')
  }




  # change from 14 254 to 15 232
  iter1 <- 20 
  iter2 <- 8 
  n1 <- nn$assignmentsDebug[nn$assignmentsDebug[,1]==iter1,][closestPoint,]
  n2 <- nn$assignmentsDebug[nn$assignmentsDebug[,1]==iter2,][closestPoint,]
  print(n1)
  print(n2)
 
  # the goal here is to get all prior points
  m1 <- nn$assignmentsDebug[ (nn$assignmentsDebug[,1] == 1) ,2] == n1[2]
  m1 <- c(rep(m1,iter1),rep(FALSE, length(m1)*(50-iter1)))
  
  m2 <- nn$assignmentsDebug[ (nn$assignmentsDebug[,1] == iter1) ,2] == n2[2]
  m2 <- c(rep(m2,iter1+1),rep(FALSE, length(m2)*(100-iter1 -1)))


  v1 <- nn$valuesDebug[m1,,drop=FALSE]
  v2 <- nn$valuesDebug[m2,,drop=FALSE]

  array( v2[,-1], 1, function(x) sum(( t(v1[,-1]) - x)^2))


  nn$neighborsDebug[nn$neighborsDebug[,1]==iter1,][closestPoint,1:10]
  nn$neighborsDebug[nn$neighborsDebug[,1]==iter2,][closestPoint,1:10]
  
  nn$weightsDebug[nn$weightsDebug[,1]==iter1,][closestPoint,1:10]
  nn$weightsDebug[nn$weightsDebug[,1]==iter2,][closestPoint,1:10]


  a1 <- (nn$assignmentsDebug[,2] == n1[2]) & (nn$assignmentsDebug[,1] == iter1)
  a2 <- (nn$assignmentsDebug[,2] == n2[2]) & (nn$assignmentsDebug[,1] == iter2)
  a3 <- (nn$assignmentsDebug[,2] == n2[2]) & (nn$assignmentsDebug[,1] == iter1)

  print(nn$valuesDebug[a1,,drop=FALSE][1,])
  print(nn$valuesDebug[a2,,drop=FALSE][1,])
  
  print(sqrt(sum((nn$valuesDebug[a1,,drop=FALSE][1,-1]/bandwidth - nn$valuesDebug[a3,,drop=FALSE][1,-1]/bandwidth)^2)))
  
  nn$neighborsDebug[a1,1:10]
  nn$neighborsDebug[a2,1:10]
  nn$weightsDebug[a1,1:10]
  nn$weightsDebug[a2,1:10]
  
  
  
  a1 <- (nn$assignmentsDebug[,2] %in% c(223) ) & (nn$assignmentsDebug[,1] %in%  7)
  a2 <- (nn$assignmentsDebug[,2] == 232) & (nn$assignmentsDebug[,1] %in% 1:17)

  nn$valuesDebug[a1,]
  nn$valuesDebug[a2,]

  records <- 1:10 

  weight <- nn$weightsDebug[nn$weightsDebug[,1]==8,][closestPoint,1:10][-1]
  neighbors <- nn$neighborsDebug[nn$neighborsDebug[,1]==8,][closestPoint,1:10][-1]

  weightCheck <- apply( dnorm( t(x.sampled[ neighbors+1,])/bandwidth - nn$valuesDebug[a1,-1]/bandwidth) , 2, prod)

  print(weightCheck) 
#
#
#
#  ###### end check ####

} 

#  print( table(nn$neighborsQuery) )
#  print( length(table(nn$neighborsQuery)) )
#
#  r.new <- r
#  values(r.new) <- nn$interpolateAssignment
#  plot(r.new)
#
#  z.new <- r
#  values(z.new) <- nn$interpolateValue[,3]
#  plot(z.new)
#

#  points(out$valueTrain[,1:2], col=as.factor(out$assignmentTrain))
    }
  }
}

if( "testShiftNN2" %in% ls() ) {
  if( testShiftNN2 ) {


    # load raster
    library(raster)
    r <- raster("~/src/flannR/data/test.gri")

    Ext <- c( 
    755774,
    755960,
    2064572,
    2064731
    )

    r <- crop(r, extent(Ext) )  

    # sample
    #sampled <- sample(1:ncell(r), round(ncell(r)/121))
    nNeighbors <- 500
    maxIter <- 300 
    bandwidth <- c(50,50,1)
    spatial <- T
    alpha <- 0
    sampled = 'grid'
    samplingRate = 11 
    epsilon= 0.0001 

    out <- sampleClassify( 
      r, 
      maxIter=maxIter, 
      sampled, 
      nNeighbors=nNeighbors, 
      bandwidth = bandwidth, 
      kernelMethod = "HYBRID",
      epsilon= epsilon, 
      spatial=spatial, 
      alpha=alpha,
      samplingRate=samplingRate
    ) 


    print(length(unique(out$assignmentTrain)))
    print(length(unique(out$assignment)))
    print(length(unique(out$valueTrain[,3])))
    print(length(unique(out$value[,3])))

    plot(out$assignment)
    points(out$valueTrain[,1:2], col=out$assignmentTrain)

  }
}



if( "testShiftNN3" %in% ls() ) {
  if( testShiftNN3 ) {

#    for(i in 2:8) {
#    for(i in 2:4) {
#      if( !( i %in% dev.list()) ) x11()
#    }

    # load raster
    library(raster)
    library(fields)
    r <- raster("~/src/flannR/data/test.gri")/100

    Ext <- c( 
    755774,
    755960,
    2064572,
    2064731
    )
   
    Ext.poly <- cbind(    
      c(rep(Ext[1:2],2), Ext[1]) , 
      c( rep(Ext[3:4],each=2), Ext[3])
    )
    Ext.poly <- Ext.poly[c(1,3,4,2,5),]
    
    png(file=sprintf("~/src/dissertation/plots/NNS.png"))
    plot(r, col=two.colors(256, 'black','white','grey'))
    points( Ext.poly, col='red',type='l')
    dev.off()

    r <- crop(r, extent(Ext) )  

    # add on coordinates
    x <- cbind(xyFromCell(r,1:ncell(r)),values(r))

    # sample
    draw <- sample(1:nrow(x), round(nrow(x)/225))
    x <- x[draw,]

  #devList <- c(2:4)
  x.list <- list()

  for( j in 1:5) {
       
   if( j == 1)  alpha = 0 
   if( j == 2)  alpha = 0.25 
   if( j == 3)  alpha = 0.5 
   if( j == 4)  alpha = 0.75 
   if( j == 5)  alpha = 1 


    #dev.set( devList[j] )
  
    png(file=sprintf("~/src/dissertation/plots/NNS%d.png",100*alpha))
    plot(r, col=two.colors(256, 'black','white','grey'))
    points(x[,1:2],col='red',cex=.3)

    x.list[[j]] <- c(x)
    nn <- getNN(
      x, 
      x, 
      nrow(x),
      kernelMethod = "HYBRID", bandwidth=c(25,25,.1),
      alpha = alpha 
    )
    for( i in 1:100 ) {
   
      nn <- getNN(
        nn$query, 
        x, 
        nrow(x),
        kernelMethod = "HYBRID", bandwidth=c(25,25,.1), 
        alpha = alpha 
      )
      x.list[[j]] <- rbind(x.list[[j]],c(nn$query))
    }
    for(i in 1:nrow(x)) points( x.list[[j]][, c(i,nrow(x) + i)],col='purple',type='l' ) 
    #points( x.list[[j]][, c(1,nrow(x) + 1)],type='l' ,col='red', lwd='2') 
    points(nn$query[,1:2], col='blue',pch=16,cex=2)
    title(sprintf("alpha = %4.2f", alpha))

    dev.off()
 
  }  


}
}

#dev.set(6)
#k <- 1
#k.offset <- c(k, nrow(x) + k, 2*nrow(x) + k )
#a <- cbind( 
#      x.list[[1]][,k.offset],
#      x.list[[2]][,k.offset],
#      x.list[[3]][,k.offset],
#      x.list[[4]][,k.offset]
#      )
#
#plot(1:nrow(a), a[,1],type='l')
#points(1:nrow(a), a[,4],type='l', col='blue')
#points(1:nrow(a), a[,7],type='l', col='purple')
#points(1:nrow(a), a[,10],type='l', col='red')
#
#dev.set(7)
#plot(1:nrow(a), a[,2],type='l')
#points(1:nrow(a), a[,5],type='l', col='blue')
#points(1:nrow(a), a[,8],type='l', col='purple')
#points(1:nrow(a), a[,11],type='l', col='red')
#
#dev.set(8)
#plot(1:nrow(a), a[,3],type='l')
#points(1:nrow(a), a[,6],type='l', col='blue')
#points(1:nrow(a), a[,9],type='l', col='purple')
#points(1:nrow(a), a[,12],type='l', col='red')

#    library(raster)
#    r <- raster("~/src/dissertation/src/rSegmentation/test.gri")
#
#    Ext <- c( 
#    755774,
#    755960,
#    2064572,
#    2064731
#    )
#
#    #r <- crop(r, extent(Ext) )  
#
#    source('~/src/dissertation/src/cSmooth/smooth.R')
#    h <- 17 
#    WMean   <- matrix(1/h^2,nrow=h,ncol=h)
#    WVar    <- matrix(1/h^2,nrow=h,ncol=h)
#    r <- rasterSmoothMoments(r, WMean, WVar)
#
#
#    r.values <- cbind(
#      #rep(1:ncol(r$mu),nrow(r$mu)),
#      #rep(1:nrow(r$mu),each=ncol(r$mu)),
#      values(r$mu),
#      values(log(r$var))
#      )
#
#     r.values.finite <- is.finite(rowSums(r.values)) 
# 
#     x <- r.values[ r.values.finite, ]
#    
#     sampleRate <- 40 
#     #bandwidth <- c(.1,.1,20,2)
#     bandwidth <- c(21,1)
# 
#     indexCol <- seq.int(from=1,to=ncol(r$mu),by=sampleRate)
#     indexRow <- seq.int(from=1,to=nrow(r$mu),by=sampleRate)
#
#     index <- c( indexCol %o% indexRow )
#     
# 
#     #index <- sample( 1:nrow(x), sampleSize ) 
#     train <- x[index,]
#    
#     nn <- getNN(train, train, 
#                 #length(index), 
#                 100,
#                 kernelMethod = "NORMAL" , 
#                 bandwidth=bandwidth)
# 
#     for( i in 1:50) {
#       nn <- getNN(nn$query, train, 
#                   #length(index), 
#                   100,
#                   kernelMethod="NORMAL", 
#                   bandwidth=bandwidth)
#     }
#    
#     # classify   
#     y <- getNN(x, nn$query, 1)
#     
# 
# #    # copy the results back
#     r.new <-r$mu
#     r.new.values <- values(r.new) 
#     r.new.values[r.values.finite] <- y$query[,1]
#     values(r.new) <- r.new.values
#     
#     r.new.var <-r$var
#     r.new.values <- values(r.new) 
#     r.new.values[r.values.finite] <- y$query[,2]
#     values(r.new.var) <- r.new.values
# 
#     par( mfrow=c(2,2) )
#     plot(r$mu, zlim = c(0,100))
#     plot(r.new, zlim=c(0,100) )
#     plot(log(r$var), zlim = c(-2,8) )
#     plot(r.new.var, zlim = c(-2, 8) )
#
#    nCheck <- 2:4
#    
#    y <- x[nCheck,]
#    
#    check <- 3 
#    
#    y.check <-  t(y) - x[check,]
#    y.check.w <- dnorm(y.check)
#    y.check.w <- y.check.w[1,] * y.check.w[2,]
#    
#    output <- colSums(y * y.check.w) / sum(y.check.w)
#    print(output)
  
#  }
#}

####################################################################################################
if( "testShiftNN4" %in% ls() ) {
  if( testShiftNN4 ) {

    # example data
    x <- matrix( c( 
               1, 1,
               1, 2,
               2, 1, 
               5, 3,
               5, 4,
               4, 5 ), byrow=T, ncol=2)


  #x1 <- matrix( rnorm( 10 ),ncol=2)
  #x2 <- matrix( rnorm( 10 ),ncol=2) + 2 
  #x <- rbind( x1, x2 ) 

  nn <- getNN(
    x, 
    x, 
    #nrow(x),
    nNeighbors=3,
    kernelMethod = "HYBRID", 
    bandwidth=c(3,3),
    alpha=0,
    iterations = 1 ,
    returnDistance= T
  ) 
  print(nn)


  for( i in 1:2 ) {
  nn <- getNN(
    nn$query, 
    x, 
    #nrow(x),
    nNeighbors=3,
    kernelMethod = "HYBRID", 
    bandwidth=c(3,3),
    alpha=0,
    iterations = 1 ,
    returnDistance= T
  ) 
  }
  print(nn)


  # neighbor restriction check (only 3 neighbors)
  w <- apply( cbind( dnorm( (x[1,1] - x[1:3,1])/3 ), dnorm( (x[1,2] - x[1:3,2])/3 )), 1, prod)
  x1 <- sum(x[1:3,1] * w)/sum(w)
  x2 <- sum(x[1:3,2] * w)/sum(w)
  print(x1)
  print(x2)
  
  # no neighbor restriction check 
  w <- apply( cbind( dnorm( (x[1,1] - x[,1])/3 ), dnorm( (x[1,2] - x[,2])/3 )), 1, prod)
  x1 <- sum(x[,1] * w)/sum(w)
  x2 <- sum(x[,2] * w)/sum(w)
  print(x1)
  print(x2)
 

  }
}


