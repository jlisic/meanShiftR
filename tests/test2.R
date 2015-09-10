  # load raster
  library(raster)
 
  years <- c(1,2,3,4,5,6)
  nLayer = length(years)

  data(checkRaster)

  r <- subset( checkRaster, c(years,years+6))
  x <- cbind(xyFromCell(r,1:ncell(r)),values(r))
  samplingRate <- 11
      
  r.row <- nrow(r)
  r.col <- ncol(r)
  
  # random starts
  # create a systematic sample
  k.row <- sample(1:r.row,size=1) 
  k.col <- sample(1:r.col,size=1) 
  
  row.values <- sort(unique(c( seq( from=k.row, to=1, by=-1*samplingRate), seq( from=k.row, to=r.row, by=samplingRate) )))
  col.values <- sort(unique(c( seq( from=k.col, to=1, by=-1*samplingRate), seq( from=k.col, to=r.col, by=samplingRate) )))
  
  # take sample
  sampled <- sapply( col.values, function(r) r + r.col * (row.values -1) )

  x.sampled <- x[sampled,]
  
  # set bandwidth 
  bandwidth=c(1000,1000,rep(0.04,nLayer), rep(2,nLayer))


  nn <- getNN(
    x.sampled, 
    x.sampled,
    500,
    kernelMethod = "HYBRID", 
    bandwidth=bandwidth,
    alpha=0.3,
    iterations = 100, 
    interpolate = x, 
    epsilon= 0.5,
    debugTrain = F 
  ) 
  
  unique.assignment <- unique(nn$assignment)
  reduced.assignment <- sapply( nn$assignment, function(x) ( which(x == unique.assignment) ) )
  colors.assignment <- sample(rainbow( length(unique.assignment)))[reduced.assignment]

  #plot initial picture
  plot( subset(r,1), col=two.colors(256, end='white',start='black', middle='grey')) 
  points( x.sampled[,1:2], col=colors.assignment, pch=16)

  # plot assignment  
  a <- subset(r,1)
  values(a) <- as.factor(nn$interpolateAssignment)
  plot(a, col=larry.colors() )

