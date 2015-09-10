\name{getNN}
\alias{getNN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Mean Shift
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
meanShift(queryData, trainData, nNeighbors = nrow(trainData), indexParameters = NULL, searchParameters = NULL, kernelMethod = "NONE", bandwidth = NULL, alpha = 0.1, iterations = 1, epsilon = 1e-04, interpolate, debugTrain = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{queryData}{
%%     ~~Describe \code{queryData} here~~
}
  \item{trainData}{
%%     ~~Describe \code{trainData} here~~
}
  \item{nNeighbors}{
%%     ~~Describe \code{nNeighbors} here~~
}
  \item{indexParameters}{
%%     ~~Describe \code{indexParameters} here~~
}
  \item{searchParameters}{
%%     ~~Describe \code{searchParameters} here~~
}
  \item{kernelMethod}{
%%     ~~Describe \code{kernelMethod} here~~
}
  \item{bandwidth}{
%%     ~~Describe \code{bandwidth} here~~
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
}
  \item{iterations}{
%%     ~~Describe \code{iterations} here~~
}
  \item{epsilon}{
%%     ~~Describe \code{epsilon} here~~
}
  \item{interpolate}{
%%     ~~Describe \code{interpolate} here~~
}
  \item{debugTrain}{
%%     ~~Describe \code{debugTrain} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (queryData, trainData, nNeighbors = nrow(trainData), 
    indexParameters = NULL, searchParameters = NULL, kernelMethod = "NONE", 
    bandwidth = NULL, alpha = 0.1, iterations = 1, epsilon = 1e-04, 
    interpolate, debugTrain = F) 
{
    trainRow <- nrow(trainData)
    queryRow <- nrow(queryData)
    queryCol <- ncol(queryData)
    if (is.null(queryRow)) {
        if (length(queryData) != queryCol) {
            stop("Error: queryData columns do not equal trainData columns")
        }
        queryData <- matrix(queryData, nrow = 1)
    }
    nNeighbors = min(nNeighbors, nrow(trainData))
    neighbors = rep(0, queryRow * nNeighbors)
    distances = rep(0, queryRow * nNeighbors)
    if (iterations < 1) {
        stop("Error: iterations are not at least 1")
    }
    prob <- rep(-1, queryRow)
    neighborsQuery = rep(-1, queryRow)
    distancesQuery = rep(-1, queryRow)
    integerSearchParameters <- 0
    doubleSearchParameters <- 0
    integerAlgorithmParameters <- 4
    doubleAlgorithmParameters <- 0
    paramKernelMethod <- .returnKernelEnum(kernelMethod)
    if (kernelMethod \%in\% c("NEWTON", "HYBRID")) {
        if (is.null(bandwidth)) 
            stop("bandwidth NULL for kernel method")
        if (length(bandwidth) != queryCol) 
            stop("bandwidth not equal to columns x")
        if (length(bandwidth) > 1) {
            queryData <- queryData \%*\% diag(1/bandwidth)
            trainData <- trainData \%*\% diag(1/bandwidth)
        }
        else {
            queryData <- queryData/bandwidth
            trainData <- trainData/bandwidth
        }
    }
    if (!missing(interpolate)) {
        print("interpolateFound")
        if (kernelMethod \%in\% c("NEWTON", "HYBRID")) {
            if (is.null(bandwidth)) 
                stop("bandwidth NULL for kernel method")
            if (length(bandwidth) != queryCol) 
                stop("bandwidth not equal to columns x")
            if (length(bandwidth) > 1) {
                interpolate <- interpolate \%*\% diag(1/bandwidth)
                interpolateIndex <- 1:nrow(interpolate)
            }
            else {
                interpolate <- interpolate/bandwidth
                interpolateIndex <- 1:nrow(interpolate)
            }
        }
    }
    else {
        interpolate <- -1
        interpolateIndex <- -1
    }
    if (debugTrain) {
        print("Warning Debugging is being performed, this may exhaust system memory for large data sets")
        assignmentsDebug <- as.integer(rep(1, NROW(trainData) * 
            iterations))
        weightsDebug <- as.double(rep(0, NROW(trainData) * iterations * 
            nNeighbors))
        neighborsDebug <- as.integer(rep(0, NROW(trainData) * 
            iterations * nNeighbors))
        valuesDebug <- as.double(rep(0, NROW(trainData) * iterations * 
            ncol(trainData)))
    }
    else {
        assignmentsDebug <- as.integer(rep(0, 1))
        weightsDebug <- as.double(rep(0, 1))
        neighborsDebug <- as.integer(rep(0, 1))
        valuesDebug <- as.double(rep(0, 1))
    }
    Cprog <- proc.time()
    r.result <- .C("R_meanShiftNN", as.double(t(queryData)), 
        as.double(t(trainData)), as.integer(neighbors), as.double(distances), 
        as.double(prob), as.integer(neighborsQuery), as.double(distancesQuery), 
        as.integer(queryRow), as.integer(trainRow), as.integer(queryCol), 
        as.integer(nNeighbors), as.integer(integerSearchParameters), 
        as.double(doubleSearchParameters), as.integer(integerAlgorithmParameters), 
        as.double(doubleAlgorithmParameters), as.integer(paramKernelMethod), 
        as.double(bandwidth), as.double(alpha), as.integer(iterations), 
        as.double(epsilon), as.double(t(interpolate)), as.integer(interpolateIndex), 
        as.integer(nrow(interpolate)), assignmentsDebug, weightsDebug, 
        neighborsDebug, valuesDebug)
    print("C running time")
    print(proc.time() - Cprog)
    prob <- matrix(r.result[[7]], ncol = 1)
    assignment <- matrix(r.result[[6]], byrow = TRUE, ncol = 1)
    value <- matrix(r.result[[1]], byrow = TRUE, ncol = queryCol)
    if (interpolateIndex[1] >= 0) {
        interpolate <- matrix(r.result[[21]], byrow = TRUE, ncol = queryCol)
        interpolateIndex <- matrix(r.result[[22]], byrow = TRUE, 
            ncol = 1)
        if (!debugTrain) {
            return(list(assignment = assignment, value = value, 
                probability = prob, interpolateAssignment = interpolateIndex, 
                interpolateValue = interpolate))
        }
    }
    if (debugTrain) {
        assignmentsDebug <- cbind(rep(1:iterations, each = trainRow), 
            matrix(r.result[[24]], byrow = TRUE, ncol = 1))
        weightsDebug <- cbind(rep(1:iterations, each = trainRow), 
            matrix(r.result[[25]], byrow = TRUE, ncol = nNeighbors))
        neighborsDebug <- cbind(rep(1:iterations, each = trainRow), 
            matrix(r.result[[26]], byrow = TRUE, ncol = nNeighbors))
        valuesDebug <- cbind(rep(1:iterations, each = trainRow), 
            matrix(r.result[[27]], byrow = TRUE, ncol = queryCol))
        if (interpolateIndex[1] >= 0) {
            return(list(assignment = assignment, value = value, 
                probability = prob, interpolateAssignment = interpolateIndex, 
                interpolateValue = interpolate, assignmentsDebug = assignmentsDebug, 
                weightsDebug = weightsDebug, neighborsDebug = neighborsDebug, 
                valuesDebug = valuesDebug))
        }
        else {
            return(list(assignment = assignment, value = value, 
                probability = prob, assignmentsDebug = assignmentsDebug, 
                weightsDebug = weightsDebug, neighborsDebug = neighborsDebug, 
                valuesDebug = valuesDebug))
        }
    }
    return(list(assignment = assignment, value = value, probability = prob))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
