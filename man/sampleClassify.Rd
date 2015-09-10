\name{sampleClassify}
\alias{sampleClassify}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sample Classify
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
sampleClassify(x, maxIter = 1, sampled = FALSE, nNeighbors = 500, bandwidth = NULL, kernelMethod = "HYBRID", epsilon = 0.01, spatial = F, alpha = 0, samplingRate = 11)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{maxIter}{
%%     ~~Describe \code{maxIter} here~~
}
  \item{sampled}{
%%     ~~Describe \code{sampled} here~~
}
  \item{nNeighbors}{
%%     ~~Describe \code{nNeighbors} here~~
}
  \item{bandwidth}{
%%     ~~Describe \code{bandwidth} here~~
}
  \item{kernelMethod}{
%%     ~~Describe \code{kernelMethod} here~~
}
  \item{epsilon}{
%%     ~~Describe \code{epsilon} here~~
}
  \item{spatial}{
%%     ~~Describe \code{spatial} here~~
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
}
  \item{samplingRate}{
%%     ~~Describe \code{samplingRate} here~~
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
function (x, maxIter = 1, sampled = FALSE, nNeighbors = 500, 
    bandwidth = NULL, kernelMethod = "HYBRID", epsilon = 0.01, 
    spatial = F, alpha = 0, samplingRate = 11) 
{
    if (is.character(sampled)) {
        if (toupper(sampled) == "GRID") {
            if (samplingRate > max(nrow(x), ncol(x))) 
                stop("can't sample rate exceeds observations")
            r.row <- nrow(x)
            r.col <- ncol(x)
            k.row <- sample(1:r.row, size = 1)
            k.col <- sample(1:r.col, size = 1)
            row.values <- sort(unique(c(seq(from = k.row, to = 1, 
                by = -1 * samplingRate), seq(from = k.row, to = r.row, 
                by = samplingRate))))
            col.values <- sort(unique(c(seq(from = k.col, to = 1, 
                by = -1 * samplingRate), seq(from = k.col, to = r.col, 
                by = samplingRate))))
            sampled <- sapply(col.values, function(x) x + r.col * 
                (row.values - 1))
        }
    }
    x.values <- values(x)
    if (spatial == T) 
        x.values <- cbind(xyFromCell(x, 1:NROW(x.values)), x.values)
    x.nona <- !is.na(rowSums(x.values))
    x.values.nona <- x.values[x.nona, ]
    if (spatial == F) {
        train <- matrix(x.values, ncol = nlayers(x))[sampled, 
            , drop = FALSE]
    }
    else {
        train <- matrix(x.values, ncol = nlayers(x) + 2)[sampled, 
            , drop = FALSE]
    }
    train <- train[rowSums(is.na(train)) == 0, , drop = FALSE]
    if (length(bandwidth) != ncol(train)) 
        stop("layers of image and bandwidth do not match")
    nNeighbors = min(nNeighbors, nrow(train))
    nn <- getNN(train, train, nNeighbors, kernelMethod = "HYBRID", 
        bandwidth = bandwidth, alpha = alpha, iterations = maxIter, 
        epsilon = 0.5, interpolate = x.values.nona, debugTrain = F)
    value <- x.values
    value[x.nona, ] <- nn$interpolateValue
    x <- subset(x, 1)
    x.values <- values(x)
    x.values[!x.nona] <- NA
    x.values[x.nona] <- nn$interpolateAssignment
    values(x) <- x.values
    return(list(assignment = x, value = value, train = train, 
        sampleTrain = sampled, assignmentTrain = nn$assignment, 
        valueTrain = nn$value, probabilityTrain = nn$probability))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
