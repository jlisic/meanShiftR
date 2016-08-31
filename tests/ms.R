
ms <- function (X, h, subset, thr = 1e-04, scaled = TRUE, iter = 200, 
    plotms = 2, or.labels = NULL, ...) 
{
    n <- dim(X)[1]
    d <- dim(X)[2]
    if (missing(subset)) {
        subset <- 1:n
    }
    s1 <- apply(X, 2, function(dat) {
        diff(range(dat))
    })
    if (missing(h)) {
        if (!scaled) {
            h <- s1/10
        }
        else {
            h <- 0.1
        }
    }
    if (length(h) == 1) {
        h <- rep(h, d)
    }
    if (scaled) {
        X <- sweep(X, 2, s1, "/")
    }
    if (d == 2 && plotms > 0) {
        if (missing(or.labels)) {
            plot(X, col = "grey70", ...)
        }
        else {
            plot(X, col = or.labels, ...)
        }
    }
    finals <- matrix(0, n, d)
    ncluster <- 0
    savecluster <- matrix(0, 0, d)
    cluster.label <- closest.label <- rep(0, n)
    if (length(h) == 1) {
        h <- rep(h, d)
    }
    for (i in subset) {
        # run mean shift for the first item 
        temp.ms <- ms.rep(X, X[i, ], h, plotms = 0, thresh = 1e-08, iter)

        # record the final lcation to the first row of finals 
        finals[i, ] <- temp.ms$final

        cluster.dist <- rep(0, ncluster)
        # check to see if there is a prior cluster that is close
        # this creates a slight bias induced by the order of the points
        if (ncluster >= 1) {
            for (j in 1:ncluster) {
                cluster.dist[j] <- enorm(savecluster[j, ] - finals[i, ])/enorm(savecluster[j, ])
            }
        }
        if (ncluster == 0 || min(cluster.dist) > thr) {
            ncluster <- ncluster + 1
            savecluster <- rbind(savecluster, finals[i, ])
            cluster.label[i] <- ncluster
        }
        else {
            cluster.label[i] <- which(cluster.dist == min(cluster.dist))
        }
        if (d == 2 && plotms == 1) {
            lines(rbind(temp.ms$start, temp.ms$Meanshift.points), 
                col = "grey30")
        }
        if (d == 2 && plotms == 2) {
            lines(rbind(temp.ms$start, temp.ms$Meanshift.points), 
                col = cluster.label[i] + 1)
        }
    }
    for (i in subset) {
        closest.label[i] <- mindist(savecluster, X[i, ])$closest.item
    }
    if (d == 2 && plotms == 1) {
        points(finals, pch = 15, col = 2)
    }
    if (d == 2 && plotms == 2) {
        points(finals, pch = 15)
    }
    if (d > 2 && plotms > 1) {
        pairs(rbind(as.matrix(X), savecluster), col = c(cluster.label + 
            1, rep(1, dim(savecluster)[1])), pch = c(rep(20, 
            dim(X)[1]), rep(24, dim(savecluster)[1])), ...)
    }
    dimnames(savecluster) <- list(1:ncluster, NULL)
    fit <- list(cluster.center = savecluster, cluster.label = cluster.label, 
        closest.label = closest.label, h = h, data = X, scaled = scaled, 
        scaled.by = if (scaled) s1 else rep(1, d))
    class(fit) <- "ms"
    return(fit)
}
