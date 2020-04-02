FindMaxHeight <- function (hc, maxK) {
    height <- rev(diff(hc$height))[1:(maxK-1)]
    i=which.max(height)
    i+1
}

GetSimilarityFromGrouping <- function(g) {
    N = length(g)
    S = matrix(0, N, N)
    colnames(S) = rownames(S) = names(g)
    for (j in unique(g)) {
        X = rep(0, N);
        X[which(g == j)] = 1
        S = S + X %*% t(X)
    }
    S
}

CalcClusterAgreement <- function(l) {
    N <- length(l[[1]])
    A <- matrix(0, N, N)
    for (group in l) {
        A <- A + GetSimilarityFromGrouping(group)
    }
    A = A / (length(l))
    ret = (sum(A == 0) + sum(A == 1)) / (N ^ 2)
    
    ret
}

ClusteringFunWrapper <- function(kMax, groupings, FUN){
    clusters = list()
    
    agreeS <- lapply(2 : kMax, FUN = function(i){
        cluster = FUN(i)
        clusters[[i]] <<- cluster
        
        CalcClusterAgreement(
            c(
                groupings,
                list(cluster)
            )
        )
    })
    
    ret <- list()
    ret$agreeS <- unlist(agreeS)
    ret$k <- which.max(agreeS) + 1
    ret$agree <- ret$agreeS[ret$k - 1]
    ret$cluster <- clusters[[ret$k]]
    ret
}

ClusterUsingPAM <- function(orig, kMax, groupings) {
    ClusteringFunWrapper(kMax, groupings, FUN = function(k) pamWrapper(1 - orig, k, diss = T))
}

ClusterUsingHierarchical <- function(orig, kMax, groupings) {
    hcO <- hclust(as.dist(1 - orig), method = "average")
    ClusteringFunWrapper(kMax, groupings, FUN = function(k) cutree(hcO, k))
}
# 
# PClusGap <- function (x, FUNcluster, K.max, B = 100, d.power = 1, spaceH0 = c("scaledPCA", "original"), ncore) {
#     
#     stopifnot(is.function(FUNcluster), length(dim(x)) == 2, 
#               K.max >= 2, (n <- nrow(x)) >= 1, ncol(x) >= 1)
#     if (B != (B. <- as.integer(B)) || (B <- B.) <= 0) 
#         stop("'B' has to be a positive integer")
#     cl. <- match.call()
#     if (is.data.frame(x)) 
#         x <- as.matrix(x)
#     ii <- seq_len(n)
#     W.k <- function(X, kk) {
#         clus <- if (kk > 1) FUNcluster(X, kk)$cluster else rep.int(1L, nrow(X))
#         0.5 * sum(vapply(split(ii, clus), function(I) {
#             xs <- X[I, , drop = FALSE]
#             sum(dist(xs)^d.power/nrow(xs))
#         }, 0))
#     }
#     logW <- E.logW <- SE.sim <- numeric(K.max)
#     
#     parCluster <- if (.Platform$OS.type == "unix") makeForkCluster(ncore) else makePSOCKcluster(ncore)
#     registerDoParallel(parCluster)
#     
#     seeds = abs(round(rnorm(K.max)*10^6))
#     k <- NULL
#     ret <- foreach (k = 1:K.max) %dopar% {
#         set.seed(seeds[k])
#         log(W.k(x, k))
#     }
#     
#     for (k in 1:K.max) logW[k] <- ret[[k]]
#     
#     spaceH0 <- match.arg(spaceH0)
#     xs <- scale(x, center = TRUE, scale = FALSE)
#     m.x <- rep(attr(xs, "scaled:center"), each = n)
#     switch(spaceH0, scaledPCA = {
#         V.sx <- svd(xs, nu = 0)$v
#         xs <- xs %*% V.sx
#     }, original = {
#     }, stop("invalid 'spaceH0':", spaceH0))
#     rng.x1 <- apply(xs, 2L, range)
#     logWks <- matrix(0, B, K.max)
#     
#     seeds = abs(round(rnorm(B)*10^6))
#     b <- NULL
#     ret <- foreach (b = 1:B) %dopar% {
#         set.seed(seeds[b])
#         z1 <- apply(rng.x1, 2, function(M, nn) runif(nn, min = M[1], max = M[2]), nn = n)
#         z <- switch(spaceH0, scaledPCA = tcrossprod(z1, V.sx), original = z1) + m.x
#         
#         lapply(1:K.max, function(k){
#             log(W.k(z, k))
#         })
#     }
#     
#     stopCluster(parCluster)
#     
#     for (b in 1:B){
#         for (k in 1:K.max) {
#             logWks[b, k] <- ret[[b]][[k]]
#         }
#     }
#     
#     E.logW <- colMeans(logWks)
#     SE.sim <- sqrt((1 + 1/B) * apply(logWks, 2, var))
#     structure(class = "clusGap", list(Tab = cbind(logW, E.logW, gap = E.logW - logW, SE.sim), call = cl., spaceH0 = spaceH0, n = n, B = B, FUNcluster = FUNcluster))
# }
