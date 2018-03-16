GetNoise <- function(data, noisePercent = "median") {
    if (is.null(noisePercent)) noisePercent = "median"
    if (noisePercent == "median") {
        sds <- apply(data, 2, sd)
        noise <- median(sds)
    } else {
        sds <- apply(data, 2, sd)
        sds <- sort(sds)
        ind <- round(length(sds) * noisePercent/100)
        noise <- sds[ind]
    }
    noise
}

AddNoisePerturb <- function(data, noise) {
    rowNum <- nrow(data)
    colNum <- ncol(data)
    epsilon <-
        matrix(
            data = rnorm(rowNum * colNum, mean = 0, sd = noise),
            nrow = rowNum,
            ncol = colNum
        )
    
    list(
        data = data + epsilon,
        ConnectivityMatrixHandler = function(connectivityMatrix, ...) {
            connectivityMatrix
        },
        MergeConnectivityMatrices = function(oldMatrix, newMatrix, iter, ...){
            return((oldMatrix*(iter-1) + newMatrix)/iter)
        }
    )
}

SubSampling <- function(data, percent = 80) {
    N = nrow(data)
    randCount = ceiling(N * percent/100)
    randOrder = sample(N, N)
    
    randMatrix = data[randOrder[1:randCount], ]
    list(
        data = randMatrix,
        ConnectivityMatrixHandler = function(connectivityMatrix, ...) {
            S = matrix(0, N, N)
            S[1:randCount, 1:randCount] <- connectivityMatrix
            for(i in 1:N){
                S[N,N] = 1
            }
            S[randOrder, randOrder]
        },
        MergeConnectivityMatrices = function(oldMatrix, newMatrix, iter, ...){
            rand = randOrder[1:randCount]
            oldMatrix[rand,rand] <- oldMatrix[rand,rand]*(iter-1)/iter
            oldMatrix + newMatrix/iter
        }
    )
}