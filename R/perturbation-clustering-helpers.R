GetClusteringAlgorithm <- function(clusteringMethod = "kmeans", clusteringFunction = NULL, clusteringOptions = NULL, ...){
    name = clusteringMethod
    
    if (!is.function(clusteringFunction)) {
        switch(
            clusteringMethod,
            kmeans = {
                clusteringFunction <- kmeansWrapper
            },
            pam = {
                clusteringFunction <- pamWrapper
            },
            hclust = {
                clusteringFunction <- hclustWrapper
            },
            {
                stop("clusteringMethod not found. Please pass a clusteringFunction instead")
            }
        )
    }
    else {
        name = "Unknow"
    }
    
    list(
        fun = function(data, k) do.call(clusteringFunction, c(list(data = data, k = k), clusteringOptions)),
        name = name
    )
}

GetPerturbationAlgorithm <- function(data = data, perturbMethod = "noise", perturbFunction = NULL, perturbOptions = NULL, ...){
    name = perturbMethod
    
    if (!is.function(perturbFunction)) {
        switch(
            perturbMethod,
            noise = {
                if (is.null(perturbOptions))
                    perturbOptions <- list()
                
                noise = perturbOptions$noise
                if (is.null(noise))
                    noise <- GetNoise(data, noisePercent = perturbOptions$noisePercent)
                
                # get perturbed similarity
                perturbFunction <- function(data, ...) {
                    AddNoisePerturb(data = data, noise = noise)
                }
            },
            subsampling = {
                perturbFunction <- SubSampling
            },
            {
                stop("perturbMethod not found. Please pass a perturbFunction instead")
            }
        )
    }
    else {
        name = "Unknow"
    }
    
    list(
        fun = function(data) do.call(perturbFunction, c(list(data = data), perturbOptions)),
        name = name
    )
}

BuildConnectivityMatrix <- function(data, clus, clusteringAlgorithm) {
    rowNum <- nrow(data)
    S <- matrix(0, rowNum, rowNum)
    cluster <- clusteringAlgorithm(data, clus)
    
    for (j in 1:clus) {
        X <- rep(0, rowNum)
        X[which(cluster == j)] <- 1
        S <- S + X %*% t(X)
    }
    
    rownames(S) <- rownames(data)
    colnames(S) <- rownames(data)
    
    list(matrix = S, groups = cluster)
}

GetOriginalSimilarity <- function(data, clusRange, clusteringAlgorithm, showProgress = F, ncore) {
    data <- prcomp(data)$x
    groupings <- list()
    origS <- list()
    
    .ncore = if (.Platform$OS.type == "unix") ncore else 1
    
    seeds = abs(round(rnorm(max(clusRange))*10^6))
    params = list(
        X = clusRange,  
        mc.cores = if (length(clusRange) > .ncore) .ncore else length(clusRange),  
        FUN = function(clus) {
            set.seed(seeds[clus])
            list(cMatrix = BuildConnectivityMatrix(data, clus, clusteringAlgorithm),clus = clus)
        }
    )
    
    if(showProgress) {
        params$mc.style = if(showProgress) "txt" else NULL
    }
    
    rets <- do.call((if(showProgress) pbmclapply else mclapply), params)
    
    for(ret in rets){
        origS[[ret$clus]] <- ret$cMatrix$matrix
        groupings[[ret$clus]] <- ret$cMatrix$group
    }
    
    list(origS = origS, groupings = groupings)
}

GetPerturbedSimilarity <- function(data, clusRange, iterMax, iterMin, clusteringAlgorithm, perturbedFunction, stoppingCriteriaHandler, showProgress = F, ncore) {
    pertS <- list()
    currentIter <- rep(0,max(clusRange)) # this iter will be used for merge connectivity matrix
    
    for (clus in clusRange) {
        pertS[[clus]] <- matrix(0, nrow(data), nrow(data))
        rownames(pertS[[clus]]) <- rownames(data)
        colnames(pertS[[clus]]) <- rownames(data)
    }
    
    jobs <- rep(clusRange, iterMax)
    maxJob = length(jobs)
    
    if (showProgress) pb <- txtProgressBar(min = 0, max = maxJob, style = 3)
    
    kProgress = rep(0, max(clusRange))
    perturbedRets <- list()
    
    parCluster <- if (.Platform$OS.type == "unix") makeForkCluster(ncore) else makePSOCKcluster(ncore)
    
    if (.Platform$OS.type != "unix") {
        clusterExport(parCluster, varlist = c("BuildConnectivityMatrix"), envir = environment())
    }
    
    registerDoParallel(parCluster)
    
    step = iterMin
    while (length(jobs) > 0) {
        jobLength = step*length(clusRange)
        step = 10

        if (jobLength < ncore * step){
            jobLength = ncore * step
        }
        
        if (jobLength > length(jobs)){
            jobLength =  length(jobs)
        }

        currentJobs = list()
        count = 1
        for (clus in jobs[1:jobLength]){
            kProgress[clus] = kProgress[clus] + 1
            currentJobs[[count]] <- list(iter = kProgress[clus], clus = clus, seed = abs(round(rnorm(1)*10^6))) # this iter will be use for get perturbed data from perturbedRets
            count = count + 1
        }

        if (jobLength < length(jobs)){
            jobs <- jobs[(jobLength+1):length(jobs)]
        }
        else jobs <- c()
        
        iters = as.numeric(unique(lapply(currentJobs, function(j) j$iter)))
        minIter = min(iters, kProgress[which(kProgress != -1)])
        maxIter = max(iters, kProgress)
        
        for (i in 1:maxIter){
            if (i < minIter){
                perturbedRets[[i]] <- 0
            }
            else if (is.null(perturbedRets[i][[1]])){
                perturbedRets[[i]] <- perturbedFunction(data = data)
            }
        }
        job <- NULL
        rets <- foreach(job = currentJobs) %dopar% {
            clus = job$clus
            set.seed(job$seed)
            perturbedRet <- perturbedRets[[job$iter]]
            cMatrix <- BuildConnectivityMatrix(data = perturbedRet$data, clus, clusteringAlgorithm)
            connectivityMatrix = perturbedRet$ConnectivityMatrixHandler(connectivityMatrix = cMatrix$matrix, iter = job$iter, k = clus)

            list(
                connectivityMatrix = connectivityMatrix,
                clus = clus,
                perturbedRet = perturbedRet
            )
        }
        
        for(ret in rets){
            clus = ret$clus
            
            currentIter[clus] <- currentIter[clus] + 1
            pertS[[clus]] <- ret$perturbedRet$MergeConnectivityMatrices(oldMatrix = pertS[[clus]], newMatrix = ret$connectivityMatrix, iter = currentIter[clus], k = clus)
            
            stop <- stoppingCriteriaHandler(iter = currentIter[clus], k = clus, pert = pertS[[clus]])
            if (stop) {
                if (showProgress) setTxtProgressBar(pb, getTxtProgressBar(pb) + length(jobs[jobs == clus]))
                jobs <- jobs[jobs != clus]
                kProgress[clus] <- -1
            }
            
            if (showProgress) setTxtProgressBar(pb, getTxtProgressBar(pb) + 1)
        }
    }
    
    if (showProgress) {
        setTxtProgressBar(pb, maxJob)
        cat("\n")
    }
    
    stopCluster(parCluster)
    
    pertS
}

CalcAUC <- function(orig, pert) {
    N <- nrow(orig)
    S <- abs(orig - pert)
    diag(S) <- 0
    # added -10^(-5) for visual purposes
    A <- c(-10^(-5), sort(unique(as.numeric(S))))
    if (max(A) < 1)
        A <- c(A, 1)
    B <- NULL
    for (i in 1:length(A)) {
        B[i] <- sum(S <= A[i])/(N * N)
    }
    
    area <- 0
    for (i in (2:length(A))) {
        area <- area + B[i - 1] * (A[i] - A[i - 1])
    }
    
    list(area = area, entry = A, cdf = B)
}

CalcPerturbedDiscrepancy <- function(origS, pertS, clusRange) {
    diff <- NULL
    for (clus in clusRange) {
        diff[clus] <- sum(abs(origS[[clus]] - pertS[[clus]]))
    }
    
    AUC <- NULL
    entries <- list()
    cdfs <- list()
    
    for (clus in clusRange) {
        ret <- CalcAUC(origS[[clus]], pertS[[clus]]);
        
        entries[[clus]] <- ret$entry
        cdfs[[clus]] <- ret$cdf
        AUC[clus] <- ret$area
    }
    
    list(Diff = round(diff, digits = 10), Entry = entries, CDF = cdfs, AUC = round(AUC, digits = 10))
}