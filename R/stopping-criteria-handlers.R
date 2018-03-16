DiffDevianceStoppingHandler <- function(kMax, origS, iterMin = 20, madMin = 1e-03, msdMin = 1e-06){
    listAUC <- list()
    isStopped <- rep(F, kMax)
    
    for (k in 2:kMax) {
        listAUC[k] <- list(c())
    }

    aucMeanSquareDeviation <- function(arr) sum((arr - mean(arr))^2)/length(arr)
    aucMeanAbsoluteDeviation <- function(arr) sum(abs(arr - mean(arr)))/length(arr)

    function(iter, k, pert, ...){
        if (isStopped[k]) return(T)
        
        auc = CalcAUC(origS[[k]], pert)$area
        listAUC[[k]] <<- c(listAUC[[k]], auc)
        
        stop = F
        if (iter > iterMin){
            lengthAUC = length(listAUC[[k]])
            start = lengthAUC - 20
            if (start < 1)
                start = 1
            lastAUCs = listAUC[[k]][start:lengthAUC]
            
            if (aucMeanAbsoluteDeviation(lastAUCs) < madMin & aucMeanSquareDeviation(lastAUCs) < msdMin) {
                aucs = as.numeric(lapply(Filter(Negate(is.null), listAUC[1:kMax != k]), function(l) l[[length(l)]]))
                stop = length(which((abs(aucs - auc) < madMin) == T)) == 0
            }
        }
        
        isStopped[k] <<- stop
        
        stop
    }
}