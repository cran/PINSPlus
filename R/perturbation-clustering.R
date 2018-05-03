#' @title Perturbation clustering
#' @description Perform subtyping using one type of high-dimensional data
#' 
#' @param data Input matrix or data frame. The rows represent items while the columns represent features.
#' @param kMax The maximum number of clusters. The algorithm runs from \code{k = 2} to \code{k = kMax}. Default value is \code{5}.
#' @param verbose Boolean value indicating the algorithm to run with or without logging. Default value is \code{TRUE}.
#' @param ncore Number of cores that the algorithm should use. Default value is \code{2}.
#' @param clusteringMethod The name of built-in clustering algorithm that PerturbationClustering will use. Currently supported algorithm are \code{kmeans}, \code{pam} and \code{hclust}. Default value is "\code{kmeans}".
#' @param clusteringOptions A list of parameter will be passed to the clustering algorithm in \code{clusteringMethod}.
#' @param clusteringFunction The clustering algorithm function that will be used instead of built-in algorithms.
#' @param perturbMethod The name of built-in perturbation method that PerturbationClustering will use, currently supported methods are \code{noise} and \code{subsampling}. Default value is "\code{noise}".
#' @param perturbOptions A list of parameter will be passed to the perturbation method in \code{perturbMethod}.
#' @param perturbFunction The perturbation method function that will be used instead of built-in ones.
#' @param iterMin The minimum number of iterations. Default value is \code{20}.
#' @param iterMax The maximum number of iterations. Default value is \code{200}.
#' @param madMin The minimum of Mean Absolute Deviation of \code{AUC} of Connectivity matrix for each \code{k}. Default value is \code{1e-03}.
#' @param msdMin The minimum of Mean Square Deviation of \code{AUC} of Connectivity matrix for each \code{k}. Default value is \code{1e-06}.
#
#' @details
#' PerturbationClustering implements the Perturbation Clustering algorithm of Nguyen, et al (2017).
#' It aims to determine the optimum cluster number and location of each sample in the clusters in an unsupervised analysis.
#' 
#' PerturbationClustering takes input as a numerical matrix or data frame of items as rows and features as columns.
#' It uses a clustering algorithm as the based algorithm.
#' Current built-in algorithms that users can use directly are \code{kmeans}, \code{pam} and \code{hclust}.
#' The default parameters for built-in \code{kmeans} are \code{nstart = 20 and iter.max = 1000}.
#' Users can change the parameters of built-in clustering algorithm by passing the value into \code{clusteringOptions}.
#' 
#' PerturbationClustering also allows users to pass their own clustering algorithm instead of using built-in ones by using \code{clusteringFunction} parameter. 
#' Once \code{clusteringFunction} is specified, \code{clusteringMethod} will be skipped.
#' The value of \code{clusteringFunction} must be a function that takes two arguments: \code{data} and \code{k}, 
#' where \code{data} is a numeric matrix or data frame containing data that need to be clustered, and \code{k} is the number of clusters.
#' \code{clusteringFunction} must return a vector of labels indicating the cluster to which each sample is allocated.
#' 
#' PerturbationClustering uses a perturbation method to perturb clustering input data.
#' There are two built-in methods are \code{noise} and \code{subsampling} that users can use directly by passing to \code{perturbMethod} parameter.
#' Users can change the default value of built-in perturbation methods by passing new value into \code{perturbOptions}:
#' 
#' 1. \code{noise} perturbation method takes two arguments: \code{noise} and \code{noisePercent}. The default values are \code{noise = NULL and noisePercent = "median"}.
#' If \code{noise} is specified. \code{noisePercent} will be skipped.\cr
#' 2. \code{subsampling} perturbation method takes one argument \code{percent} which has default value of \code{80}
#' 
#' Users can also use their own perturbation methods by passing them into \code{perturbFunction}. 
#' Once \code{perturbFunction} is specified, \code{perturbMethod} will be skipped.
#' The value of \code{perturbFunction} must be a function that takes one argument \code{data}
#' - a numeric matrix or data frame containing data that need to be perturbed.
#' \code{perturbFunction} must return an object list which is as follows:
#' 
#' 1. \code{data}: the perturbed data\cr
#' 2. \code{ConnectivityMatrixHandler}: a function that takes three arguments:
#' \code{connectivityMatrix} - the connectivity matrix generated after clustering returned \code{data}, 
#' \code{iter} - the current iteration and \code{k} - the number of cluster. 
#' This function must return a compatible connectivity matrix with the original connectivity matrix. 
#' This function aims to correct the \code{connectivityMatrix} if needed and returns the corrected version of it.\cr
#' 3. \code{MergeConnectivityMatrices}: a function that takes four arguments: \code{oldMatrix}, \code{newMatrix}, \code{k} and \code{iter}. 
#' The \code{oldMatrix} and \code{newMatrix} are two connectivity matrices that need to be merged,
#' \code{k} is the cluster number and \code{iter} is the current number of iteration.
#' This function must returns a connectivity matrix that is merged from \code{oldMatrix} and \code{newMatrix}
#' 
#' @return
#' \code{PerturbationClustering} returns a list with at least the following components:
#' \item{k}{The optimal number of clusters}
#' \item{cluster}{A vector of labels indicating the cluster to which each sample is allocated}
#' \item{origS}{A list of original connectivity matrices}
#' \item{pertS}{A list of perturbed connectivity matrices}
#' 
#' 
#' @references
#' 
#' 1. T Nguyen, R Tagett, D Diaz, S Draghici. A novel method for data integration and disease subtyping. Genome Research, 27(12):2025-2039, 2017.
#' 
#' 2. T. Nguyen, "Horizontal and vertical integration of bio-molecular data", PhD thesis, Wayne State University, 2017.
#' 
#' @seealso \code{\link{kmeans}}, \code{\link{pam}}
#' 
#' @examples
#' \donttest{
#' # Load the dataset AML2004
#' data(AML2004)
#' 
#' # Perform the clustering
#' result <- PerturbationClustering(data = AML2004$Gene)
#' 
#' # Plot the result
#' condition = seq(unique(AML2004$Group[, 2]))
#' names(condition) <- unique(AML2004$Group[, 2])
#' plot(
#'     prcomp(AML2004$Gene)$x,
#'     col = result$cluster,
#'     pch = condition[AML2004$Group[, 2]],
#'     main = "AML2004"
#' )
#' legend(
#'     "bottomright",
#'     legend = paste("Cluster ", sort(unique(result$cluster)), sep = ""),
#'     fill = sort(unique(result$cluster))
#' )
#' legend("bottomleft", legend = names(condition), pch = condition)
#' 
#' # Change kmeans parameters
#' result <- PerturbationClustering(
#'     data = AML2004$Gene,
#'     clusteringMethod = "kmeans",
#'     clusteringOptions = list(
#'         iter.max = 500,
#'         nstart = 50
#'     )
#' )
#' 
#' # Change to use pam
#' result <- PerturbationClustering(data = AML2004$Gene, clusteringMethod = "pam")
#' 
#' # Change to use hclust
#' result <- PerturbationClustering(data = AML2004$Gene, clusteringMethod = "hclust")
#' 
#' # Pass a user-defined clustering algorithm
#' result <- PerturbationClustering(data = AML2004$Gene, clusteringFunction = function(data, k){
#'     # this function must return a vector of cluster
#'     kmeans(x = data, centers = k, nstart = k*10, iter.max = 2000)$cluster
#' })      
#' 
#' # Use noise as the perturb method
#' result <- PerturbationClustering(data = AML2004$Gene, 
#'                                  perturbMethod = "noise", 
#'                                  perturbOptions = list(noise = 0.3))
#' # or
#' result <- PerturbationClustering(data = AML2004$Gene, 
#'                                  perturbMethod = "noise", 
#'                                  perturbOptions = list(noisePercent = 10))
#' 
#' # Change to use subsampling
#' result <- PerturbationClustering(data = AML2004$Gene, 
#'                                  perturbMethod = "subsampling", 
#'                                  perturbOptions = list(percent = 90))
#'
#' # Users can pass their own perturb method
#' result <- PerturbationClustering(data = AML2004$Gene, perturbFunction = function(data){
#'    rowNum <- nrow(data)
#'    colNum <- ncol(data)
#'    epsilon <-
#'        matrix(
#'            data = rnorm(rowNum * colNum, mean = 0, sd = 1.234),
#'            nrow = rowNum,
#'            ncol = colNum
#'        )
#'    
#'    list(
#'        data = data + epsilon,
#'        ConnectivityMatrixHandler = function(connectivityMatrix, ...) {
#'            connectivityMatrix
#'        },
#'        MergeConnectivityMatrices = function(oldMatrix, newMatrix, iter, ...){
#'            return((oldMatrix*(iter-1) + newMatrix)/iter)
#'        }
#'    )
#' })
#' }
#' @import stats utils cluster entropy parallel pbmcapply doParallel
#' @importFrom foreach foreach %dopar%
#' @export
PerturbationClustering <- function(data, kMax = 5, verbose = T, ncore = 2, # Algorithm args
                                   clusteringMethod = "kmeans", clusteringFunction = NULL, clusteringOptions = NULL, # Based clustering algorithm args
                                   perturbMethod = "noise", perturbFunction = NULL, perturbOptions = NULL, # Perturbed function args
                                   iterMin = 20, iterMax = 200, madMin = 1e-03, msdMin = 1e-06# Stopping condition for generating perturbation matrix
) {
    now = Sys.time()
    
    # defined log function
    log <- if(!verbose) function(...){} else function(...){
        message(...)
        flush.console()
    }
    
    clusteringAlgorithm = GetClusteringAlgorithm(clusteringMethod = clusteringMethod, clusteringFunction = clusteringFunction, clusteringOptions = clusteringOptions)
    perturbationAlgorithm = GetPerturbationAlgorithm(data = data, perturbMethod = perturbMethod, perturbFunction = perturbFunction, perturbOptions = perturbOptions)
    
    log("Clustering method: ", clusteringAlgorithm$name)
    log("Perturbation method: ", perturbationAlgorithm$name)

    seed = round(rnorm(1)*10^6)
    pca <- prcomp(data, rank. = min(nrow(data), 200))

    # get the partitioning from simply clustering the real data, for consecutive k
    log("Building original connectivity matrices")
    set.seed(seed)
    origPartition <- GetOriginalSimilarity(data = pca$x, clusRange = 2 : kMax, clusteringAlgorithm = clusteringAlgorithm$fun, showProgress = verbose, ncore)
    origS <- origPartition$origS

    listAUC <- list()
    for (k in 2:kMax){
        listAUC[k] <- list(c()) 
    }
    # declare stopping criteria for Perturbed Similarity computing
    stoppingCriteriaHandler <- DiffDevianceStoppingHandler(kMax = kMax, origS = origS, iterMin = iterMin, madMin = madMin, msdMin = msdMin, onExcute = function(k, AUCs){
        listAUC[[k]] <<- AUCs
    })

    log("Building perturbed connectivity matrices")
    set.seed(seed)
    pertS <- GetPerturbedSimilarity(data = pca$x, clusRange = 2 : kMax, iterMax = iterMax, iterMin = iterMin, origS = origS,
                                    clusteringAlgorithm = clusteringAlgorithm$fun, perturbedFunction = perturbationAlgorithm$fun,
                                    stoppingCriteriaHandler = stoppingCriteriaHandler,
                                    showProgress = verbose, ncore = ncore)
    
    for (k in 2:kMax){
        AUCs <- listAUC[[k]]
        pert <- matrix(0, nrow(data), nrow(data))
        for (i in 1:length(AUCs)){
            pert <- pert + pertS[[k]][[i]]*length(which(AUCs == AUCs[i]))
        }
        pertS[[k]] <- pert/sum(as.matrix(table(AUCs))[,1]^2)
    }

    # get discrepancy message('Calculate discrepancy between original and perturbed connectivity matrices')
    Discrepancy <- CalcPerturbedDiscrepancy(origS, pertS, clusRange = 2 : kMax)
    
    Discrepancy$AUC = round(Discrepancy$AUC, digits = 4)#meanAUCs
    clus <- min(which(Discrepancy$AUC == max(Discrepancy$AUC[2:kMax])))
    
    timediff = Sys.time() - now;
    log("Done in ", timediff, " ", units(timediff), ".\n")
    
    list(k = clus, cluster = origPartition$groupings[[clus]], origS = origS, pertS = pertS, Discrepancy = Discrepancy)
}