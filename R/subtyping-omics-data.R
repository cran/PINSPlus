#' @title Subtyping multi-omics data
#' @description Perform subtyping using multiple types of data
#' 
#' @param dataList a list of data matrices or data frames. Each matrix represents a data type where the rows are items and the columns are features. The matrices must have the same set of items.
#' @param kMax the maximum number of clusters. Default value is \code{5}.
#' @param agreementCutoff agreement threshold to be considered consistent. Default value is \code{0.5}.
#' @param verbose set it to \code{TRUE} of \code{FALSE} to get more or less details respectively.
#' @param ... these arguments will be passed to \code{PerturbationClustering} algorithm. See details for more information
#' 
#' @details
#' 
#' \code{SubtypingOmicsData} implements the Subtyping multi-omic data that are based on Perturbaion clustering algorithm of Nguyen, et al (2017).
#' The input is  a list of data matrices where each matrix represents the molecular measurements of a data type. The input matrices must have the same number of rows. 
#' \code{SubtypingOmicsData} aims to find the optimum number of subtypes and location of each sample in the clusters from integrated input data \code{dataList} through two processing stages:
#' 
#' 1. Stage I: The algorithm first partitions each data type using the function \code{PerturbationClustering}.
#' It then merges the connectivities across data types into similarity matrices.
#' Both kmeans and similarity-based clustering algorithms - partitioning around medoids \code{pam} are used to partition the built similarity.
#' The algorithm returns the partitioning that agrees the most with individual data types.\cr
#' 2. Stage II: The algorithm attempts to split each discovered group if there is a strong agreement between data types,
#' or if the subtyping in Stage I is very unbalanced.
#'
#' @return
#' 
#' \code{SubtypingOmicsData} returns a list with at least the following components:
#' \item{cluster1}{A vector of labels indicating the cluster to which each sample is allocated in Stage I}
#' \item{cluster2}{A vector of labels indicating the cluster to which each sample is allocated in Stage II}
#' \item{dataTypeResult}{A list of results for individual data type. Each element of the list is the result of the \code{PerturbationClustering} for the corresponding data matrix provided in dataList.}
#' 
#' 
#' @references
#' 
#' 1. T Nguyen, R Tagett, D Diaz, S Draghici. A novel method for data integration and disease subtyping. Genome Research, 27(12):2025-2039, 2017.
#' 
#' 2. T. Nguyen, "Horizontal and vertical integration of bio-molecular data", PhD thesis, Wayne State University, 2017.
#' 
#' @seealso \code{\link{PerturbationClustering}}
#' 
#' @examples
#' \donttest{
#' # Load the kidney cancer carcinoma data
#' data(KIRC)
#' 
#' # Perform subtyping on the multi-omics data
#' dataList <- list (KIRC$GE, KIRC$ME, KIRC$MI) 
#' names(dataList) <- c("GE", "ME", "MI")
#' result <- SubtypingOmicsData(dataList = dataList)
#' 
#' # Change Pertubation clustering algorithm's arguments
#' result <- SubtypingOmicsData(
#'     dataList = dataList, 
#'     clusteringMethod = "kmeans", 
#'     clusteringOptions = list(nstart = 50)
#' )
#'
#' # Plot the Kaplan-Meier curves and calculate Cox p-value
#' library(survival)
#' cluster1=result$cluster1;cluster2=result$cluster2
#' a <- intersect(unique(cluster2), unique(cluster1))
#' names(a) <- intersect(unique(cluster2), unique(cluster1))
#' a[setdiff(unique(cluster2), unique(cluster1))] <- seq(setdiff(unique(cluster2), unique(cluster1))) 
#'                                                   + max(cluster1)
#' colors <- a[levels(factor(cluster2))]
#' coxFit <- coxph(
#'  Surv(time = Survival, event = Death) ~ as.factor(cluster2),
#'  data = KIRC$survival,
#'  ties = "exact"
#' )
#' mfit <- survfit(Surv(Survival, Death == 1) ~ as.factor(cluster2), data = KIRC$survival)
#' plot(
#'  mfit, col = colors,
#'  main = "Survival curves for KIRC, level 2",
#'  xlab = "Days", ylab = "Survival",lwd = 2
#' )
#' legend("bottomright", 
#'     legend = paste(
#'         "Cox p-value:", 
#'         round(summary(coxFit)$sctest[3], digits = 5), 
#'         sep = ""
#'     )
#' )
#' legend(
#'     "bottomleft",
#'     fill = colors,
#'     legend = paste(
#'         "Group ",
#'         levels(factor(cluster2)),": ", table(cluster2)[levels(factor(cluster2))], 
#'         sep =""
#'     )
#' )
#' 
#' }
#' @export
SubtypingOmicsData <- function (dataList, kMax = 5, agreementCutoff = 0.5, verbose = T, ...) {
    now = Sys.time()
    
    # defined log function
    mlog <- if(!verbose) function(...){} else function(...){
        message(...)
        flush.console()
    }
    
    seed = round(rnorm(1)*10^6)
    
    runPerturbationClustering <- function(dataList, kMax, stage = 1){
        dataTypeResult <- lapply(dataList, function(data) {
            set.seed(seed)
            PerturbationClustering(data, kMax, verbose = verbose,...)
        })
        origList <- lapply(dataTypeResult, function(r) r$origS[[r$k]])
        orig = Reduce('+', origList)/length(origList)
        PW = Reduce('*', origList)
        agreement = (sum(orig == 0) + sum(orig == 1) - nrow(orig)) / (nrow(orig) ^ 2 - nrow(orig))
        
        pert =  Reduce('+', lapply(dataTypeResult, function(r) r$pertS[[r$k]]))/length(dataList)
        
        groups <- NULL
        
        mlog("STAGE : ", stage, "\t Agreement : ", agreement)
        
        if (agreement >= agreementCutoff){
            hcW <- hclust(dist(PW))
            maxK = min(kMax*2, dim(unique(PW, MARGIN = 2))[2] - (stage - 1))
            maxHeight = FindMaxHeight(hcW, maxK = min(2*maxK, 10))
            groups <- cutree(hcW, maxHeight)
        }
        
        list(dataTypeResult = dataTypeResult, orig = orig, pert = pert, PW = PW, groups = groups, agreement = agreement)
    }
    
    pResult <- runPerturbationClustering(dataList, kMax)
    
    groups <- pResult$groups
    groups2 <- NULL
    
    if (!is.null(groups)) {
        groups2 <- groups
        
        for (g in sort(unique(groups))) {
            miniGroup <- names(groups[groups == g])
            if (length(miniGroup) > 30) {
                groupsM <- runPerturbationClustering(dataList = lapply(dataList, function(d) d[miniGroup, ]), kMax = min(kMax, 5), stage = 2)$groups
                if (!is.null(groupsM))
                    groups2[miniGroup] <- paste(g, groupsM, sep = "-")
            }
        }
    }
    else{
        set.seed(seed)
        
        orig <- pResult$orig
        dataTypeResult <- pResult$dataTypeResult
        clusteringAlgorithm = GetClusteringAlgorithm(...)$fun
        
        groupings <- lapply(dataTypeResult, function(r) clusteringAlgorithm(data = r$origS[[r$k]], k = r$k))
        
        pGroups <- ClusterUsingPAM(orig = orig, kMax = kMax*2, groupings = groupings)
        hGroups <- ClusterUsingHierarchical(orig = orig, kMax = kMax*2, groupings = groupings)
        
        pAgree  = pGroups$agree; hAgree  = hGroups$agree;
        
        groups <- (if (pAgree > hAgree) pGroups else if (hAgree > pAgree) hGroups else {
            pAgree = ClusterUsingPAM(orig = pResult$pert, kMax = kMax, groupings = groupings)$agree
            hAgree = ClusterUsingHierarchical(orig = pResult$pert, kMax = kMax, groupings = groupings)$agree
            if (hAgree - pAgree >= 1e-3) hGroups else pGroups
        })$cluster
        
        names(groups) <- rownames(orig)
        
        mlog("Check if can proceed to stage II")
        groups2 <- groups
        normalizedEntropy = entropy::entropy(table(groups)) / log(length(unique(groups)), exp(1))
        
        if (normalizedEntropy < 0.5) {
            for (g in sort(unique(groups))) {
                miniGroup <- names(groups[groups == g])
                #this is just to make sure we don't split a group that is already very small
                if (length(miniGroup) > 30) {
                    #this is to check if the data types in this group can be split
                    groupsM <- runPerturbationClustering(dataList = lapply(dataList, function(d) d[miniGroup, ]), kMax = min(kMax, 5), stage = 2)$groups
                    if (!is.null(groupsM))
                        groups2[miniGroup] <- paste(g, groupsM, sep = "-")
                }
            }
        }
    }
    
    timediff = Sys.time() - now;
    mlog("Done in ", timediff, " ", units(timediff), ".\n")
    
    list(
        cluster1 = groups,
        cluster2 = groups2,
        dataTypeResult = pResult$dataTypeResult
    )
    
}