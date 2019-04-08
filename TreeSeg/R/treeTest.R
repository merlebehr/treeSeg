#' Multiscale test for association between tips and tree structure
#' 
#' For a given tree with observations at the tips, 
#' tests for association between the tree structure and the observational distribution.
#' The method is based on a multiscale approach. The null hypothesis assumes no association, which corresponds to
#' i.i.d. tree-ordered observations. 
#' In general, this is computationally faster than deriving the full \code{\link{treeSeg}} segmentation.
#' 
#' @param y observations, binary if \code{fam} equals "binomial" and numeric if \code{fam} equals "gauss". Missing data should be NA's and is silently ignored.
#' @param tree rooted tree object with tip labels \code{1:length(y)} of class phylo. 
#' The tip labels are assumed to be of consecutive order, such that for any inner node the offspring tips are of the form i, i+1, ..., j-1, j.
#' If missing, it is assumed that \code{y} is already ordered with respect to the tree-tips.
#' @param q threshold parameter corresponding to \code{1-alpha} quantile of multiscale statistic 
#' @param alpha confidence level (if specified \code{q} is silently ignored)
#' @param fam specifies distribution of data and can take either of the values "gauss" or "binomial" (default) where the letter assumes n=1, i.e. Bernoulli distribution
#' @param tipOrder specifies ordering of tips and can take either of the values "cladewise" (default) or "unchanged". 
#' The former applies \code{ape::reorder.phylo(tree)} prior to analysis, to make outcome independent of tip ordering.
#' @return A single integer number, where 0 corresponds to acceptance of null hypothesis at given level \code{alpha}. 
#' Any integer lager than 0 corresponds to rejection of null hypothesis at given level \code{alpha} 
#' and gives a lower bound on the minimal number of active nodes in the full \code{\link{treeSeg}} segmentation.
#' 
#' @examples
#' library(ape)
#' n <- 2 * 50
#' tree <- rtree(n, rooted = TRUE)
#' aktiveNode <- n + 50
#' offAktiveNode <- getOffspringTip(aktiveNode,tree)
#' prob1 <- 0.1
#' prob2 <- 0.75
#' p <- prob1 + prob2 * is.element(1:n, offAktiveNode)
#' y <- rbinom(n, 1, p)
#' ans <- treeTest(y, tree, alpha = 0.1)
#' yGauss <- rnorm(n, p, sd = 0.1)
#' ansGauss <- treeTest(yGauss, tree, alpha = 0.1, fam = "gauss")
#' @references 
#' Frick, K., Munk, A., Sieling, H. (2014) 
#' Multiscale change-point inference. 
#' With discussion and rejoinder by the authors. 
#' Journal of the Royal Statistical Society, Series B 76(3), 495:580.
#' @note If neither \code{q} nor \code{alpha} is specified, \code{alpha = 0.1} is selected. 
#' @seealso \code{\link{treeSeg}}
#' @export
#' @import stepR

treeTest<- function(y, tree, q, alpha, fam, tipOrder){
  
  if(!missing(tree)){
    if(missing(tipOrder)){
      tipOrder = "cladewise"
    }
    
    ### check if tips are of consecutive order
    for(i in (length(tree$tip.label)+1):(length(tree$tip.label) + tree$Nnode)){
      offspringLength <- length(getOffspringTip(i, tree))
      indexRange <- getOffspringTipB(i, tree)
      indexLength <- max(indexRange) - min(indexRange) + 1
      if(offspringLength != indexLength){
        warning("Indexes of tip labels are not of consecutive order.") 
      }
    }
    
    if(tipOrder == "cladewise"){
      #### reordering the tips
      tree <- ape::reorder.phylo(tree)
      orderT <- tree$edge[tree$edge[,2] <= length(tree$tip.label)  ,2]
      edgeVec <- as.numeric(tree$edge)
      ind <- lapply(orderT,function(x) which(edgeVec == x))
      for(i in 1:length(ind)){
        edgeVec[ind[[i]]] <- i
      }
      tree$edge <- matrix(edgeVec, ncol = 2)
      y <- y[orderT]
    }
    
    
    #### check if there is missing data
    missingData <- !all(!is.na(y))
    if(missingData){
      indMiss <- which(is.na(y))
      orderTips <- which(!is.na(y))
      treeFull <- tree
      yFull <- y
      tree <- ape::drop.tip(tree, indMiss, collapse.singles = TRUE)
      y <- y[-indMiss]
      if(is.null(treeFull$ancM)){
        treeFullNP <- treeFull
        treeFull <- prepTree(treeFull)
      }
    }
    
    #### prepare tree
    if(is.null(tree$ancM)){
      tree <- prepTree(tree)
    }
    
  }
  
  if(missing(fam)){
    fam = "binomial"
  }
  
  if(fam == "binomial"){
    if(!all(is.element(y, c(0,1)))){
      warning("y must be binomial vector for fam = binomial")
    }
  }
  
  
  if(missing(q)){
    if(missing(alpha)){
      alpha = 0.1
    }
    q <- stats::quantile(stepR::MRC.asymptotic, 1-alpha)
  }
  
  #remove missing data
  if(any(is.na(y))){
    ind <- is.na(y)
    y <- y[-ind]
  }

  #compute smuce estimate
  if(fam == "binomial"){
    ans <- length(smuceR(y = y, q = q, lengths = 1:length(y), family = "binomial", param = 1)$value)
  }else{
    ans <- length(smuceR(y = y, q = q, lengths = 1:length(y))$value)
  }
  
  #compute minimal number of active nodes
  ans <- ceiling((ans - 1)/2)

  return(ans)
}