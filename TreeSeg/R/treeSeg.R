#' Segmentation of tip-observations from tree structure
#' 
#' For a given tree with observations at the tips, 
#' performs segmentation which corresponds to a set of inner nodes in the tree. 
#' The method is based on a multiscale approach. This controls the probability of 
#' overestimating the number of active nodes at level alpha and yields confidence statements 
#' for the set of active nodes and the observational means (success probabilities for binary data).
#' 
#' @param y observations, binary if \code{fam} equals "binomial" and numeric if \code{fam} equals "gauss". 
#' Missing data should be NA's for which predictions are computed.
#' @param tree rooted tree object with tip labels \code{1:length(y)} of class phylo. 
#' The tip labels are assumed to be of consecutive order, such that for any inner node the offspring tips are of the form i, i+1, ..., j-1, j.
#' @param q threshold parameter corresponding to \code{1-alpha} quantile of multiscale statistic 
#' @param alpha confidence level (if specified \code{q} is silently ignored)
#' @param fam specifies distribution of data and can take either of the values "gauss" or "binomial" (default) where the letter assumes n=1, 
#' i.e. Bernoulli distribution
#' @param tipOrder specifies ordering of tips and can take either of the values "cladewise" (default) or "unchanged". 
#' The former applies \code{ape::reorder.phylo(tree)} prior to analysis, to make outcome independent of tip ordering.
#' @return A list with values \code{numbAN}, \code{mlAN}, \code{mlP}, \code{confSetAN}, 
#' and \code{confBandP}.
#' 
#' \item{numbAN}{single integer value which contains the number of estimated active nodes constituting a \code{1-alpha} lower confidence bound}
#' \item{mlAN}{vector of integer values containing the estimated active nodes (maximum likelihood solution under multiscale constraint)}
#' \item{mlP}{vector of numeric values between 0 and 1 with estimated success probabilities for tips when \code{fam = "binomial"} and vector of numeric values with estimated observational means for tips when \code{fam = "gauss"} }
#' \item{confSetAN}{vector of integer values constituting \code{1-alpha} confidence set for active nodes}
#' \item{confBandP}{matrix with two columns and number of rows equal to \code{length(y)} constituting \code{1-alpha} confidence band for success probabilites and mean values, respectively, of observations \code{y}}
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
#' ans <- treeSeg(y, tree, alpha = 0.1)
#' yGauss <- rnorm(n, p, sd = 0.1)
#' ansGauss <- treeSeg(yGauss, tree, alpha = 0.1, fam = "gauss")
#' @references 
#' Frick, K., Munk, A., Sieling, H. (2014) 
#' Multiscale change-point inference. 
#' With discussion and rejoinder by the authors. 
#' Journal of the Royal Statistical Society, Series B 76(3), 495:580.
#' @note If neither \code{q} nor \code{alpha} is specified, \code{alpha = 0.1} is selected.
#' @seealso \code{\link{getOffspringTip}}
#' @export
#' @import stepR

treeSeg<- function(y, tree, q, alpha, fam, tipOrder){
  
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
 
  if(missing(fam)){
    fam = "binomial"
  }
  
  if(fam == "binomial"){
    fam = 0
  }else{
    fam = 1
  }
  
  if(fam == 0){
    if(!all(is.element(y, c(0,1)))){
      warning("y must be binomial vector for fam = binomial")
    }
  }
    

  if(missing(q)){
    if(missing(alpha)){
      alpha = 0.1
    }
      ans <- segTree(y = y, tree = tree, alpha = alpha, fam = fam)
    }
  else{
    ans <- segTree(y = y, tree = tree, q = q, fam = fam)
  }
  
  if(length(ans) == 0){
    print("Problem size too large.")
    return(0)
  }
  
  #### prepare results in case of no change-point
  if(!is.null(ans$noAN)){
    numbAN = 0
    mlAN = numeric(0)
    mean = mean(y)
    if(mean > ans$upper){
      mlP = rep(ans$upper, length(y))
    }
    if(mean < ans$lower){
      mlP = rep(ans$lower, length(y))
    }
    if(mean >= ans$lower & mean <= ans$upper){
      mlP = rep(mean, length(y))
    }
    confSetAN = numeric(0)
    confBandP = cbind(rep(ans$lower, length(y)), rep(ans$upper, length(y)))
    
    return(list(numbAN = numbAN, mlAN = mlAN, mlP = mlP, confSetAN = confSetAN, confBandP = confBandP))
  }
  
  ans$minI <- ans$minI[-1]
  ans$comb <- ans$comb[-1]
  combCS <- ans$comb
  ans$optCost <- ans$optCost[-1]
  ans <- prepResult(ans, y, tree)
  
  
  #### modify results for missing data
  if(missingData){
    #transform mlAN
    for(i in 1:length(ans$mlAN)){
      off <- getOffspringTip(ans$mlAN[i], tree)
      if(length(off)>1){
        ans$mlAN[i] <- ape::getMRCA(treeFullNP, orderTips[off])
      }else{
        ans$mlAN[i] <- orderTips[off]
      }
    }
    #transform mlP
    mlPFull <- rep(NA, length(yFull))
    mlPFull[orderTips] <- ans$mlP
    
    #predict missing data
    ancMFull <- treeFull$ancM
    
    for(i in indMiss){
      anc <- ancMFull[,i]
      count <- length(anc)-1
      while(is.na(mlPFull[i])){
        siblings <- mlPFull[getOffspringTip(anc[count], treeFull)]
        if(all(is.na(siblings))){
          count <- count - 1
        }else{
          mlPFull[i] <- siblings[ which(!is.na(siblings))[1] ]
        }
      }
    }
    ans$mlP <- mlPFull
    
    #transform confSetAN
    for(i in 1:length(ans$confSetAN)){
      off <- getOffspringTip(ans$confSetAN[i],tree)
      if(length(off)>1){
        ans$confSetAN[i] <- ape::getMRCA(treeFullNP, orderTips[off])
      }else{
        ans$confSetAN[i] <- orderTips[off]
      }
    }
    
    #transform confBandP
    confBandPFull <- matrix(NA, ncol = 2, nrow = length(yFull))
    confBandPFull[orderTips, ] <- ans$confBandP
    confBandPFull[indMiss, ] <- matrix(rep(c(0,1), each = length(indMiss)), ncol = 2)
    ans$confBandP <- confBandPFull
  }
  
  
  if(tipOrder == "cladewise"){
    
  #modify results for original ordering of the tips
  #transform mlAN
  for(i in 1:length(ans$mlAN)){
    if(ans$mlAN[i] <= length(tree$tip.label)){
      ans$mlAN[i] <- orderT[ans$mlAN[i]]
    }
  }
  
  #transform mlP
  ans$mlP[orderT] <- ans$mlP
  
  #transform confSetAN
  ind <- which(ans$confSetAN <= length(tree$tip.label))
  ans$confSetAN[ind] <- orderT[ans$confSetAN[ind]]
  
  #transform confBandP
  ans$confBandP[orderT, ] <- ans$confBandP
  
  }
  
  #attr(ans, "combCS") <- combCS
  return(ans)
}