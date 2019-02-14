prepResult <- function(ans, y, tree){
  if(length(tree$ancM) == 0){
    tree <- prepTree(tree)
  }
  ancM <- tree$ancM
  n <- length(tree$tip.label)
  comb <- ans$comb
  optcost <- ans$optcost
  root <- ans$root
  minI <- ans$minI[root]
  
  if(minI > length(comb[root][[1]][[1]]$comb)){
    mlSolNodes <- unique(get.ml.sol(comb[root][[1]][[1]]$comb, comb, minI))
  }else{
    mlSolNodes <- unique(as.numeric(unlist(comb[root][[1]][[1]]$comb)))
  }
  confSet <- sort(get.conf(unlist(sapply(comb[root][[1]], function(x) x$comb)), comb))
  confBand <- get.conf.band(root, cbind(rep(1,n), rep(0,n)), comb, tree)
  confBand <- unname(confBand)
  
  mlSol <- apply(ancM, 2, function(x) sum(is.element(mlSolNodes, x) * 2^(0:(minI-1)))  )
  mlSolP <- numeric(n)
  for(i in 1 : length(unique(mlSol))){
    ind <- which(mlSol == unique(mlSol)[i])
    #indB <- which(apply(bou, 1, function(x) length(setdiff(x[1] : x[2], ind)) == 0))
    #mlSolP[ind] <- max(min(mean(y[ind]), min(bou[indB,4])), max(bou[indB,3]))
    mlSolP[ind] <- max(min(mean(y[ind]), min(confBand[ind,2])), max(confBand[ind,1]))
  }
  return(list(numbAN = minI, mlAN = mlSolNodes, mlP = mlSolP, confSetAN = confSet, confBandP = confBand))
}


#function to extract the maximum likelihood solution from the output of 
# seg Tree [node] is the set of active nodes which is clostest to the root
# comb[x][[1]] are all multiscale solutions of inner node [x]
# and comb[x][[1]][[1]] is the solution with minimal cost
# maxN is the final number of active nodes
get.ml.sol <- function(node, comb, maxN){
  new.node <- unlist(sapply(node, function(x) comb[x][[1]][[1]]$comb))
  if(length(c(node, new.node)) < maxN){
    nnew.node <- get.ml.sol(new.node, comb, maxN - length(node))
    return(c(node, new.node, nnew.node))
  }else{
    return(c(node, new.node))
  }
}

# function to extract the confidence sets from the output of segTree
# it seccessively combines the different entries of comb
get.conf <- function(node, comb){
  new.node <- unique(unlist(sapply(node, function(x) 
    sapply(comb[x][[1]], function(y) y$comb))))
  if(length(new.node) > 0){
    nnew.node <- get.conf(new.node, comb)
    return(unique(c(node, new.node, nnew.node)))
  }else{
    return(unique(c(node, new.node)))
  }
}


#extracts the confidence bands from the output of segTree
get.conf.band <- function(node, p, comb, tree){
  for(j in 1:length(comb[node][[1]])){
    p[offspringDiffTip(node, comb[node][[1]][[j]]$comb, tree), 1] <- 
      pmin(p[offspringDiffTip(node, comb[node][[1]][[j]]$comb, tree), 1], 
           comb[node][[1]][[j]]$maxB)
    p[offspringDiffTip(node, comb[node][[1]][[j]]$comb, tree), 2] <- 
      pmax(p[offspringDiffTip(node, comb[node][[1]][[j]]$comb, tree), 2], 
           comb[node][[1]][[j]]$minB)
  }
  
  new.node <- unlist(sapply(comb[node][[1]], function(x) x$comb))
  
  if(length(new.node)> 0){
    P <- lapply(new.node, function(x) get.conf.band(x, p, comb, tree))
    plower <- apply( sapply(P, function(x) x[,1]), 1, min )
    pupper <- apply( sapply(P, function(x) x[,2]), 1, max )
    return(cbind(plower, pupper))
  }else{
    return(p)
  }
}

