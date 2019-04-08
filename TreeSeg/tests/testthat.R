library(treeSeg)
library(ape)

n <- 2 * 50
tree <- rtree(n, rooted = TRUE)
tree <- as.phylo(tree)

aktiveNode <- n + 50
offAktiveNode <- getOffspringTip(aktiveNode,tree)
prob1 <- 0.1
prob2 <- 0.75
p <- prob1 + prob2 * is.element(1:n, offAktiveNode)

yBin <- rbinom(n, 1, p)
ansBin <- treeSeg(yBin, tree, alpha = 0.1)
ansBin$numbAN
ansBinTest <- treeTest(yBin, tree, alpha = 0.1)
ansBinTest <- treeTest(yBin)


yNormal <- rnorm(n, p, sd = 0.1)
ansNormal <- treeSeg(yNormal, tree, alpha = 0.1, fam = "gauss")
ansNormal$numbAN
ansNormalTest <- treeTest(yNormal, tree, alpha = 0.1, fam = "gauss")
ansNormalTest <- treeTest(yNormal, fam = "gauss")


#add a second active node

aktiveNode2 <- n + 62
offAktiveNode2 <- getOffspringTip(aktiveNode2,tree)
p[is.element(1:n, offAktiveNode2)] <- 0.5

yBin <- rbinom(n, 1, p)
ansBin <- treeSeg(yBin, tree, alpha = 0.1)
ansBin$numbAN

yNormal <- rnorm(n, p, sd = 0.1)
ansNormal <- treeSeg(yNormal, tree, alpha = 0.1, fam = "gauss")
ansNormal$numbAN

#add missing data
yBin[sample(1:n, 10)] <- NA
ansBin <- treeSeg(yBin, tree, alpha = 0.1)
ansBin$numbAN

yNormal[sample(1:n, 10)] <- NA
ansNormal <- treeSeg(yNormal, tree, alpha = 0.1, fam = "gauss")
ansNormal$numbAN




