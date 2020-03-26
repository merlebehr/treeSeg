library(treeSeg)
library(ape)

set.seed(1)
n <- 2 * 50
tree <- rtree(n, rooted = TRUE)
tree <- as.phylo(tree)

aktiveNode <- n + 3
offAktiveNode <- getOffspringTip(aktiveNode,tree)
prob1 <- 0.01
prob2 <- 0.95
p <- prob1 + prob2 * is.element(1:n, offAktiveNode)

yBin <- rbinom(n, 1, p)
ansBin <- treeSeg(yBin, tree, alpha = 0.1)
ansBin$numbAN
ansBin <- treeSeg(yBin, tree, checkOrder = F)
ansBin$numbAN
ansBin <- treeSeg(yBin, tree, lengths = "dyadic")
ansBin$numbAN
ansBinTest <- treeTest(yBin, tree, alpha = 0.1)
ansBinTest
ansBinTest <- treeTest(yBin)
ansBinTest


yNormal <- rnorm(n, p, sd = 0.01)
ansNormal <- treeSeg(yNormal, tree, alpha = 0.1, fam = "gauss")
ansNormal$numbAN
ansNormal <- treeSeg(yNormal, tree, alpha = 0.1, fam = "gauss", checkOrder = F)
ansNormal$numbAN
ansNormal <- treeSeg(yNormal, tree, alpha = 0.1, fam = "gauss", lengths = "dyadic")
ansNormal$numbAN

ansNormalTest <- treeTest(yNormal, tree, alpha = 0.1, fam = "gauss")
ansNormalTest

ansNormalTest <- treeTest(yNormal, fam = "gauss")
ansNormalTest



#add a second active node

aktiveNode2 <- n + 62
offAktiveNode2 <- getOffspringTip(aktiveNode2,tree)
p[is.element(1:n, offAktiveNode2)] <- 0.5

yBin <- rbinom(n, 1, p)
ansBin <- treeSeg(yBin, tree, alpha = 0.1)
ansBin$numbAN


set.seed(99)
yNormal <- rnorm(n, p, sd = 0.001)
ansNormal <- treeSeg(yNormal, tree, alpha = 0.1, fam = "gauss")
ansNormal$numbAN



#add missing data
yBin[sample(1:n, 10)] <- NA
ansBin <- treeSeg(yBin, tree, alpha = 0.1)
ansBin <- treeSeg(yBin, tree, alpha = 0.1, lengths = "dyadic")
ansBin$numbAN

yNormal[sample(1:n, 10)] <- NA
ansNormal <- treeSeg(yNormal, tree, alpha = 0.1, fam = "gauss")
ansNormal <- treeSeg(yNormal, tree, alpha = 0.1, fam = "gauss", lengths = "dyadic")
ansNormal$numbAN




