library("seventyGeneData")
library(treeSeg)
library(ape)
library(stepR)
library(treeio)
library(data.tree)
library(Biobase)

# load the microarray data from the "seventyGeneData" package.
# --------------------
data(vantVeer)
our_pheno_data = pData(vantVeer)

# extract the primary breast cancers
# -------------------------
our_dataset = vantVeer[,vantVeer$DataSetType != '19samples.gz']

# extract genes that were used in the original paper.
# ---------------------------------------------
pvals = assayData(our_dataset)$pValue
temp_data = Biobase::exprs(our_dataset)
idx = rowSums(temp_data > log10(2) | temp_data <= log10(0.5) ,na.rm=T)>5 & rowSums(pvals<0.01)>=5
our_genes = our_dataset[idx,]

# cluster the samples using correlation as distance metric
# -----------------------------
our_counts = Biobase::exprs(our_genes)

# let's get the heirachical tree
x = as.dist(1- cor((our_counts),use='pairwise.complete.obs'))
temp_tree = hclust(x,method='comp')

# convert tree to APE format.
tree = as.phylo(temp_tree)

# extract the phenotypes for the tree tips and order them correctly for the treeSeg
# -------------------------------

# let's get the order at which the tips will be plotted.
tree_data =  as.Node(tree)
tip_order = Get(tree_data$leaves, "name")

# let's get the order of the tips as they are plotted an internally so that we can be consistant for plotting and analysis.
pheno_idx = match(tip_order, our_genes$SampleName)
tree_idx = match(tip_order, tree$tip.label)

#rename indices of tree tips accordingly
tree$tip.label = tree$tip.label[tree_idx]
edgeVec <- as.numeric(tree$edge)
ind <- lapply(tree_idx,function(x) which(edgeVec == x))
for(i in 1:length(ind)){
  edgeVec[ind[[i]]] <- i
}
tree$edge <- matrix(edgeVec, ncol = 2)




# For each of the phenotypes of interest run treeSeg with alpha = 0.05.
# ---------------------------------------

# BRCA mutation phenotype
all(our_genes$SampleName[pheno_idx] == tip_order)
pheno = our_genes$Brca1.mutation[pheno_idx]
pheno[pheno == 'BRCA2'] = '1'
y_brca = as.numeric(pheno)
brca <- treeSeg(y_brca, tree, alpha = 0.05) #compute estimate

# ERP
all(our_genes$SampleName[pheno_idx] == tip_order)
pheno = our_genes$ERp[pheno_idx]
pheno[pheno >0] = 1
y_erp = pheno
erp <- treeSeg(y_erp, tree, alpha = 0.05) #compute estimate

# cancer grades
all(our_genes$SampleName[pheno_idx] == tip_order)
pheno = our_genes$grade[pheno_idx]
pheno[pheno <=2] = 0
pheno[pheno ==3] = 1
y_grade = pheno
grade <- treeSeg(y_grade, tree, alpha = 0.05) #compute estimate

# lymphatic infiltration
all(our_genes$SampleName[pheno_idx] == tip_order)
pheno = our_genes$Lymphocytic.Infiltrate[pheno_idx]
y_li = pheno
li <- treeSeg(y_li, tree, alpha = 0.05) #compute estimate

# angioinvasion
all(our_genes$SampleName[pheno_idx] == tip_order)
pheno = our_genes$Angioinvasion[pheno_idx]
y_angio = pheno
angio <- treeSeg(y_angio, tree, alpha = 0.05) #compute estimate

# metastasis
all(our_genes$SampleName[pheno_idx] == tip_order)
pheno = our_genes$metastases[pheno_idx]
y_meta = pheno
meta <- treeSeg(y_meta, tree, alpha = 0.05) #compute estimate


# make the plots for each of the phenotypes.
# --------------------------------
n = length(tree$tip.label)
lwdEdge <- rep(1.5, dim(tree$edge)[1])
layout(mat=matrix(1:18,ncol = 2),heights = rep(c(3,0.1, 1),6))
par(mar = c(0.1,4,0.1,0.1))

# BRCA mutation
# -----

# plot tree
plot(tree, type = "phylogram", show.tip.label =F , use.edge.length = F, 
     node.pos = 1, direction = "downwards", show.node.label = F, edge.width = lwdEdge) 

# plot maximum likelihood estimate for the nodes where the distribution of phenotype changes.
for(i in 1:length(brca$mlAN)){
  nodelabels("", brca$mlAN[i], pch = 18, col = "red", frame = "none",cex = 2)
}
text(80,70,'BRCA mutation')
text(-0,90,'a',cex=2)

# plot tip phenos
plot(1:n, rep(0,n), axes = F, col = c("lightgray", "black")[y_brca+1], pch = "|", cex = 2,xlab = '',ylab = '') 

# plot the maximum likelihood estimates for the phenotype distributions for each segment.
par(mar = c(2,4,0.1,0.1))
plot(1:n, brca$mlP,type = "l", lwd = 1, axes = F,ylim = c(-0.1,1.1),col='red',xlab = '',ylab = '')

# plot the confidence interval for the phenotype distributions.
polygon(c(1:n, rev(1:n)), c(brca$confBandP[,2], rev(brca$confBandP[,1])),
        col = rgb(255, 165, 0, alpha = 100, maxColorValue = 255), border = NA)

#axis(1, at = c(1,50,100,150,200), cex = 2)
axis(2, at = c(0,0.5,1),labels = c('0','0.5','1') )
mtext(side=2,'Prob.',line=2.5,cex=0.8)
points(1:n, rep(0,n), ylim = c(-0.1,1.1), type = "l", lwd = 1, axes = F)

# ERP plot.
# -------------
par(mar = c(0.1,4,0.1,0.1))

# plot tree
plot(tree, type = "phylogram", show.tip.label =F , use.edge.length = F, 
     node.pos = 1, direction = "downwards", show.node.label = F, edge.width = lwdEdge) 
for(i in 1:length(erp$mlAN)){
  nodelabels("", erp$mlAN[i], pch = 18, col = "red", frame = "none",cex = 2)
}
text(80,70,'ER expression')
text(-0,90,'b',cex=2)

# plot tip phenos
plot(1:n, rep(0,n), axes = F, col = c("lightgray", "black")[y_erp+1], pch = "|", cex = 2,xlab = '',ylab = '') 

# plot maximum likelihood estimate for the nodes where the distribution of phenotype changes.
par(mar = c(2,4,0.1,0.1))
plot(1:n, erp$mlP,type = "l", lwd = 1, axes = F,ylim = c(-0.1,1.1),col='red',xlab = '',ylab = '')

# plot the confidence interval for the phenotype distributions.
polygon(c(1:n, rev(1:n)), c(erp$confBandP[,2], rev(erp$confBandP[,1])),
        col = rgb(255, 165, 0, alpha = 100, maxColorValue = 255), border = NA)
axis(2, at = c(0,0.5,1) ,labels = c('0','0.5','1') )
mtext(side=2,'Prob.',line=2.5,cex=0.8)
points(1:n, rep(0,n), ylim = c(-0.1,1.1), type = "l", lwd = 1, axes = F)

# plot cancer grades
# --------------------
par(mar = c(0.1,4,0.1,0.1))

# plot tree
plot(tree, type = "phylogram", show.tip.label =F , use.edge.length = F, 
     node.pos = 1, direction = "downwards", show.node.label = F, edge.width = lwdEdge) 


for(i in 1:length(grade$mlAN)){
  nodelabels("", grade$mlAN[i], pch = 18, col = "red", frame = "none",cex = 2)
}
text(80,70,'Tumour grade > 2')
text(-0,90,'c',cex=2)

# plot tip phenos
plot(1:n, rep(0,n), axes = F, col = c("lightgray", "black")[y_grade+1], pch = "|", cex = 2,xlab = '',ylab = '') 

# plot the maximum likelihood estimates for the distributions.
par(mar = c(2,4,0.1,0.1))
plot(1:n, grade$mlP,type = "l", lwd = 1, axes = F,ylim = c(-0.1,1.1),col='red',xlab = '',ylab = '')


polygon(c(1:n, rev(1:n)), c(grade$confBandP[,2], rev(grade$confBandP[,1])),
        col = rgb(255, 165, 0, alpha = 100, maxColorValue = 255), border = NA)
axis(2, at = c(0,0.5,1),labels = c('0','0.5','1')  )
mtext(side=2,'Prob.',line=2.5,cex=0.8)
points(1:n, rep(0,n), ylim = c(-0.1,1.1), type = "l", lwd = 1, axes = F)

# plot lymphatic infiltration
# -------------------------
par(mar = c(0.1,4,0.1,0.1))

# plot tree
plot(tree, type = "phylogram", show.tip.label =F , use.edge.length = F, 
     node.pos = 1, direction = "downwards", show.node.label = F, edge.width = lwdEdge) 

for(i in 1:length(li$mlAN)){
  nodelabels("", li$mlAN[i], pch = 18, col = "red", frame = "none",cex = 2)
}
text(80,70,'Lymphocytic infiltrate')
text(-0,90,'d',cex=2)

# plot tip phenos
plot(1:n, rep(0,n), axes = F, col = c("lightgray", "black")[y_li+1], pch = "|", cex = 2,xlab = '',ylab = '') 

# plot the maximum likelihood estimates for the distributions.
par(mar = c(2,4,0.1,0.1))
plot(1:n, li$mlP,type = "l", lwd = 1, axes = F,ylim = c(-0.1,1.1),col='red',xlab = '',ylab = '')

polygon(c(1:n, rev(1:n)), c(li$confBandP[,2], rev(li$confBandP[,1])),
        col = rgb(255, 165, 0, alpha = 100, maxColorValue = 255), border = NA)
axis(2, at = c(0,0.5,1) ,labels = c('0','0.5','1') )
mtext(side=2,'Prob.',line=2.5,cex=0.8)
points(1:n, rep(0,n), ylim = c(-0.1,1.1), type = "l", lwd = 1, axes = F)


# plot angio invasion
# -----------
par(mar = c(0.1,4,0.1,0.1))

# plot tree
plot(tree, type = "phylogram", show.tip.label =F , use.edge.length = F, 
     node.pos = 1, direction = "downwards", show.node.label = F, edge.width = lwdEdge) 

for(i in 1:length(angio$mlAN)){
  nodelabels("", angio$mlAN[i], pch = 18, col = "red", frame = "none",cex = 2)
}
text(80,70,'Angioinvasion')
text(-0,90,'e',cex=2)

# plot tip phenos
plot(1:n, rep(0,n), axes = F, col = c("lightgray", "black")[y_angio+1], pch = "|",  cex = 2,xlab = '',ylab = '') 

# plot the maximum likelihood estimates for the distributions.
par(mar = c(2,4,0.1,0.1))
plot(1:n, angio$mlP,type = "l", lwd = 1, axes = F,ylim = c(-0.1,1.1),col='red',xlab = '',ylab = '')

polygon(c(1:n, rev(1:n)), c(angio$confBandP[,2], rev(angio$confBandP[,1])),
        col = rgb(255, 165, 0, alpha = 100, maxColorValue = 255), border = NA)
axis(2, at = c(0,0.5,1) ,labels = c('0','0.5','1') )
mtext(side=2,'Prob.',line=2.5,cex=0.8)
points(1:n, rep(0,n), ylim = c(-0.1,1.1), type = "l", lwd = 1, axes = F)

# plot metastasis
# -------------
par(mar = c(0.1,4,0.1,0.1))

# plot tree
plot(tree, type = "phylogram", show.tip.label =F , use.edge.length = F, 
     node.pos = 1, direction = "downwards", show.node.label = F, edge.width = lwdEdge) 

for(i in 1:length(meta$mlAN)){
  nodelabels("", meta$mlAN[i], pch = 18, col = "red", frame = "none",cex = 2)
}
text(80,70,'Metastasis status')
text(-0,90,'f',cex=2)

# plot tip phenos
plot(1:n, rep(0,n), axes = F, col = c("lightgray", "black")[y_meta+1], pch = "|", cex = 2,xlab = '',ylab = '') 

# plot the maximum likelihood estimates for the distributions.
par(mar = c(2,4,0.1,0.1))
plot(1:n, rep(0,n), ylim = c(-0.1,1.1), type = "l", lwd = 1, axes = F,xlab = '',ylab = '')

polygon(c(1:n, rev(1:n)), c(meta$confBandP[,2], rev(meta$confBandP[,1])),
        col = rgb(255, 165, 0, alpha = 100, maxColorValue = 255), border = NA)
points(1:n, meta$mlP,type = "l", lwd = 1, axes = F,ylim = c(-0.1,1.1),col='red')
axis(2, at = c(0,0.5,1) ,labels = c('0','0.5','1') )
mtext(side=2,'Prob.',line=2.5,cex=0.8)


