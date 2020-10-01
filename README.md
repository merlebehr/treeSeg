# treeSeg
Tree structures, showing hierarchical relationships and the latent structures between samples, are ubiquitous. A common question in many studies is whether there is an association between a response variable measured on each sample and the latent group structure represented by the tree. *treeSeg* is a statistical method with statistical guarantees that tests for association between the response variable and the tree structure across all levels of the tree hierarchy with high power, while accounting for the overall false positive error rate. The method is based on multiscale change point approach (Frick et. al 2014) applied to the tree structures. The *treeSeg* algorithm is implemented as an R package.

For more details about how *treeSeg* works please see the manuscript 

*Testing for dependence on tree structures*, Behr M, Ansari M, Munk A, Holmes C. (2020) PNAS, doi:10.1073/pnas.1912957117.
([Link to paper](https://www.pnas.org/content/117/18/9787))


## Using treeSeg as an R package

## Installation
The package can be installed in R using the commands:
```{r}
install.packages('devtools')
library(devtools)
devtools::install_github("merlebehr/treeSeg", subdir="TreeSeg")
```

You should then be able to load the package with:
```{r}
library(treeSeg)
```
Please see our R code in the testData directory to see an example of how to use the method and produce figures for an example transcriptomic micro-array dataset.

An example on how to use *treeSeg* on simulated data can be found in the jupyter notebook *illustrationTreeSeg.ipynb*.


# Help
Do contact us if you need any help with using the software or when you have suggestion on how to improve the implementation.

Email address: behr (at) berkeley (dot) edu
