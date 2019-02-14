#include "treeseg.h"

//' Offspring tips of given inner node
//' 
//' Gives back all tips of a \code{tree} which are an offspring of 
//' a given \code{node} of the tree
//' @param node single integer value which contrains the index of a node from \code{tree}
//' @param tree rooted tree object
//' @return Integer vector with all indeces of offspring tips of \code{node} in \code{tree}.
//' @examples 
//' library(ape)
//' n <- 2 * 50
//' tree <- rtree(n, rooted = TRUE)
//' aktiveNode <- n + 50
//' offAktiveNode <- getOffspringTip(aktiveNode,tree)
//' @export
// [[Rcpp::export]]
IntegerVector getOffspringTip(int node, List tree){
  // get offspring tip
  
  if(!tree.containsElementNamed("ancM")){
    tree = prepTree(tree);
  }
  
  
  int i,j;
  IntegerMatrix ancM=tree["ancM"];
  IntegerVector currentCol(ancM.nrow());
  IntegerVector sol(0);
  
  //look in each column of ancM
  for(i=0;i<ancM.ncol();i++){
    currentCol=ancM(_,i);
    for(j=0;j<ancM.nrow();j++){
      if(currentCol[j]==node){
        sol.push_back(currentCol[ancM.nrow()-1]);
        break;
      }
    }
  }
  return sol;
}