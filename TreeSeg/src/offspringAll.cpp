#include "treeseg.h"
//' All offspring nodes of given inner node
//' 
//' Gives back all nodes of a \code{tree} which are an offspring of 
//' a given \code{node} of the tree
//' @param node single integer value which contrains the index of a node from \code{tree}
//' @param tree rooted tree object
//' @return Integer vector with all indeces of offspring nodes of \code{node} in \code{tree}.
//' @examples 
//' library(ape)
//' n <- 2 * 50
//' tree <- rtree(n, rooted = TRUE)
//' aktiveNode <- n + 50
//' offAktiveNode <- offspringAll(aktiveNode,tree)
//' @export
// [[Rcpp::export]]
IntegerVector offspringAll(int node, List tree){
  // all the offspring
  
  if(!tree.containsElementNamed("ancM")){
    tree = prepTree(tree);
  }
  
  int i,j;
  NumericMatrix ancM=tree["ancM"];
  NumericVector currentCol(ancM.nrow());
  IntegerVector sol(0);
  bool flag=false;
  
  //look in each column of ancM
  for(i=0;i<ancM.ncol();i++){
    currentCol=ancM(_,i);
    for(j=0;j<ancM.nrow();j++){
      if(flag&&(std::find(sol.begin(),sol.end(),currentCol[j])==sol.end())){
        sol.push_back(currentCol[j]);
      }
      if(currentCol[j]==node){
        flag=true;
      }
    }
    flag=false;
  }
  std::sort(sol.begin(),sol.end());
  return sol;
}
