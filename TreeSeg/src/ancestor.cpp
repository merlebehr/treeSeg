#include "treeseg.h"

// [[Rcpp::export]]
int ancestor(int node, List tree){
  //ancestor function
  
  int i;
  
  int root=tree["root"];
  NumericMatrix edge=tree["edge"];
  
  if(node==root){
    return node;
  }
  else{
    for(i=0; i<edge.nrow();i++){
      if(edge(i,1)==node){
        return edge(i,0);
      }
    }
  }
  return -1;
}