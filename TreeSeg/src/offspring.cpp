#include "treeseg.h"

// [[Rcpp::export]]
IntegerVector offspring(int node, List tree){
  //offspring function
  
  int i;
  IntegerVector off(0);
  
  IntegerMatrix edge=tree["edge"];
  for(i=0; i<edge.nrow();i++){
    if(edge(i,0)==node){
      off.push_back(edge(i,1));
    }
  }
  if(off.size()==0){
    off=node;
  }
  return off;
}