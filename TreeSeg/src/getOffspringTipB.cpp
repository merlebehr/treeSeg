#include "treeseg.h"

// [[Rcpp::export]]
IntegerVector getOffspringTipB(int node, List tree){
  // get offspring tip bounds
  
  IntegerVector tip=getOffspringTip(node,tree);
  tip.sort();
  
  return IntegerVector::create(tip[0],tip[tip.length()-1]);
}