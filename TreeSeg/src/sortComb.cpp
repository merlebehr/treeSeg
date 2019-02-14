#include "treeseg.h"

// [[Rcpp::plugins("cpp11")]]

// [[Rcpp::export]]
List sortComb(List& ncomb, NumericVector& cost){
  // sorts the combinations [ncomb] according to [cost]
  IntegerVector idx = seq_along(cost) - 1;
  // Then sort that vector by the values of y
  std::sort(idx.begin(), idx.end(),[&](int i, int j){return cost[i] < cost[j];});
  // And return x in that order
  ncomb = ncomb[idx];
  cost = cost[idx];
  return List::create(ncomb,cost);
}