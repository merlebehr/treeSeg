#include "treeseg.h"

// [[Rcpp::export]]
NumericVector colSums(const arma::mat & X){
  int nCols = X.n_cols;
  NumericVector out(nCols);
  for(int i = 0; i < nCols; i++){
    out(i) = sum(X.col(i));
  }
  return(out);
}