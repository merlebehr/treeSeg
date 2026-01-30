#include "treeseg.h"

// [[Rcpp::export]]
arma::umat mycombn(double n, double k) {
  // similar to combn
  
  double n_subsets = Rf_choose(n, k);
  umat out = zeros<umat>(k, n_subsets);
  uvec a = linspace<uvec>(1, k, k);  
  out.col(0) = a;
  int m = 0;
  int h = k;
  uvec j;
  
  for(long long i = 1; i < n_subsets; i++){
    if(m < (n - h)){  
      h = 1;
      m = a(k - 1);
      j = linspace<uvec>(1, 1, 1);
    }
    else{
      m = a(k - h - 1);
      ++h;
      j = linspace<uvec>(1, h, h);
    }
    a.elem(k - h + (j-1) ) = m + j;
    out.col(i) = a;
  }
  return(out);
}