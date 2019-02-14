#include "treeseg.h"

// [[Rcpp::export]]
arma::mat allComb(IntegerVector lengths){
  arma::mat M(lengths[0], lengths.length(), fill::ones);
  
  for(int i=0; i<lengths[0]; i++){
    M(i,0)=i+1;
  }
  
  for(int i=1; i<lengths.length(); i++){
    arma::mat auxM(M.n_rows, M.n_cols);
    auxM=M;
    for(int j=1; j<lengths[i]; j++){
      for(int k=0; (unsigned)k<(auxM.n_rows); k++){
        auxM(k,i)=j+1;
      }
      M=join_cols(M,auxM);
    }
  }
  
  return M;
}
