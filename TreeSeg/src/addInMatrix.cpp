#include "treeseg.h"

// [[Rcpp::export]]
NumericMatrix addInMatrix(NumericVector vec, NumericMatrix M){
  int i;
  NumericMatrix Mnew(M.nrow()+1,M.ncol());
  
  Mnew(0,_)=vec;
  for(i=1; i<M.nrow()+1; i++){
    Mnew(i,_)=M(i-1,_);
  }
  
  return Mnew;
}