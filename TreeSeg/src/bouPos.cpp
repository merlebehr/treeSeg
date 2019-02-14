#include "treeseg.h"


// [[Rcpp::export]]
double bouPos(int ny,int li,int ri){
  // extracts the position of li and ri from the output of bounds_call. [ny] is the length of the observations.
  int liStart = (li-1)*ny-((li-1)*(li-2))/2;
  return(liStart+ri-li);
}