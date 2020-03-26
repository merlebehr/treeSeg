#include "treeseg.h"


// [[Rcpp::export]]
double bouPos(int li,int ri, int allInt, int liStart){
  // extracts the position of li and ri from the output of bounds_call. [liStart] is the start position of li in the bounds object.
  int ans;
  if(allInt == 1){
    // all intervals
    ans = liStart + ri - li;
  }
  else{
    // dyadic intervals
    int logrl = (int) log2(ri - li + 1);
    
    if(pow(2, logrl) < ri - li + 1){
      ans = -1;
    }
    else{
      ans = liStart + logrl;
    }
  }
  
  //printf("li is %d and ri is %d and pos is %d \n", li, ri, ans);
  
  return(ans);
}