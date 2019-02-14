#include "treeseg.h"

// [[Rcpp::export]]
double mlCost(NumericVector y, int fam){
  //cost of y, without an ordered unique vector
  
  int n = y.length();
  double mean = 0, cost = 0;
  for(int i = 0; i < n; i++){
    mean+=y[i];
  }
  mean/=n;
  if(fam == 0){
    for(int i = 0; i < n; i++){
      cost+=-log(pow(mean,y[i]) * pow(1-mean,1-y[i]));
    }
  }
  if(fam == 1){
    for(int i = 0; i < n; i++){
      cost+= -pow(y[i] - mean,2);
    }
  }
  
  
  return cost;
}