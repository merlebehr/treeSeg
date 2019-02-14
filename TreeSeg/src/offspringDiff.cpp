#include "treeseg.h"

// [[Rcpp::export]]
NumericVector offspringDiff(int Node, NumericVector node, List tree){
  //offspring difference function
  
  int i,j;
  NumericMatrix ancM=tree["ancM"];
  NumericVector currentCol(ancM.nrow());
  NumericVector sol(0);
  bool flag=false;
  
  //look in each column of ancM
  for(i=0;i<ancM.ncol();i++){
    currentCol=ancM(_,i);
    for(j=0;j<ancM.nrow();j++){
      if(currentCol[j]==Node){
        flag=true;
      }
      if(flag&&(std::find(node.begin(),node.end(),currentCol[j])!=node.end())){
        flag=false;
      }
      if(flag&&currentCol[j]!=Node&&(std::find(sol.begin(),sol.end(),currentCol[j])==sol.end())){
        sol.push_back(currentCol[j]);
      }
    }
    flag=false;
  }
  return sol;
}
