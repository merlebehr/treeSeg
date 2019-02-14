#include "treeseg.h"


// [[Rcpp::export]]
NumericVector offspringDiffTip(int Node, NumericVector node, List tree){
  //offspring difference tips
  
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
    }
    if(flag){
      sol.push_back(currentCol[ancM.nrow()-1]);
    }
    flag=false;
  }
  return sol;
}