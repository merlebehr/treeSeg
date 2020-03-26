#include "treeseg.h"

//' @useDynLib treeSeg
//' @importFrom Rcpp evalCpp

// [[Rcpp::export]]
List prepTree(List& tree){
  // tree preparation
  
  int i;
  int root=-1;
  
  StringVector tipLabel=tree["tip.label"];
  int Nnode=tree["Nnode"];
  NumericMatrix edge=tree["edge"];
  NumericVector childs=edge( _, 1);
  
  //add root
  for(i=Nnode+1;i<=Nnode+ tipLabel.size();i++){
    if(std::find(childs.begin(),childs.end(),i)==childs.end()){
      root=i;
    }
  }
  tree["root"]=root;
  
  //add ancestor matrix
  NumericMatrix ancM(1,tipLabel.size());
  NumericVector anc(tipLabel.size());
  bool flag=false;
  
  for(i=0;i<tipLabel.size();i++){
    ancM(0,i)=i+1;
  }
  
  while(!flag){
    flag=true;
    for(i=0; i<tipLabel.size(); i++){
      anc[i]= ancestor(ancM(0,i), tree);
      flag=flag&&(anc[i]==root);
    }
    ancM=addInMatrix(anc,ancM);
  }
  tree["ancM"]=ancM;
  
  return tree;
}
