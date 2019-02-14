#include "treeseg.h"

// [[Rcpp::export]]
List getNewCandidates(int subroot, int aN, List tree, IntegerVector minIVar, LogicalMatrix isOffMat, 
                      bool without){
  int i,j;
  IntegerVector news= offspringAll(subroot, tree);
  news.push_front(subroot);
  arma::umat numbers(news.length(),news.length()+2);
  IntegerVector target(news.length()+1);
  target[0]=aN;
  
  //printf("get new candidate with number of active nodes equal %d \n", aN);
  //int helpNewsLength = news.length();
  //printf("number of all active nodes in this subtree is %d \n", helpNewsLength);
  
  int maxSize = 500000;
  
  for(i=0; i<news.length(); i++){
    numbers(i,0)=minIVar[news[i]]+1;
    numbers(i, news.length()+1)= i+1;
    target[i+1]=1;
    for(j=1; j<news.length()+1; j++){
      numbers(i,j)=isOffMat(news[i],news[j-1]);
    }
  }
  //printf("now calling subsetSum New \n");
  
  //List ans = subsetSum(numbers, target, mat(0,0), List::create(), maxSize);
  List ans = subsetSum(numbers, target, maxSize);
  
  //printf("now finished subsetSum New \n");
  
  
  List final(0);
  
  // if(ans.size() >= maxSize){
  //   printf("problem size too large \n");
  //   return(final);
  // }
  
  //int help = ans.length();
  //printf("Finished with %d candidates \n", help);
  
  for(i=0; i<ans.length(); i++){
    mat partialAns=ans[i];
    IntegerVector partialFinal(partialAns.n_rows);
    bool within=0;
    
    for(j=0; (unsigned)j<partialAns.n_rows; j++){
      partialFinal[j]=news[partialAns(j, partialAns.n_cols-1)-1];
      if(without&&(partialFinal[j]==subroot)){
        within=1;
        break;
      }
    }
    if(!within){
      final.push_back(partialFinal);
    }
  }
  
  return final;
}