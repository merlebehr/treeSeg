#include "treeseg.h"

// [[Rcpp::export]]
List getNewCandidates(int subroot, int aN, List tree, IntegerVector minIVar, LogicalMatrix isOffMat, 
                      bool without){
  // get all candidate solutions for given subroot with aN number of active nodes
  // minimal number of active nodes at other nodes is given in minIVar
  // logical matrix isOffMat indicates whether node i is offspring of node j
  
  int i,j;
  IntegerVector news= offspringAll(subroot, tree); // all offsprings nodes of subroot (candidates for new active nodes)
  news.push_front(subroot);
  arma::umat numbers(news.length(), news.length()+2); // numbers that have to add up to target in subsetsub
  IntegerVector target(news.length()+1); //target for subset sum
  target[0]=aN;
  
  
  int maxSize = 500000; // maximal number of cadidate solutions feasible in terms of memory
  
  for(i=0; i<news.length(); i++){
    numbers(i,0)=minIVar[news[i]]+1;
    numbers(i, news.length()+1)= i+1;
    target[i+1]=1;
    for(j=1; j<news.length()+1; j++){
      numbers(i,j)=isOffMat(news[i],news[j-1]);
    }
  }

  List ans = subsetSum(numbers, target, maxSize);
  
  
  List final(0); //feasible candidates
  
  
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