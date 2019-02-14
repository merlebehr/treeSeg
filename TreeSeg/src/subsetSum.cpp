#include "treeseg.h"

// [[Rcpp::export]]
List subsetSum (arma::umat numbers, IntegerVector target, int maxSize) {

  List partial; //list of partial row indeces of numbers which is succesively updated
  List ans; //list with row index sets of numbers which sum up to target
  unsigned i;
  int ncols= numbers.n_cols;
  IntegerVector S(ncols); //sum of partial rows 
  int alls;
  List help(maxSize);
  

  //initialize list of partial indices
  for(i=0; i<numbers.n_rows; i++){
    uvec news(1);
    news(0) = i;
    partial.push_back(news);
  }
  
  // int pSize = partial.size();
  // printf("Initial length of partial is %d \n", pSize);
  
  while(partial.size() > 0){
    
    //int pSize = partial.size();
    //printf("\r Current partial list has size %d", pSize);
    //fflush(stdout);
    //printf(" Current partial List has size %d \n", pSize);
   
    
    uvec parNodes = partial(0);
    
    arma::umat par = numbers.rows(parNodes);
    partial.erase(partial.begin());
    
    
    //pSize = partial.size();
    //printf("After erase partial List has size %d \n", pSize);
    
    
    if(par.n_rows > 1){
      
      S = sum(par,0);
      

      alls = !(S[0]==target[0]);
      for(i=0; i<parNodes.size(); i++){
        alls= alls + !(S[parNodes[i]+1]==target[parNodes[i]+1]);
        if(S[parNodes[i]+1]>target[parNodes[i]+1]){
          continue;
        }
      }
      if(alls==0){
        ans.push_back(par);
        if(ans.size() >= maxSize){
          return(help);
        }
        continue;
      }
      else if(S[0]>=target[0]){
        continue;
      }
      
    //add new partials
    for(i = parNodes[parNodes.size()-1] + 1; i<numbers.n_rows; i++){
      uvec news = parNodes;
      news.resize(news.n_rows+1);
      news(news.n_rows-1) = i;
      
      // printf("\n We are adding news:");
      // for(j = 0; j < news.n_rows; j++){
      //   printf("add %d ", news(j));
      // }
      // printf("Done with adding \n");
      
      
      partial.push_back(news);
    }
      
    }else{
      if(par[0]==target[0]){
        ans.push_back(par);
        if(ans.size() >= maxSize){
          return(help);
        }
        continue;
      }
      else if(par[0] > target[0]){
        continue;
      }
      
      //add new partials
      for(i = parNodes(0) + 1; i<numbers.n_rows; i++){
        uvec news = parNodes;
        news.resize(news.n_rows + 1);
        news(news.n_rows - 1) = i; 
        
        // printf("\n We are adding news:");
        // for(j = 0; j < news.n_rows; j++){
        //   printf("add %d ", news(j));
        // }
        // printf("Done with adding \n");
        
        partial.push_back(news);
      }
    }
  }
  return ans;
}
