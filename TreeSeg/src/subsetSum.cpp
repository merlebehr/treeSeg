#include "treeseg.h"

// [[Rcpp::export]]
List subsetSum (arma::umat numbers, IntegerVector target, int maxSize) {
  //return all subset of rows from numbers which sum up to target

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
  
  int pSize = partial.size();
  int anSize = ans.size();

  while(partial.size() > 0){
    
    uvec parNodes = partial(0);
    
    arma::umat par = numbers.rows(parNodes);
    partial.erase(partial.begin());
    
    
    if( pSize < partial.size() ){
      pSize = partial.size();
      anSize = ans.size();
    }
    
    
    if(par.n_rows > 1){
      
      S = sum(par,0);
      

      alls = !(S[0]==target[0]);
      for(i=0; i<parNodes.size(); i++){
        alls= alls + !(S[parNodes[i]+1]==target[parNodes[i]+1]);
        if(S[parNodes[i]+1]>target[parNodes[i]+1]){
          continue;
        }
      }
      
      //check whether non-related nodes would be violated
      if( (S[0] == target[0]) && (alls == 0) ){
        for(i = 0; i < (S.size() - 2); i++){
          if((S[i + 1] == 0) && ( numbers(i,0) > 1)){
            //printf("Violation neighboring condition \n");
            alls = 1;
            continue;
          }
        }
      }
      
      
      if(alls==0){
        ans.push_back(par);
        
        anSize = ans.size();
        //printf("Max ans List has size %d \n", anSize);
        
        
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
