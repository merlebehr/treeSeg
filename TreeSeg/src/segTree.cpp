#include "treeseg.h"

// [[Rcpp::export]]

List segTree(NumericVector y, List& tree, double q = NA_REAL, double alpha = 0.05, int fam = 0){
  
  int i, j, k, r, ri, li, indN, indI;
  
  if(!tree.containsElementNamed("ancM")){
    tree = prepTree(tree);
  }

  int n = y.length();
  
  //printf("family is %d \n", fam);
  

  List bou = boundsCall(y, alpha, q, fam);

  
  NumericVector lower = bou["lower"];                  
  NumericVector upper = bou["upper"];
  
  // double minUp = min(upper);
  // double maxLow = max(lower);
  // printf("minupper is %f \n", minUp);
  // printf("maxlower is %f \n", maxLow);
  

  if(min(upper) >= max(lower)){
    bool noAN = 1;
    double lowerMax = max(lower);
    double upperMin = min(upper);
    return List::create(Named("noAN")= noAN, Named("lower")=lowerMax, Named("upper")=upperMin);
  }
    
  
  int Nnode= tree["Nnode"];
  int Nn=Nnode+n;                                       //total number of nodes (internal + leaves)
  
  IntegerVector minI(Nn+1, IntegerVector::get_na());    //minimal number of important nodes 
  List comb(Nn+1);                                      //valid suboptimal solutions
  NumericVector optCost(Nn+1, NumericVector::get_na()); //ML cost of suboptimal combinations
  IntegerVector empty(0);
  
  LogicalMatrix isOffMat(Nn+1,Nn+1);
  IntegerMatrix ancM=tree["ancM"];
  
  for(i=0; i<Nn; i++){
    isOffMat(i+1,i+1)=1;
  }
  
  for(i=0; i<ancM.ncol(); i++){
    for(j=0; j<ancM.nrow()-1; j++){
      for(k=j+1; k<ancM.nrow(); k++){
        isOffMat(ancM(k,i), ancM(j,i))=1;
      }
    }
  }
  
  //Initialize solution for leave nodes (1,.., n)
  for(i=1;i<n+1;i++){
    minI[i]=0;
    optCost[i]=0;
    comb[i]= List::create(List::create(Named("comb")=empty, Named("minB")=upper[bouPos(n,i,i)], 
                                       Named("maxB")=lower[bouPos(n,i,i)]));
  }
  
  //Initialize ancestor vector
  IntegerVector anc(n);
  for(i=0;i<n;i++){
    anc[i]=i+1;
  }
  
  int calcSol=n;                                        //calculated solutions
  
  while(calcSol<Nn){
    //for(r=0; r<10; ++r){
    //printf("solution computed for %d out of %d inner nodes. \n", calcSol, Nn);
    
    IntegerVector aanc(anc.length());       //new ancestors
    IntegerVector aancMiss(Nn+1);           //aux vector to check if all the off of a node appear in a level
    LogicalVector aancFlags(Nn+1);          //aux logical vector
    LogicalVector ancRemove(anc.length());  //remove the repeated nodes in anc
    IntegerVector Anc(0);
    
    //calculate new ancestors (uniquely) + merge subtrees
    for(i=0;i<anc.length();i++){
      int auxAnc = ancestor(anc[i],tree);
      aanc[i]=auxAnc;
      
      if(aancFlags[auxAnc]==0){
        IntegerVector off=offspring(auxAnc, tree);
        aancMiss[auxAnc]=off.length();
        aancFlags[auxAnc]=1;
      }
      
      aancMiss[auxAnc]=aancMiss[auxAnc]-1;
      if(aancMiss[auxAnc]==0){
        Anc.push_back(auxAnc);
        ancRemove[i]=1;
        aancFlags[auxAnc]=0;
      }
    }
    
    //update current ancester vector 
    for(i=0; i<anc.length(); i++){
      if(ancRemove[i]==1){
        anc[i]=aanc[i];
      }
      else if(aancFlags[aanc[i]]==1){
        ancRemove[i]=1;
      }
    }
    anc=anc[ancRemove];
    
    //int lengthAnc = Anc.length();
    //printf("next compute %d solutions.\n", lengthAnc);
    
    for(i=0;i<Anc.length();i++){
      
      //suboptimal solutions for offsprings of Anc[i]
      IntegerVector off=offspring(Anc[i],tree);
      List combI(off.length());
      IntegerVector minII(off.length());
      List oftI(off.length());
      IntegerVector oftb(2*off.length());
      
      for(j=0; j<off.length(); j++){
        IntegerVector offTip=getOffspringTipB(off[j], tree);
        combI[j]=comb[off[j]];
        minII[j]=minI[off[j]];
        oftI[j]=offTip;
        oftb[2*j]=offTip[0];
        oftb[2*j+1]=offTip[1];
      }
      
      int maxoftb = max(oftb);
      int minoftb = min(oftb);
      
      if(is_true(all(minII==0))){
        double maxB = min(lower);
        double minB = max(upper);
        
        for(j=0; j < off.length(); j++){
          List auxcomb=combI[j];
          List onlyauxcomb=auxcomb[0];
          double maxaux=onlyauxcomb["maxB"];
          double minaux=onlyauxcomb["minB"];
          maxB=std::max(maxB, maxaux);
          minB=std::min(minB, minaux);
          
          if(j < off.length()-1){
            for(li = oftb[2*j]; li <= oftb[2*j+1]; li++){
              for(ri = oftb[2*j+2]; ri <= maxoftb; ri++){
                if(li <= ri){
                  int ind=bouPos(n,li,ri);
                  maxB=std::max(maxB, lower[ind]);
                  minB=std::min(minB, upper[ind]);
                }
              }
            }
          }
          
        }
        
        if(minB>=maxB){
          minI[Anc[i]]=0;
          IntegerVector ncomb(0);
          comb[Anc[i]]=List::create(List::create(Named("comb")=ncomb, Named("minB")=minB, Named("maxB")=maxB));
          LogicalVector range(n);
          for(j=minoftb-1; j<maxoftb; j++){
            range[j]=1;
          }
          optCost[Anc[i]]=mlCost(y[range], fam);
        }
        else{
          int minItest=0;
          
          while(minI[Anc[i]]==NA){
            minItest = minItest+1;
            List ncomb(0);
            NumericVector cost(0);
            
            IntegerVector minIVar=clone(minI);
            minIVar[Anc[i]]=0;
            
            List combNew = getNewCandidates(Anc[i], minItest, tree, minIVar, isOffMat,1); 
            
            if(combNew.length() == 0){
              //printf("Problem size too large \n");
              return(combNew);
            }
            
            //return(List::create(combNew));
            
            //int lengthHelp = combNew.length();
            //printf("next check combNew of length %d .\n", lengthHelp);
            
            for(j=0; j<combNew.length(); j++){
              IntegerVector auxComb = combNew[j];
              IntegerMatrix ofItb(auxComb.length(),2); //boundary of non-important offspring area of Anc[i]
              
              double maxB = min(lower);
              double minB = max(upper);
              
              // int lengthHelp = auxComb.length();
              // printf("next check auxComb of length %d .\n", lengthHelp);
              
              for(k=0; k<auxComb.length(); k++){
                IntegerVector auxB = getOffspringTipB(auxComb[k], tree);
                ofItb(k,0)=auxB[0];
                ofItb(k,1)=auxB[1];
                
             
                
                for(li=minoftb;li <= lower.length();li++){
                  for(ri = 0; ri <= maxoftb; ri++){
                    if((li<=ri) && !(li <= auxB[1] && ri >= auxB[0])){
                      int ind = bouPos(n,li,ri);
                      
                      // if(auxComb[k] == 151){
                      //   printf("we check li %d \n", li);
                      //   printf("we check ri %d \n", ri);
                      // }
                      
                      maxB = std::max(maxB, lower[ind]);
                      minB = std::min(minB, upper[ind]);
                    }
                  }
                }
              }
              
              if(minB>=maxB){
                auxComb.sort();
                ncomb.push_back(List::create(Named("comb")=auxComb, Named("minB")=minB, Named("maxB")=maxB));
                
                LogicalVector range(n);
                
                ofItb.sort();
                
                for(k=minoftb-1;k<ofItb[0]-1;k++){
                  range[k]=1;
                }
                
                for(r=0; r<(ofItb.length()-2)/2;r++){
                  for(k=ofItb[2*r+1]; k<ofItb[2*r+2]-1; k++){
                    range[k]=1;
                  }
                }
                
                for(k=ofItb[ofItb.length()-1]; k<maxoftb; k++){
                  range[k]=1;
                }
                double auxCost = mlCost(y[range], fam);
                for(k=0; k<auxComb.length(); k++){
                  auxCost = auxCost + optCost[auxComb[k]];
                }
                
                cost.push_back(auxCost);
                
              }
            }
            if(ncomb.length()>0){
              minI[Anc[i]] = minItest;
              sortComb(ncomb,cost);
              comb[Anc[i]]=ncomb;
              optCost[Anc[i]]=cost[0];
              
            }
            
            
          }
        }
      }
      else{
        List ncomb(0);
        NumericVector cost(0);
        //printf("check feasiblility of all combinations of suboptimal solutions from off\n");
        
        IntegerVector auxLengths(combI.length());
        
        
        // int lengthHelp = combI.length();
        // printf("next check comIlength of length %d .\n", lengthHelp);
        
        for(j=0; j<combI.length(); j++){
          List auxcomb=combI[j];
          auxLengths[j]=auxcomb.length();
        }
        
        arma:: mat newsI = allComb(auxLengths);
        
        // lengthHelp = newsI.n_rows;
        // printf("next check newsI.n_rows of length %d .\n", lengthHelp);
        
        for(indN=0; (unsigned)indN<newsI.n_rows; indN++){
          IntegerVector news(0);
          IntegerVector ofItb(0);
          double maxB = min(lower);
          double minB = max(upper);
          

          for(j=0; (unsigned)j<newsI.n_cols; j++){
            List auxcombs=combI[j];
            List auxcomb = auxcombs[newsI(indN,j)-1];
            List onlyAuxComb = auxcomb[0];
            
            for(k=0; k<onlyAuxComb.length(); k++){
              IntegerVector auxBounds = getOffspringTipB(onlyAuxComb[k], tree);
              
              news.push_back(onlyAuxComb[k]);
              ofItb.push_back(auxBounds[0]);
              ofItb.push_back(auxBounds[1]);
            }
            

            maxB = std::max(maxB, as<double>(auxcomb[2]));
            minB = std::min(minB, as<double>(auxcomb[1]));
          }
          ofItb.sort();
          
          for(j=0; j<off.length()-1; j++){
            IntegerVector subaux1 = ofItb[ofItb <= oftb[2*j+1]];
            IntegerVector subaux2 = ofItb[ofItb >= oftb[2*j+2]];
            
            
            
            
            int liStart;
            if(subaux1.size() == 0){
              liStart = oftb[2*j];
            }else{
              liStart = std::max(max(subaux1)+1, oftb[2*j]);
            }
            
            int riEnd;
            if(subaux2.size() == 0){
              riEnd = oftb[2*j+3];
            }else{
              riEnd = std::min(min(subaux2) - 1, oftb[2*j+3]);
            }
            
        
            
            for(li = liStart; li <= oftb[2*j+1]; li++){
              for(ri = oftb[2*j+1]+1; ri <= riEnd; ri++){
                if(li <= ri){
                  int ind=bouPos(n,li,ri);
                  
                  maxB=std::max(maxB, lower[ind]);
                  minB=std::min(minB, upper[ind]);
                }
              }
            }
          }
          
          if(minB >= maxB){
            news.sort();
            ncomb.push_back(List::create(Named("comb")=news, Named("minB")=minB, Named("maxB")=maxB));
            LogicalVector range(n);
            double auxCost=0;
            
            for(j=minoftb-1; j<maxoftb; j++){
              range[j]=1;
            }
            
            for(j=0; j<news.length(); j++){
              IntegerVector auxOff = getOffspringTip(news[j],tree);
              for(k=0; k<auxOff.length(); k++){
                range[auxOff[k]-1]=0;
              }
              auxCost = auxCost + optCost[news[j]];
            }
            
            auxCost = auxCost+mlCost(y[range], fam);
            cost.push_back(auxCost);
            
          }
        }
        if(ncomb.length()>0){
          int auxMinI = 0;
          
          for(j=0; j<off.length(); j++){
            auxMinI = auxMinI+minI[off[j]]; 
          }
          minI[Anc[i]]=auxMinI;
          sortComb(ncomb,cost);
          comb[Anc[i]]=ncomb;
          optCost[Anc[i]]=cost[0];
        }
        else{
          //printf("check solutions with additional change point \n");
          
          int minITest = 0;
          int addNodes = 0;
          
          for(j=0; j<off.length(); j++){
            minITest = minITest+minI[off[j]]; 
          }
          while(minI[Anc[i]]==NA){
            
            addNodes = addNodes + 1;
            minITest=minITest + addNodes;
            
            umat indAddNodes = mycombn(off.length(), addNodes);
            
            //lengthHelp = indAddNodes.n_cols;
            //printf("next check indAddNodes.n_cols of length %d .\n", lengthHelp);
            
            for(indI=0; (unsigned)indI<indAddNodes.n_cols; indI++){
              int auxInd = 0;
              List newsComb(off.length());
              IntegerVector auxLengths(off.length());
              for(j=0; j<off.length(); j++){
                
                if(((unsigned)auxInd<indAddNodes.n_rows)&&(((unsigned)j+1)==indAddNodes(auxInd, indI))){
                  
                  List auxNewsComb=getNewCandidates(off[j], minI[off[j]]+1, tree, minI, isOffMat);
                  
                  if(auxNewsComb.length() == 0){
                    //printf("Problem size too large! \n");
                    return(auxNewsComb);
                  }
                  
                  
                  newsComb[j]=auxNewsComb;
                  auxInd = auxInd+1;
                  auxLengths[j]=auxNewsComb.length();
                }
                else{
                  List auxComb=combI[j];
                  List auxNewsComb(auxComb.length());
                  
                  for(k=0; k<auxComb.length(); k++){
                    List auxComblist = auxComb[k];
                    auxNewsComb[k]=auxComblist[0];
                  }
                  
                  auxLengths[j]=auxComb.length();
                  newsComb[j]=auxNewsComb;
                }
              }
              
              mat newsI = allComb(auxLengths);
              
              //printf("newsI has size %d .\n", newsI.n_rows);
              
              for(indN=0; (unsigned)indN<newsI.n_rows; indN++){
                IntegerVector news(0);
                double maxB = min(lower);
                double minB = max(upper);
                LogicalVector rm(lower.length());
                
                for(j=0; (unsigned)j<newsI.n_cols; j++){
                  List auxcombs=newsComb[j];
                  IntegerVector auxcomb = auxcombs[newsI(indN,j)-1];
                  
                  for(k=0; k<auxcomb.length(); k++){
                    news.push_back(auxcomb[k]);
                  }
                  
                }
                
                
                for(j=0; j<news.length(); j++){
                  IntegerVector auxB = getOffspringTipB(news[j], tree);
                  
                  for(li = 1;li<=auxB[1];li++){
                    for(ri=auxB[0]; ri<=n; ri++){
                      if(li<=ri){
                        int ind=bouPos(n,li,ri);
                        rm[ind]=1;
                      }
                    }
                  }
                  
                }
                
                if(oftb[0]>1){
                  for(li = 1;li<=oftb[0]-1;li++){
                    for(ri=1; ri<=n; ri++){
                      if(li<=ri){
                        int ind=bouPos(n,li,ri);
                        rm[ind]=1;
                      }
                    }
                  }
                }
                
                if(oftb[3]<n){
                  for(li = 1;li<=n;li++){
                    for(ri=oftb[4]+1; ri<=n; ri++){
                      if(li<=ri){
                        int ind=bouPos(n,li,ri);
                        rm[ind]=1;
                      }
                    }
                  }
                }
                
                minB=vecmin(upper[!rm]);
                maxB=vecmax(lower[!rm]);
                
                
                if(minB >= maxB){
                  news.sort();
                  ncomb.push_back(List::create(Named("comb")=news, Named("minB")=minB, Named("maxB")=maxB));
                  LogicalVector range(n);
                  double auxCost=0;
                  
                  for(j=minoftb-1; j<maxoftb; j++){
                    range[j]=1;
                  }
                  
                  for(j=0; j<news.length(); j++){
                    IntegerVector auxOff = getOffspringTip(news[j],tree);
                    for(k=0; k<auxOff.length(); k++){
                      range[auxOff[k]-1]=0;
                    }
                    auxCost = auxCost + optCost[news[j]];
                  }
                  
                  auxCost = auxCost+mlCost(y[range], fam);
                  cost.push_back(auxCost);
                }
              }
              
            }
            
            if(ncomb.length()>0){
              minI[Anc[i]]= minITest;
              sortComb(ncomb,cost);
              comb[Anc[i]]=ncomb;
              optCost[Anc[i]]=cost[0];
            }
            
          }
          
        }
        
      }
    }
    
    
    calcSol=calcSol+Anc.length();
  }
  
  return List::create(Named("minI")=minI,Named("comb")=comb, Named("optCost")=optCost, Named("root")=tree["root"]);
}

