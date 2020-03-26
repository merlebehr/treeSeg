#include "treeseg.h"

// [[Rcpp::export]]


List segTree(NumericVector y, NumericVector lengths, List& tree, double q = NA_REAL, double alpha = 0.05, int fam = 0){
  
  int i, j, k, r, li, indN, indI, allInt, iRiS, iRiE, riS, riE; //ri
  IntegerVector startLi(0);

  int n = y.length();
  
  if(n == lengths.length()){
    // interval system with all lengths
    allInt = 1;
    for(int li = 1; li <= n; li++){
      startLi.push_back( (li-1)*n - ((li-1)*(li-2))/2 );
    }
  }
  else{
    // interval system with dyadic lengths
    allInt = 0;
    for(int li = 1; li <= n; li++){
      int listart = li - 1;
      for(int i = 0; i < li - 1; i++){
        listart = listart + (int) log2(n - i);
      }
      startLi.push_back(listart);
    }
  }
  
  // generate multiscale bounds from stepR package
  List bou = boundsCall(y, lengths, alpha, q, fam);
  
  NumericVector lower = bou["lower"];                  
  NumericVector upper = bou["upper"];
  

  if(min(upper) >= max(lower)){
    // multiscale constraints can be satisfied without active nodes, return constant solution
    bool noAN = 1;
    double lowerMax = max(lower);
    double upperMin = min(upper);
    return List::create(Named("noAN")= noAN, Named("lower")=lowerMax, Named("upper")=upperMin);
  }
  
  if(!tree.containsElementNamed("ancM")){
    // check whether tree contrains ancestor matrix
    tree = prepTree(tree);
  }
    
  //start running treeSeg
  int Nnode= tree["Nnode"];
  int Nn=Nnode+n;                                       //total number of nodes (internal + leaves)
  
  IntegerVector minI(Nn+1, IntegerVector::get_na());    //minimal number of important nodes for sub-problems
  List comb(Nn+1);                                      //valid solutions for sub-problems
  NumericVector optCost(Nn+1, NumericVector::get_na()); //ML cost of sub-problems
  IntegerVector empty(0);
  
  LogicalMatrix isOffMat(Nn+1,Nn+1);                    //logical matrix which indicating whether one node is offspring of another
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
  for(i=1; i <= n; i++){
    minI[i]=0;
    optCost[i]=0;
    int ind = startLi[i - 1];
    //bouPos(i,i,allInt,startLi[i-1])
    comb[i]= List::create(List::create(Named("comb")=empty, Named("minB")=upper[ind], Named("maxB")=lower[ind]));
  }
  
  //Initialize ancestor vector
  IntegerVector anc(n);
  for(i=0;i<n;i++){
    anc[i]=i+1;
  }
  
  int calcSol=n;                                        //calculated solutions
  
  while(calcSol<Nn){
    //printf("solution computed for %d out of %d inner nodes. \n", calcSol, Nn);
    
    IntegerVector aanc(anc.length());       //new ancestors
    IntegerVector aancMiss(Nn+1);           //aux vector to check if all the off of a node appear in a level
    LogicalVector aancFlags(Nn+1);          //aux logical vector
    LogicalVector ancRemove(anc.length());  //remove the repeated nodes in anc
    IntegerVector Anc(0);                   //nodes for which optimal solution is computed next
    
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
    
    
    for(i=0;i<Anc.length();i++){
      
      //printf("Currently consider AN %d \n", Anc[i]);
      
      //compute optimal solutions for root node Anc[i]
      IntegerVector off=offspring(Anc[i],tree);         //direct offsprings of Anc[i]
      List combI(off.length());                         //valid solutions for those offsprings
      IntegerVector minII(off.length());                //minimal active nodes of offsprings
      List oftI(off.length());                          //offspring tip intervals
      IntegerVector oftb(2*off.length());               //individual offspring tip bounds
      
      for(j=0; j<off.length(); j++){
        IntegerVector offTip=getOffspringTipB(off[j], tree);
        combI[j]=comb[off[j]];
        minII[j]=minI[off[j]];
        

        oftI[j]=offTip;
        oftb[2*j]=offTip[0];
        oftb[2*j+1]=offTip[1];
      }
      
      int maxoftb = max(oftb);                          //most right offspring of Anc[i]
      int minoftb = min(oftb);                          //most left offspring of Anc[i]
      

      if(is_true(all(minII==0))){
        //both offsprings have no active nodes

        double maxB = min(lower);         //initialize maximum of lower bounds
        double minB = max(upper);         //initialize minimum of upper bounds
        
        for(j=0; j < off.length(); j++){
          List auxcomb=combI[j];
          List onlyauxcomb=auxcomb[0];
          double maxaux=onlyauxcomb["maxB"];
          double minaux=onlyauxcomb["minB"];
          maxB=std::max(maxB, maxaux);
          minB=std::min(minB, minaux);
          
          if(j < off.length()-1){
            //intersect with overlapping bounds from two offsprings
            
            for(li = oftb[2*j]; li <= oftb[2*j+1]; li++){
              
              riS = std::max(oftb[2*j+2], li);        //start value for right bound
              riE = maxoftb;                          //end value for right bound
              
              IntegerVector iRi = getiRi(riS, riE, li, allInt, lengths); //indeces in lengths vector for right and left bound
              
              iRiS = iRi[0];
              iRiE = iRi[1];
              
              for(int i1 = iRiS; i1 <= iRiE; i1++){
                //only loop over bounds of length in lengths
                
                //ri = li - 1 + lengths[i1];
                //int ind = bouPos(li, ri, allInt, startLi[li-1]);
      
                int ind = startLi[li - 1] + i1;
                
                if(ind >= 0){ //when no bounds exists for this interval ind = -1
                  maxB=std::max(maxB, lower[ind]);
                  minB=std::min(minB, upper[ind]);
                }
                
              }
            }
          }
          
        }
        
        if(minB>=maxB){
          //multiscale constraint satisfied, no active node needed
          
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
          //multscale constaint not satisfied, need to add active node(s) (at most one active node for binary trees)

          int minItest=0;
          
          while(minI[Anc[i]] == NA){
            minItest = minItest+1;          //current number of active nodes
            

            List ncomb(0);
            NumericVector cost(0);
            
            IntegerVector minIVar=clone(minI);
            minIVar[Anc[i]]=0;
            
            List combNew = getNewCandidates(Anc[i], minItest, tree, minIVar, isOffMat,1); //get all possible candiates for new active nodes
            
            if(combNew.length() == 0){
              //Problem size too large (too many possible candidates) 
              // to avoid memory issue stop calculations. 
              // Change maxSize in getNewCandidates to increase feasible problem size.
              return(combNew);
            }
            
            
            for(j=0; j<combNew.length(); j++){
              //check for all new candiates whether multiscale constraint can be fulfilled
              
              IntegerVector auxComb = combNew[j];
              IntegerMatrix ofItb(auxComb.length(),2); //boundary of non-important offspring area of Anc[i]
              
              double maxB = min(lower);
              double minB = max(upper);
              
              
              for(k=0; k<auxComb.length(); k++){
                //intersect bounds of non-important offsping area of Anc[i] (which is not in influence region of potential active node)
                IntegerVector auxB = getOffspringTipB(auxComb[k], tree);

                ofItb(k,0)=auxB[0];
                ofItb(k,1)=auxB[1];
                
                for(li = minoftb; li <= maxoftb; li++){
                  
                  int riMax = maxoftb;
                  if(li <= auxB[1]){
                    riMax = std::min(riMax, auxB[0] - 1);
                  }
                  
                  riS = li;
                  riE = riMax;
                  
                  IntegerVector iRi = getiRi(riS, riE, li, allInt, lengths);
                  
                  iRiS = iRi[0];
                  iRiE = iRi[1];
                  
                  
                  for(int i1 = iRiS; i1 <= iRiE; i1++){

                    //ri = li - 1 + lengths[i1];
                    //int ind = bouPos(li, ri, allInt, startLi[li-1]);
                    
                    int ind = startLi[li - 1] + i1;
         
                    if(ind >= 0){ //when no bounds exists for this interval ind = -1
                      maxB = std::max(maxB, lower[ind]);
                      minB = std::min(minB, upper[ind]);

                    }
                  }
                }
              }
              
              if(minB >= maxB){
                //candidiate gives valid solution, add this to list of all valid solutions
                
                auxComb.sort();
                ncomb.push_back(List::create(Named("comb")=auxComb, Named("minB")=minB, Named("maxB")=maxB));
                
                LogicalVector range(n);
                
                ofItb.sort();
                
                for(k = minoftb-1; k < ofItb[0]-1; k++){
                  range[k]=1;
                }
                
                for(r=0; r<(ofItb.length()-2)/2;r++){
                  for(k = ofItb[2*r+1]; k < ofItb[2*r+2]-1; k++){
                    range[k]=1;
                  }
                }
                
                for(k = ofItb[ofItb.length()-1]; k < maxoftb; k++){
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
              //Valid solutions found.

              minI[Anc[i]] = minItest;
              sortComb(ncomb,cost);
              comb[Anc[i]]=ncomb;
              optCost[Anc[i]]=cost[0];
              
            }
          }
        }
      }
      else{
        //offsprings do have some active nodes
        //first, check feasibility of all combinations of valid solutions from offsprings
        
        List ncomb(0);
        NumericVector cost(0);
        
        IntegerVector auxLengths(combI.length());
        
        
        for(j=0; j<combI.length(); j++){
          List auxcomb = combI[j];  //valid solutions of offspring j
          auxLengths[j] = auxcomb.length(); //number of valid solutions of offspring j
        }
        
        arma:: mat newsI = allComb(auxLengths); //possible combinations of valid solutions of individual offsprings
        
        for(indN = 0; (unsigned)indN < newsI.n_rows; indN++){
          IntegerVector news(0);
          IntegerVector ofItb(0);
          double maxB = min(lower);
          double minB = max(upper);
          

          for(j = 0; (unsigned)j < newsI.n_cols; j++){
            List auxcombs = combI[j];
            List auxcomb = auxcombs[newsI(indN,j)-1];
            List onlyAuxComb = auxcomb[0]; //current combination of active nodes under consideration
            
            for(k=0; k<onlyAuxComb.length(); k++){
              IntegerVector auxBounds = getOffspringTipB(onlyAuxComb[k], tree);
              
              news.push_back(onlyAuxComb[k]);
              ofItb.push_back(auxBounds[0]);
              ofItb.push_back(auxBounds[1]);
            }
            

            maxB = std::max(maxB, as<double>(auxcomb[2])); //get multiscale lower bound for ANs of current combination
            minB = std::min(minB, as<double>(auxcomb[1])); //get multiscale upper bound for ANs of current combination
          }
          ofItb.sort();
          
          
          for(j=0; j < off.length()-1; j++){
            //get indeces of intervals of non-important region
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
              //intersect with bounds on non-important region (not in influence region of ANs)
              
              riS = std::max(oftb[2*j+1]+1, li);
              riE = riEnd;
              
              IntegerVector iRi = getiRi(riS, riE, li, allInt, lengths);
              
              iRiS = iRi[0];
              iRiE = iRi[1];
              
              for(int i1 = iRiS; i1 <= iRiE; i1++){
                
                //ri = li - 1 + lengths[i1];
                //int ind = bouPos(li, ri, allInt, startLi[li-1]);
                
                int ind = startLi[li - 1] + i1;
 
                if(ind >= 0){ //when no bounds exists for this interval ind = -1
                  maxB=std::max(maxB, lower[ind]);
                  minB=std::min(minB, upper[ind]);
                }
              }
            }
          }
          
          if(minB >= maxB){
            //current candidate is valid solution, add to list of valid solutions
            
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
        if(ncomb.length() > 0){
          //Valid solutions were found, no additional active node needs to be added

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
          //No valid solutions were found, additional active node needs to be added
          
          //printf("No valid solutions were found, additional active node.\n");
          
          int minITest = 0;
          int addNodes = 0;
          
          for(j = 0; j < off.length(); j++){
            minITest = minITest + minI[off[j]]; 
          }
          while(minI[Anc[i]] == NA){
            
            addNodes = addNodes + 1;
            minITest = minITest + addNodes;               //current number of active nodes to be tested

            umat indAddNodes = mycombn(off.length(), addNodes);   //all combinations at which offspring node active node gets added
            

            for(indI = 0; (unsigned)indI < indAddNodes.n_cols; indI++){
              //check solutions with AN added at indI
              

              int auxInd = 0; //number of active nodes added
              List newsComb(off.length()); //all combinations for different offsprings
              IntegerVector auxLengths(off.length());
              
              for(j = 0; j < off.length(); j++){
                
                if(((unsigned) auxInd < indAddNodes.n_rows) && (((unsigned) j+1) == indAddNodes(auxInd, indI))){ 
                  // if  offspring j is candidiate for which active node gets add get all new candidates
                  
                  List auxNewsComb=getNewCandidates(off[j], minI[off[j]]+1, tree, minI, isOffMat);
                  
                  if(auxNewsComb.length() == 0){
                    //Problem size too large (too many possible candidates) 
                    // to avoid memory issue stop calculations. 
                    // Change maxSize in getNewCandidates to increase feasible problem size.
                    return(auxNewsComb);
                  }
                  
                  newsComb[j] = auxNewsComb;
                  auxInd = auxInd + 1;
                  auxLengths[j] = auxNewsComb.length();
                }
                else{
                  // if offspring j is not the caditate for which active node gets added take all current valid solutions
                  List auxComb = combI[j];
                  List auxNewsComb(auxComb.length()); 
                  
                  for(k = 0; k < auxComb.length(); k++){
                    List auxComblist = auxComb[k];
                    auxNewsComb[k]=auxComblist[0]; //get active nodes for k-th candidate
                  }
                  
                  auxLengths[j]=auxComb.length();
                  newsComb[j]=auxNewsComb;
                }
              }
              
              mat newsI = allComb(auxLengths);
              
              for(indN = 0; (unsigned) indN < newsI.n_rows; indN++){
                //check feasibility of all new candidates
                
                IntegerVector news(0);
                double maxB = min(lower);
                double minB = max(upper);
                LogicalVector rm(lower.length());
                
                for(j=0; (unsigned)j < newsI.n_cols; j++){
                  
                  List auxcombs = newsComb[j];
                  IntegerVector auxcomb = auxcombs[newsI(indN,j)-1];
                  
                  for(k = 0; k < auxcomb.length(); k++){
                    news.push_back(auxcomb[k]);
                  }
                }
                
              
                for(j=0; j<news.length(); j++){
                  //for each AN remove bounds which intersect with their influence region

                  IntegerVector auxB = getOffspringTipB(news[j], tree);

                  for(li = 1;li <= auxB[1]; li++){
                    
                    riS = std::max(auxB[0], li);
                    riE = n;
                    
                    IntegerVector iRi = getiRi(riS, riE, li, allInt, lengths);
                    
                    iRiS = iRi[0];
                    iRiE = iRi[1];

                    
                    for(int i1 = iRiS; i1 <= iRiE; i1++){
                      //ri = li - 1 + lengths[i1];
                      //int ind = bouPos(li, ri, allInt, startLi[li-1]);
                      
                      int ind = startLi[li - 1] + i1;

                      if(ind >= 0){ //when no bounds exists for this interval ind = -1
                        rm[ind] = 1;
                      }
                    }
                  }
                }
                
                if(minoftb > 1){
                  //remove bounds left of influence region of AN[j]

                  for(li = 1;li <= minoftb - 1;li++){
                    
                    riS = li;
                    riE = n;
                    
                    IntegerVector iRi = getiRi(riS, riE, li, allInt, lengths);
                    
                    iRiS = iRi[0];
                    iRiE = iRi[1];
                    


                    for(int i1 = iRiS; i1 <= iRiE; i1++){
                      //ri = li - 1 + lengths[i1];
                      //int ind = bouPos(li, ri, allInt, startLi[li-1]);
                      
                      int ind = startLi[li - 1] + i1;

                      if(ind >= 0){ //when no bounds exists for this interval ind = -1
                        rm[ind]=1;
                      }
                      
                    }
                  }
                }
                
                if(maxoftb < n){
                  //remove bounds right of influence region of AN[j]
                  
                  
                  for(li = 1;li <= n; li++){
                    
                    riS = std::max(maxoftb + 1, li);
                    riE = n;
                    
                    IntegerVector iRi = getiRi(riS, riE, li, allInt, lengths);
                    
                    iRiS = iRi[0];
                    iRiE = iRi[1];
                    

                    for(int i1 = iRiS; i1 <= iRiE; i1++){
                      //ri = li - 1 + lengths[i1];
                      //int ind = bouPos(li, ri, allInt, startLi[li-1]);
                      
                      int ind = startLi[li - 1] + i1;
                      
                      if(ind >= 0){//when no bounds exists for this interval ind = -1
                        rm[ind] = 1;
                      }
                    }
                  }
                }
                
                minB=vecmin(upper[!rm]);
                maxB=vecmax(lower[!rm]);
                
                
                if(minB >= maxB){
                  //candidate is valid solution, add to list of valid solutions

                  news.sort();
                  ncomb.push_back(List::create(Named("comb")=news, Named("minB")=minB, Named("maxB")=maxB));
                  LogicalVector range(n);
                  double auxCost = 0;
                  
                  for(j = minoftb-1; j < maxoftb; j++){
                    range[j] = 1;
                  }
                  
                  for(j = 0; j < news.length(); j++){
                    IntegerVector auxOff = getOffspringTip(news[j],tree);
                    for(k=0; k < auxOff.length(); k++){
                      range[auxOff[k]-1] = 0;
                    }
                    auxCost = auxCost + optCost[news[j]];
                  }
                  
                  auxCost = auxCost + mlCost(y[range], fam);
                  cost.push_back(auxCost);
                }
              }
            }
            
            if(ncomb.length() > 0){
              //Valid solutions found, no need to add more active nodes

              minI[Anc[i]]= minITest;
              sortComb(ncomb,cost);
              comb[Anc[i]]=ncomb;
              optCost[Anc[i]]=cost[0];
            }
          }
        }
      }
    }
    
    calcSol = calcSol + Anc.length();         //update number of solved sub-problems
  }
  
  return List::create(Named("minI") = minI,Named("comb") = comb, Named("optCost") = optCost, Named("root") = tree["root"]);
}







IntegerVector getiRi(int riS, int riE, int li, int allInt, const NumericVector& lengths){
  //input: start and end point of right index ri, riS and riE; left index li; interval system allInt with length lengths
  //output: respective indeces of riS and riE in lengths + li -1
  
  int iRiS, iRiE;
  
  if(allInt == 1){
    iRiS = riS - li;
    iRiE = riE - li;
  }
  else{
    iRiS = lengths.length();
    iRiE = -1;
    for(int i1 = 0; i1 < lengths.length(); i1++){
      if(i1 == 0){
        if(li - 1 + lengths[i1] >= riS){
          iRiS = i1;
        }
      }
      else{
        if(li - 1 + lengths[i1 - 1] < riS && li - 1 + lengths[i1] >= riS){
          iRiS = i1;
        }
      }
      if(i1 == lengths.length() - 1){
        if(li - 1 + lengths[i1] <= riE){
          iRiE = i1;
        }
      }
      else{
        if(li - 1 + lengths[i1] <= riE && li - 1 + lengths[i1 + 1] > riE){
          iRiE = i1;
        }
      }
    }
  }
  return IntegerVector::create(iRiS, iRiE);
}



