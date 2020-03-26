#include<RcppArmadillo.h>
#include <algorithm>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;
using namespace std;

NumericMatrix addInMatrix(NumericVector vec, NumericMatrix M);
arma::mat allComb(IntegerVector lengths);
int ancestor(int node, List tree);
List boundsCall(const Rcpp::NumericVector& x, const Rcpp::NumericVector& lengths, double alpha = 0.05, double q = NA_REAL, int fam = 0);
double bouPos(int li,int ri, int allInt, int liStart);
NumericVector colSums(const arma::mat & X);
List getNewCandidates(int subroot, int aN, List tree, IntegerVector minIVar, LogicalMatrix isOffMat, bool without = 0);
IntegerVector getOffspringTip(int node, List tree);
IntegerVector getOffspringTipB(int node, List tree);
double mlCost(NumericVector y, int fam);
arma::umat mycombn(double n, double k);
IntegerVector offspring(int node, List tree);
IntegerVector offspringAll(int node, List tree);
NumericVector offspringDiff(int Node, NumericVector node, List tree);
NumericVector offspringDiffTip(int Node, NumericVector node, List tree);
List prepTree(List& tree);
List sortComb(List& ncomb, NumericVector& cost);
double vecmin(NumericVector x);
double vecmax(NumericVector x);
List subsetSum(arma::umat numbers, IntegerVector target, int maxSize);
IntegerVector getiRi(int riS, int riE, int li, int allInt, const NumericVector& lengths);


//' @useDynLib treeSeg
//' @importFrom Rcpp evalCpp