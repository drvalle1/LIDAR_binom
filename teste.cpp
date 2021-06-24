// [[Rcpp::depends("RcppArmadillo")]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <iostream>
#include <ctime>
#include <fstream>
using namespace Rcpp;

/***************************************************************************************************************************/
/*********************************                      UTILS          *****************************************************/
/***************************************************************************************************************************/

// This function makes multinomial draws
//// [[Rcpp::export]]
// IntegerVector rmultinom_1(int size, NumericVector probs, int N) {
//   IntegerVector outcome(N);
//   rmultinom(size, probs.begin(), N, outcome.begin());
//   return outcome;
// }

// This function makes multinomial draws
//Got this from https://stackoverflow.com/questions/24618370/using-rmultinom-with-rcpp
// [[Rcpp::export]]
IntegerVector rmultinom_1(NumericVector probs, int size) {
  int k = probs.size();
  IntegerVector ans(k);
  rmultinom(size, probs.begin(), k, ans.begin());
  return(ans);
}