#include "vlmcmethods.h"

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>

vlmcTree* generate_pst(List context_list, 
                       List probs,
                       unsigned int m){
  
  IntegerVector empty(0);
  List emptyList(m);
  unsigned int Hmax = 0;
  vlmcTree* tau = new vlmcTree(m, Hmax, empty, emptyList);
  
  for(int s=0; s<context_list.size(); s++){
    IntegerVector path = context_list[s];
    vlmcNode* node = tau->root;
    for(int d=0; d<path.size(); d++){
      if(node->children.size() == 0){
        node->growChildren(m, empty, emptyList);
      }
      node = node->children[path[d]];
    }
    node->vlmcLeaf = true;
    node->prob = as<NumericVector>(probs[s]);
  }
  if(tau->root->children.size() > 0) tau->root->vlmcLeaf = false;
  return(tau);
}

//' @export
// [[Rcpp::export]]
IntegerVector rvlmc_cpp(unsigned int n, 
                        List context_list, 
                        List probs){
  IntegerVector z(n);
  int m = as<IntegerVector>(probs[0]).size();
  IntegerVector A = seq_len(m) - 1;
  
  unsigned int H = 0;
  for(int s = 0; s<context_list.size(); s++){
    unsigned int suffix_len = as<IntegerVector>(context_list[s]).size();
    H = max(H, suffix_len); 
  }
  
  vlmcTree* tau = generate_pst(context_list, probs, m);
  
  for(unsigned int t=0; t<H; t++){
    z[t] = RcppArmadillo::sample(A, 1, false)[0];
  }
  for(unsigned int t=H; t<n; t++){
    vlmcNode* node = tau->root;
    unsigned int d = 1;
    while(node->children.size() > 0){
      node = node->children[z[t-d]];
      d++;
    }
    NumericVector p = node->prob;
    z[t] = RcppArmadillo::sample(A, 1, false, p)[0];
  }
  
  return(z);
}

// rvlmc_cpp(100, list(0L, c(1,0), c(1,1)), list(c(0.5,0.5), c(0.8, 0.2), c(0.5, 0.5)))
