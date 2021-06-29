#include "vlmcmethods.h"

//' @export
// [[Rcpp::export]]
void ibf(IntegerVector z_test,
         List z_train,
         IntegerVector renewal,
         double alpha = 1/2,
         unsigned int Hmax = 5,
         unsigned int alphlen = 2,
         unsigned int burnin = 10000){
  // Allocate maximal tree
  vlmcTree* tau = new vlmcTree(alphlen, Hmax, renewal);
  
  // Training: Metropolis-Hastings
  
  // Add data
  for(unsigned int i = 0; i<z_train.size(); i++){
    IntegerVector z = as<IntegerVector>(z_train[i]);
    bool is_first = i==0;
    tau->addData_train(z, is_first);
  }
  
  tau->addData_test(z_test, true);
  
  tau->cacheQ_train(alpha);
  tau->cacheQ_test(alpha);
  
  vector<vlmcNode*> nodes = tau->root->getNodes();
  for(int i=0; i<nodes.size(); i++){
    Rcout << nodes[i]->getPath() << "\n";
    Rcout << nodes[i]->node_logq_train << "\n";
    Rcout << nodes[i]->node_logq_test << "\n";
  }
  
  delete tau;
}

// ibf(binchain[[1]], binchain[2], renewal = 0)
