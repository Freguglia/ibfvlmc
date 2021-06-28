#include "vlmcmethods.h"

//' @export
// [[Rcpp::export]]
void ibf(IntegerVector z, 
         IntegerVector renewal, 
         unsigned int Hmax = 5,
         unsigned int alphlen = 2,
         unsigned int burnin = 10000){
  // Allocate maximal tree
  vlmcTree* tau = new vlmcTree(alphlen, Hmax, renewal);
  tau->addData(z, true);
  
  vector<vlmcNode*> nodes = tau->root->getNodes();
  for(int i=0; i<nodes.size(); i++){
    Rcout << nodes[i]->getPath() << "\n";
    Rcout << nodes[i]->cnts[0] << " " << nodes[i]->cnts[1] << "\n";
  }
  
  delete tau;
}
