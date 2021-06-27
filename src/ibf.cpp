#include "vlmcmethods.h"

//' @export
// [[Rcpp::export]]
void ibf(IntegerVector renewal, unsigned int Hmax = 5,
                            unsigned int alphlen = 2){
  // Allocate maximal tree
  vlmcTree* tau = new vlmcTree(alphlen, Hmax, renewal);
  vector<vlmcNode*> leaves = tau->root->getNodes();
  for(int i = 0; i<leaves.size(); i++){
    Rcout << leaves[i]->getPath() << "\n";
  }
}
