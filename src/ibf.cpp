#include "vlmcmethods.h"

//' @export
// [[Rcpp::export]]
void ibf(IntegerVector z_test,
         List z_train,
         IntegerVector renewal,
         double alpha = 1/2,
         double logprior_penalty = 2,
         unsigned int Hmax = 5,
         unsigned int alphlen = 2,
         unsigned int burnin = 10000,
         unsigned int nsamples = 100000){
  unsigned int all_samples = burnin + nsamples;
  int m = alphlen;
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
 
  // Run Metropolis-Hastings on Training data.
  bool move_is_grow;
  vector<vlmcNode*> growableLeaves;
  vlmcNode* nodeToGrow;
  unsigned int n_growable;
  vector<vlmcNode*> prunnableLeaves;
  vlmcNode* nodeToPrune;
  unsigned int n_prunnable;
  double u, logratio;
  
  GetRNGstate();
  for(unsigned int t=0; t<all_samples; t++){
    // Choose move
    move_is_grow = Rcpp::runif(1,0,1)[0] < 0.5;
    if(move_is_grow){ // Propose grow
      growableLeaves = tau->getGrowableLeaves();
      n_growable = growableLeaves.size();
      if(n_growable > 0){
        int toGrow = floor(runif(1,0,n_growable)[0]);
        vlmcNode* nodeToGrow = growableLeaves[toGrow];
        tau->growLeaf(nodeToGrow);
        n_prunnable = tau->getPrunnableLeaves().size();
        logratio = -nodeToGrow->node_logq_diff + 
          log(n_prunnable) - log(m) - log(n_growable) -
          logprior_penalty;
        u = runif(1,0,1)[0];
        if(log(u)>logratio){ //Rejected proposed tree, so we roll back.
          tau->pruneLeaf(nodeToGrow->children[0]);
        }
      }
    } else {          // Propose prune
      prunnableLeaves = tau->getPrunnableLeaves();
      n_prunnable = prunnableLeaves.size();
      if(n_prunnable > 0){
        int toPrune = floor(runif(1,0,n_prunnable)[0]);
        vlmcNode* nodeToPrune = prunnableLeaves[toPrune];
        tau->pruneLeaf(nodeToPrune);
        n_growable = tau->getGrowableLeaves().size();
        logratio = nodeToPrune->parent->node_logq_diff +
          log(n_growable) - log(n_prunnable) + log(m) +
          logprior_penalty;
        u = runif(1,0,1)[0];
        if(log(u)>logratio){ //Rejected proposed tree, so we roll back.
          tau->growLeaf(nodeToPrune->parent);
        }
      }
    }
  }
  PutRNGstate();
  
  Rcout << tau->concatLeaves() << "\n";
  
  //Rcout << tau->g;
  
  delete tau;
}

// ibf(binchain[[1]], binchain[2], renewal = 0)
