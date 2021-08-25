#include "vlmcmethods.h"

// [[Rcpp::export]]
List ibf(List z_test,
         IntegerVector z_train,
         IntegerVector renewal,
         List prohibited,
         double alpha = 1/2,
         double logprior_penalty = 2,
         unsigned int Hmax = 5,
         unsigned int alphlen = 2,
         unsigned int burnin = 10000,
         unsigned int nsamples = 100000){
  unsigned int all_samples = burnin + nsamples;
  int m = alphlen;
  // Allocate maximal tree
  Rcout << "Entrei ibf \n";
  vlmcTree* tau = new vlmcTree(alphlen, Hmax, renewal, prohibited);
  
  // Training: Metropolis-Hastings
  
  // Add data
  for(unsigned int i = 0; i<z_test.size(); i++){
    IntegerVector z = as<IntegerVector>(z_test[i]);
    bool is_first = i==0;
    tau->addData_test(z, is_first);
  }
  
  tau->addData_train(z_train, true);
  
  tau->cacheQ_train(alpha);
  tau->cacheQ_test(alpha);
 
  // Run Metropolis-Hastings on Training data.
  bool move_is_grow;
  vector<vlmcNode*> growableLeaves;
  unsigned int n_growable;
  vector<vlmcNode*> prunnableLeaves;
  unsigned int n_prunnable;
  double u, logratio;
  
  std::unordered_map<string, int> posterior;
  std::unordered_map<string, double> logQ;
  
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
        n_prunnable = tau->getPrunnableLeaves(false).size();
        logratio = -nodeToGrow->node_logq_diff + 
          log(n_prunnable) - log(m) - log(n_growable) -
          logprior_penalty;
        u = runif(1,0,1)[0];
        if(log(u)>logratio){ //Rejected proposed tree, so we roll back.
          tau->pruneLeaf(nodeToGrow->children[0]);
        }
      }
    } else {          // Propose prune
      prunnableLeaves = tau->getPrunnableLeaves(false);
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
    if(t > burnin){
      string this_tree = tau->concatLeaves();
      auto search = logQ.find(this_tree);
      if(search == logQ.end()){
        double logq_test = 0.0;
        for(const auto node : tau->getVlmcLeaves()){
          logq_test += node->node_logq_test;
        }
        logQ[this_tree] = logq_test;
      }
      posterior[this_tree]++;
    }
  }
  PutRNGstate();
  
  vector<string> key_posterior;
  vector<int> value_posterior;
  for(std::unordered_map<string,int>::iterator it = posterior.begin(); it != posterior.end(); ++it) {
    key_posterior.push_back(it->first);
    value_posterior.push_back(it->second);
  }
  vector<string> key_logQ;
  vector<double> value_logQ;
  for(std::unordered_map<string,double>::iterator it = logQ.begin(); it != logQ.end(); ++it) {
    key_logQ.push_back(it->first);
    value_logQ.push_back(it->second);
  }
  
  delete tau;
  
  return Rcpp::List::create(
    Rcpp::Named("posterior") = 
      Rcpp::List::create(Rcpp::Named("tree") = key_posterior, Rcpp::Named("count") = value_posterior),
    Rcpp::Named("logQ") = 
      Rcpp::List::create(Rcpp::Named("tree") = key_logQ, Rcpp::Named("logq") = value_logQ)
  );
}
