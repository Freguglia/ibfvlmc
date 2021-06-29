#include "vlmcmethods.h"

void vlmcTree::clear(){
  vector<vlmcNode*> allNodes = this->root->getNodes();
  for(unsigned int i=0; i<allNodes.size(); i++){
    delete allNodes[i];
  }
}

vlmcTree::vlmcTree(unsigned int alphlen, unsigned int Hmax, IntegerVector renewal){
  root = new vlmcNode(-1);
  root->h = 0;
  root->growPerfect(alphlen, Hmax, renewal);
  H = Hmax;
  m = alphlen;
  n_train = 0;
  n_test = 0;
  vector<vlmcNode*> nodes = root->getLeaves();
  for(unsigned int i = 0; i<nodes.size(); i++){
    nodes[i]->vlmcLeaf = true;
  }
}

vlmcTree::~vlmcTree(){
  delete root;
}

void vlmcTree::addData_train(IntegerVector z, bool reset = true){
  unsigned int H = this->H;
  unsigned int m = this->m;
  unsigned int n = z.size();
  vector<vlmcNode*> nodes = this->root->getNodes();
  vlmcNode* root = this->root;
  if(reset){
    for(unsigned int i = 0; i<nodes.size(); i++){
      for(unsigned int j = 0; j<m; j++){
        nodes[i]->cnts_train.push_back(0);
      }
    }   
  }
  vlmcNode* node;
  for(unsigned int t=H; t<n; t++){
    node = root;
    node->cnts_train[z[t]]++;
    unsigned int d = 0;
    while(node->children.size() > 0){
      d++;
      node = node->children[z[t-d]];
      node->cnts_train[z[t]]++;
    }
  }
  this->n_train += (n - H);
}

void vlmcTree::addData_test(IntegerVector z, bool reset = true){
  unsigned int H = this->H;
  unsigned int m = this->m;
  unsigned int n = z.size();
  vector<vlmcNode*> nodes = this->root->getNodes();
  vlmcNode* root = this->root;
  if(reset){
    for(unsigned int i = 0; i<nodes.size(); i++){
      for(unsigned int j = 0; j<m; j++){
        nodes[i]->cnts_test.push_back(0);
      }
    }   
  }
  vlmcNode* node;
  for(unsigned int t=H; t<n; t++){
    node = root;
    node->cnts_train[z[t]]++;
    unsigned int d = 0;
    while(node->children.size() > 0){
      d++;
      node = node->children[z[t-d]];
      node->cnts_test[z[t]]++;
    }
  }
  this->n_test += (n - H);
}

void vlmcTree::cacheQ_train(double alpha){
  vector<vlmcNode*> nodes = this->root->getNodes();
  unsigned int m = this->m;
  double logq;
  unsigned int n_tot;
  vector<double> ns(m);
  for(unsigned int i=0; i<nodes.size(); i++){
    logq = 0;
    n_tot = 0;
    for (unsigned int j=0; j<m; j++){
      ns[j] = (double)nodes[i]->cnts_train[j];
      n_tot += ns[j];
    }
    for (unsigned int j=0; j<m; j++){
      logq += lgamma(ns[j] + alpha);
    }
    logq = logq - lgamma(n_tot + m*alpha);
    nodes[i]->node_logq_train = logq;
  }
}

vector<vlmcNode*> vlmcTree::getVlmcLeaves(){
  vector<vlmcNode*> allNodes = this->root->getNodes();
  vector<vlmcNode*> tau;
  for(unsigned int i=0; i<allNodes.size(); i++){
    if(allNodes[i]->vlmcLeaf){
      tau.push_back(allNodes[i]);
    }
  }
  return(tau);
}

vector<vlmcNode*> vlmcTree::getPrunnableLeaves(){
  vector<vlmcNode*> currentLeaves = this->getVlmcLeaves();
  vector<vlmcNode*> res;
  vector<vlmcNode*> sibs;
  int m = this->m;
  int sibCount;
  for(long unsigned int i=0; i<currentLeaves.size(); i++){
    sibs = currentLeaves[i]->getSiblings();
    sibCount = 0;
    for(long unsigned int j = 0; j<sibs.size(); j++){
      if(sibs[j]->vlmcLeaf) sibCount++;
    }
    if(sibCount == m) res.push_back(currentLeaves[i]);
  }
  return(res);
}

vector<vlmcNode*> vlmcTree::getUncheckedLeaves(){
  vector<vlmcNode*> allLeaves = this->getPrunnableLeaves();
  vector<vlmcNode*> tau;
  for(unsigned int i=0; i<allLeaves.size(); i++){
    if(!allLeaves[i]->tested){
      tau.push_back(allLeaves[i]);
    }
  }
  return(tau);
}

void vlmcTree::pruneLeaf(vlmcNode* leaf){
  vector<vlmcNode*> sibs = leaf->getSiblings();
  if(leaf->parent == NULL){
    throw invalid_argument("Error: cannot prune root node.");
  }
  for(vlmcNode* sib : sibs){
    if(!sib->vlmcLeaf){
      throw invalid_argument("Error: Tried to prune a leaf which itself or one of its siblings is not a leaf.");
    }
    sib->vlmcLeaf = false;
  }
  leaf->parent->vlmcLeaf = true;
}

string vlmcTree::concatLeaves(){
  vector<vlmcNode*> leaves = this->getVlmcLeaves();
  vector<string> paths;
  for(long unsigned int i=0; i<leaves.size(); i++){
    paths.push_back(leaves[i]->getPath());
  }
  string s = "{";
  for(const auto &piece : paths) s += piece + ',';
  s.pop_back();
  return s + '}';
}