#ifndef VLMCMETHODS_H
#define VLMCMETHODS_H

#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include <cmath>
#include <unordered_map>
#include <algorithm>
//#include <Rcpp.h>
#include <RcppArmadillo.h>

using namespace std;
using namespace Rcpp;

List contextAlgorithm(Rcpp::List data, unsigned int Hmax, unsigned int alphlen, double cutoff);

class vlmcNode{
public:
  vlmcNode(int lbl);
  ~vlmcNode();
  unsigned int label;
  unsigned int headLabel;
  vector<unsigned int> cnts_train;
  vector<unsigned int> cnts_test;
  double node_logq_train;
  double node_logq_test;
  double node_logq_diff;
  NumericVector prob;
  vector<vlmcNode*> children;
  vlmcNode* parent = NULL;
  vector<vlmcNode*> getNodes();
  vector<vlmcNode*> getLeaves();
  vector<vlmcNode*> getSiblings();
  unsigned int h;
  bool vlmcLeaf = false;
  bool tested = false;
  void growPerfect(unsigned int m, unsigned int H, IntegerVector renewal, LogicalMatrix allowedMatrix);
  string getPath();
  unsigned int getN_train();
  unsigned int getN_test();
  void growChildren(unsigned int m, IntegerVector renewal, LogicalMatrix allowedMatrix);
  bool renewal_limit = false;
  bool is_prohibited = false;
private:
  bool isLeaf();
  
};

class vlmcTree{
public:
  vlmcTree(unsigned int alphlen, unsigned int Hmax, IntegerVector renewal, LogicalMatrix allowedMatrix);
  ~vlmcTree();
  void addData_train(IntegerVector z, bool reset);
  void addData_test(IntegerVector z, bool reset);
  void cacheQ_train(double alpha, LogicalMatrix allowedMatrix);
  void cacheQ_test(double alpha, LogicalMatrix allowedMatrix);
  string concatLeaves();
  unsigned int n_train;
  unsigned int n_test;
  void clear();
  void pruneLeaf(vlmcNode* leaf);
  void growLeaf(vlmcNode* leaf);
  vector<vlmcNode*> getVlmcLeaves();
  vector<vlmcNode*> getPrunnableLeaves(bool is_c);
  vector<vlmcNode*> getGrowableLeaves();
  unsigned int H;
  unsigned int m;
  vlmcNode* root;
  LogicalMatrix allowedMatrix;
  void assignLimits(IntegerVector renewal);
};

#endif
