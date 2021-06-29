#ifndef VLMCMETHODS_H
#define VLMCMETHODS_H

#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include <cmath>
#include <unordered_map>
#include <algorithm>
#include <Rcpp.h>

using namespace std;
using namespace Rcpp;

List contextAlgorithm(Rcpp::List data, unsigned int Hmax, unsigned int alphlen, double cutoff);

class vlmcNode{
public:
  vlmcNode(int lbl);
  ~vlmcNode();
  unsigned int label;
  vector<unsigned int> cnts_train;
  vector<unsigned int> cnts_test;
  double node_logq_train;
  double node_logq_test;
  vector<vlmcNode*> children;
  vlmcNode* parent = NULL;
  vector<vlmcNode*> getNodes();
  vector<vlmcNode*> getLeaves();
  vector<vlmcNode*> getSiblings();
  unsigned int h;
  bool vlmcLeaf = false;
  bool tested = false;
  void growPerfect(unsigned int m, unsigned int H, IntegerVector renewal);
  string getPath();
  unsigned int getN_train();
  unsigned int getN_test();
private:
  bool isLeaf();
  void growChildren(unsigned int m, IntegerVector renewal);
};

class vlmcTree{
public:
  vlmcTree(unsigned int alphlen, unsigned int Hmax, IntegerVector renewal);
  ~vlmcTree();
  void addData_train(IntegerVector z, bool reset);
  void addData_test(IntegerVector z, bool reset);
  void cacheQ_train(double alpha);
  string concatLeaves();
  unsigned int n_train;
  unsigned int n_test;
  void clear();
  void pruneLeaf(vlmcNode* leaf);
  vector<vlmcNode*> getVlmcLeaves();
  vector<vlmcNode*> getUncheckedLeaves();
  vector<vlmcNode*> getPrunnableLeaves();
  unsigned int H;
  unsigned int m;
  vlmcNode* root;
};

#endif