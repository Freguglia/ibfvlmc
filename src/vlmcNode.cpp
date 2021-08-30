#include "vlmcmethods.h"

vlmcNode::vlmcNode(int lbl){
  label = lbl;
}

vlmcNode::~vlmcNode(){
  if(this->children.size() > 0){
    for(unsigned int i = 0; i<this->children.size(); i++){
      delete this->children[i];
    }
  }
}

void vlmcNode::growChildren(unsigned int m, IntegerVector renewal, 
                            LogicalMatrix allowedMatrix){
  if(this->children.size() == 0 && !this->is_prohibited){
    vector<vlmcNode*> ch;
    unsigned int this_label = this->label;
    LogicalVector allowedPast(m);
    if(this->parent != NULL){
      allowedPast = allowedMatrix(_,this_label); 
    } else {
      for(unsigned int i=0; i<m; i++){
        allowedPast(i) = true;
      }
    }
    
    // Create m child nodes
    for(unsigned int i=0; i<m; i++){
      bool is_renewal = false;
      for(int j=0; j<renewal.size(); j++){
        if(renewal[j] == this_label){
          is_renewal = true;
        }
      }
      if(this->parent == NULL || !is_renewal){
        bool i_prohibited = !allowedPast(i);
        
        vlmcNode* a = new vlmcNode(i);
        a->parent = this;
        a->h = this->h + 1;
        if(a->h == 1){
          a->headLabel = i;
        } else {
          a->headLabel = this->headLabel;
        }
        if(i_prohibited){
          a->is_prohibited = true;
        } else {
          a->is_prohibited = false;
        }
        ch.push_back(a);
      }
    }
    this->children = ch;
  }
}

void vlmcNode::growPerfect(unsigned int m, unsigned int H, IntegerVector renewal,
                           LogicalMatrix allowedMatrix){
  // Only allowed in an empty tree.
  if(this->parent == NULL && this->children.size() == 0){
    unsigned int l;
    vector<vlmcNode*> currentLeaves;
    for(unsigned int h=0; h<H; h++){
      currentLeaves = this->getLeaves();
      l = currentLeaves.size();
      for(unsigned int i=0; i<l; i++){
        currentLeaves[i]->growChildren(m, renewal, allowedMatrix);
      }
    }
  }
}

bool vlmcNode::isLeaf(){
  return(this->children.size() == 0);
}

string vlmcNode::getPath(){
  vector<int> path;
  vlmcNode* node = this;
  path.push_back(node->label);
  while(node->parent != NULL){
    node = node->parent;
    path.push_back(node->label);
  }
  
  std::stringstream ss;
  for(size_t i = 0; i < path.size(); ++i){
    if(path[i] > -1){
      ss << path[i];
    } else {
      ss << 'r';
    }
  }
  std::string s = ss.str();
  s = std::string(s.rbegin(), s.rend());
  return(s);
}

vector<vlmcNode*> vlmcNode::getNodes(){
  vector<vlmcNode*> nodes;
  nodes.push_back(this);
  if(this->children.size() > 0){
    unsigned int i = 0;
    while(i < nodes.size()){
      for(unsigned int j=0; j<nodes[i]->children.size(); j++){
        nodes.push_back(nodes[i]->children[j]);
      }    
      i++;
    }
  }
  return(nodes);
}

vector<vlmcNode*> vlmcNode::getLeaves(){
  vector<vlmcNode*> leaves;
  vector<vlmcNode*> allNodes = this->getNodes();
  unsigned int l = allNodes.size();
  for(unsigned int i=0; i<l; i++){
    if(allNodes[i]->isLeaf()){
      leaves.push_back(allNodes[i]);
    }
  }
  return(leaves);
}

vector<vlmcNode*> vlmcNode::getSiblings(){
  if(this->parent == NULL){
    vector<vlmcNode*> rootVec;
    rootVec.push_back(this);
    return(rootVec);
  } else {
    return(this->parent->children);
  }
}

unsigned int vlmcNode::getN_train(){
  unsigned int n_tot = 0;
  for(unsigned int n : this->cnts_train){n_tot += n;}
  return(n_tot);
}

unsigned int vlmcNode::getN_test(){
  unsigned int n_tot = 0;
  for(unsigned int n : this->cnts_train){n_tot += n;}
  return(n_tot);
}
