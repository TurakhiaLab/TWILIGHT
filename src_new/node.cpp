#ifndef PHYLO_HPP
#include "phylogeny.hpp"
#endif

#include <functional>

phylogeny::Node::Node(std::string id, float len){
    identifier = id;
    level = 1;
    branchLength = len;
    parent = nullptr;
    placed = false;
}

phylogeny::Node::Node(std::string id, Node* par, float len){
    identifier = id;
    branchLength = len;
    parent = par;
    level = par->level + 1;
    par->children.push_back(this);
    placed = false;
}

// Copy node
phylogeny::Node::Node(Node* node){
    this->identifier = node->identifier;
    this->branchLength = node->branchLength;
    this->level = node->level;
    this->weight = node->weight;
    this->numLeaves = node->numLeaves;
    this->grpID = node->grpID;
    this->partitionParent = node->partitionParent;
    this->seqsIncluded = node->seqsIncluded;
    this->alnLen = node->alnLen;
    this->alnWeight = node->alnWeight;
    this->alnNum = node->alnNum;
    // this->parent = node->parent;
    // this->children = node->children;
}

size_t phylogeny::Node::getNumLeaves(){
    size_t num_leaves = 0;
    if (children.size() == 0) return num_leaves;
    for (auto ch: children){
        if (ch->is_leaf()) num_leaves += 1;
        else num_leaves += ch->getNumLeaves();
    }
    return num_leaves;
}

size_t phylogeny::Node::getNumNodes(){
    size_t num_nodes = 1;
    if (children.size() == 0) return num_nodes;
    for (auto ch: children) num_nodes += ch->getNumNodes();
    return num_nodes;
}

void phylogeny::Node::collectPostOrder(std::stack<Node*>& postStack) {
    std::stack<Node*> s1;
    s1.push(this); 
    Node* current; 
    while (!s1.empty()) {
        current = s1.top(); 
        postStack.push(current);
        s1.pop(); 
        for (int i = current->children.size()-1; i >= 0; --i) {
            if (current->children[i]->grpID == current->grpID) s1.push(current->children[i]);     
        }
    } 
    return;
}