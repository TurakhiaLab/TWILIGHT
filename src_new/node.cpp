#ifndef PHYLO_HPP
#include "phylogeny.hpp"
#endif

phylogeny::Node::Node(std::string id, float len){
    identifier = id;
    level = 1;
    branchLength = len;
    parent = nullptr;
}

phylogeny::Node::Node(std::string id, Node* par, float len){
    identifier = id;
    branchLength = len;
    parent = par;
    level = par->level + 1;
    par->children.push_back(this);
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