#ifndef PYYLO_HPP
#include "phylogeny.hpp"
#endif


// Helper functions for Partitioning the tree
size_t getNumLeaves(phylogeny::Node* node, int grpID) {
    size_t children = 0;
    if (node->grpID != grpID) return 0;
    if (node->children.size() == 0) return 1;
    for (auto &child: node->children)
        children += getNumLeaves(child, grpID);
    return children;
}

void updateCentroidEdge (phylogeny::Node* node, phylogeny::Node* root, size_t halfTaxa, size_t& Imbalance, phylogeny::Node*& breakEdge){
    if (node->grpID != root->grpID || node->children.size() == 0) return;
    for (auto ch: node->children){
        updateCentroidEdge(ch, root, halfTaxa, Imbalance, breakEdge);
    }
    size_t numDescendants = getNumLeaves(node, root->grpID);
    size_t tempImbalance = (halfTaxa > numDescendants) ? (halfTaxa - numDescendants) : (numDescendants - halfTaxa);
    if (tempImbalance < Imbalance) {
        breakEdge = node;
        Imbalance = tempImbalance;
    }
    return;
}

phylogeny::Node* getCentroidEdge(phylogeny::Node* node, int rootID){
    phylogeny::Node* centroidEdge = node;
    size_t numLeaves = getNumLeaves(node, node->grpID);
    size_t centroidImbalance = numLeaves;
    size_t halfTaxa = numLeaves/2;
    halfTaxa = (halfTaxa == 0) ? 1 : halfTaxa;
    updateCentroidEdge (node, node, halfTaxa, centroidImbalance, centroidEdge);
    return centroidEdge;
}

phylogeny::Node* getBreakingEdge(phylogeny::Node* root, int minSize){
    return getCentroidEdge(root, root->grpID);
}

void setChildrenGrpID(phylogeny::Node*& node, int ID_org, int ID) {
    if (node->grpID != ID_org) return;
    node->grpID = ID;
    if (node->children.size() == 0) return;
    for (auto ch: node->children) {
        setChildrenGrpID(ch, ID_org, ID);
    }
    return;
}

void phylogeny::PartitionInfo::bipartition(Node* root, Node* edge, Node*& tree1Root, Node*& tree2Root) {
    
    size_t tree1ID = (root->grpID == -1) ? 0 : root->grpID;
    size_t tree2ID = (root->grpID == -1) ? 1 : this->numPartitions + 1;
    this->numPartitions += 1;
    Node* head = edge->parent;
    int headID = edge->parent->grpID;
    while (true){
        if (head->parent == nullptr) break;
        if (head->parent->grpID != headID) break;
        head = head->parent;
    }
    tree1Root = head;
    tree2Root = edge;
    size_t tree2ID_org = tree2Root->grpID;
    size_t tree1ID_org = tree1Root->grpID;
    setChildrenGrpID(tree2Root, tree2ID_org, tree2ID);
    if (tree1Root->grpID == -1) {
        setChildrenGrpID(tree1Root, tree1ID_org, tree1ID);
    }
}

void phylogeny::PartitionInfo::partitionTree(Node* root) {
    size_t totalLeaves = getNumLeaves(root, root->grpID);
    if (totalLeaves <= this->maxPartitionSize) {
        if (this->partitionsRoot.empty()) {
            setChildrenGrpID(root, root->grpID, 0);
            size_t numLeaves = getNumLeaves(root, root->grpID);
            this->partitionsRoot[root->identifier] = std::make_pair(root, numLeaves);
        }
        return;
    }
    
    Node* breakEdge = getBreakingEdge(root, this->minPartitionSize);
    if (breakEdge->identifier == root->identifier) {
        return;
    }

    Node* tree1 = nullptr;
    Node* tree2 = nullptr;
    bipartition(root, breakEdge, tree1, tree2);
    
    size_t numTree1Leaves = getNumLeaves(tree1, tree1->grpID);
    size_t numTree2Leaves = getNumLeaves(tree2, tree2->grpID);
    if (root->parent == nullptr) {
        this->partitionsRoot[tree1->identifier] = std::make_pair(tree1, numTree2Leaves);
    }
    this->partitionsRoot[tree2->identifier] = std::make_pair(tree2, numTree2Leaves);
    this->partitionsRoot[tree1->identifier].second = numTree1Leaves;
    if (numTree2Leaves > this->maxPartitionSize) {
        partitionTree(tree2);
    }
    if (numTree1Leaves > this->maxPartitionSize) {
        partitionTree(tree1);
    }
    return;
}


phylogeny::PartitionInfo::~PartitionInfo() {
    this->partitionsRoot.clear();
}
