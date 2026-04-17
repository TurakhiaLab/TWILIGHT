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

// --- Centroid Decomposition (Fallback Strategy) ---
void updateCentroidEdge (phylogeny::Node* node, phylogeny::Node* root, size_t totalLeaves, size_t& minMaxPartition, phylogeny::Node*& breakEdge){
    if (node->grpID != root->grpID) return;

    for (auto ch: node->children){
        updateCentroidEdge(ch, root, totalLeaves, minMaxPartition, breakEdge);
    }

    if (node != root) {
        size_t numDescendants = getNumLeaves(node, root->grpID);
        if (numDescendants > 0 && numDescendants < totalLeaves) {
            size_t maxPartition = std::max(numDescendants, totalLeaves - numDescendants);
            if (maxPartition < minMaxPartition) {
                breakEdge = node;
                minMaxPartition = maxPartition;
            }
        }
    }
}

phylogeny::Node* getCentroidEdge(phylogeny::Node* node, int rootID){
    size_t numLeaves = getNumLeaves(node, node->grpID);
    if (numLeaves <= 1) return node;

    phylogeny::Node* centroidEdge = nullptr;
    for(auto child : node->children){
        if(child->grpID == node->grpID){
            centroidEdge = child;
            break;
        }
    }
    if (!centroidEdge) return node;

    size_t minMaxPartition = numLeaves;
    updateCentroidEdge (node, node, numLeaves, minMaxPartition, centroidEdge);
    return centroidEdge;
}

// --- Greedy Strategy ---
void findGreedyEdge(phylogeny::Node* node, phylogeny::Node* root, size_t maxPartitionSize, int& bestSize, phylogeny::Node*& breakEdge) {
    if (node->grpID != root->grpID) return;

    for (auto ch : node->children) {
        findGreedyEdge(ch, root, maxPartitionSize, bestSize, breakEdge);
    }

    if (node != root) {
        size_t numDescendants = getNumLeaves(node, root->grpID);
        if (numDescendants <= maxPartitionSize && (int)numDescendants > bestSize) {
            bestSize = numDescendants;
            breakEdge = node;
        }
    }
}

phylogeny::Node* getGreedyEdge(phylogeny::Node* node, size_t maxPartitionSize) {
    phylogeny::Node* breakEdge = nullptr;
    int bestSize = -1; 
    findGreedyEdge(node, node, maxPartitionSize, bestSize, breakEdge);
    return breakEdge;
}

void setChildrenGrpID(phylogeny::Node* node, int ID_org, int ID) {
    if (node->grpID != ID_org) return;
    node->grpID = ID;
    if (node->children.size() == 0) return;
    for (auto ch: node->children) {
        setChildrenGrpID(ch, ID_org, ID);
    }
    return;
}

void phylogeny::PartitionInfo::bipartition(phylogeny::Node* rootOfPartition, phylogeny::Node* edgeToCut, int newGrpID) {
    setChildrenGrpID(edgeToCut, rootOfPartition->grpID, newGrpID);
}

void phylogeny::PartitionInfo::partitionTree(Node* root) {
    std::stack<phylogeny::Node*> treesToPartition;
    
    setChildrenGrpID(root, -1, 0);
    treesToPartition.push(root);
    this->numPartitions = 0;

    while (!treesToPartition.empty()) {
        phylogeny::Node* currentRoot = treesToPartition.top();
        treesToPartition.pop();

        size_t totalLeaves = getNumLeaves(currentRoot, currentRoot->grpID);

        if (totalLeaves <= this->maxPartitionSize) {
            this->partitionsRoot[currentRoot->identifier] = std::make_pair(currentRoot, totalLeaves);
            continue;
        }

        Node* breakEdge = getGreedyEdge(currentRoot, this->maxPartitionSize);

        if (breakEdge == nullptr) {
            breakEdge = getCentroidEdge(currentRoot, currentRoot->grpID);
        }
        
        if (breakEdge == nullptr || breakEdge == currentRoot) {
            this->partitionsRoot[currentRoot->identifier] = std::make_pair(currentRoot, totalLeaves);
            continue;
        }

        this->numPartitions++;
        bipartition(currentRoot, breakEdge, this->numPartitions);

        this->partitionsRoot[breakEdge->identifier] = std::make_pair(breakEdge, getNumLeaves(breakEdge, breakEdge->grpID));
        treesToPartition.push(currentRoot);
    }
    return;
}


phylogeny::PartitionInfo::~PartitionInfo() {
    this->partitionsRoot.clear();
}
